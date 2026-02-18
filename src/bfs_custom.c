// bfs_custom.c - Version Hybride Corrigée (Top-Down / Bottom-Up)
// Beamer et al. direction-optimizing BFS pour Graph500
#include "common.h"
#include "csr_reference.h"
#include "bitmap_reference.h"
#include <stdint.h>
#include <stdlib.h>
#include <string.h>
#include <mpi.h>

#undef ulong_mask
#undef ulong_shift
#define ulong_bits 64
#define ulong_mask 63
#define ulong_shift 6

#define SET_BIT(bm, v)  do { (bm)[(v) >> 6] |= (1ULL << ((v) & 63)); } while (0)
#define TEST_BIT(bm, v) (((bm)[(v) >> 6] & (1ULL << ((v) & 63))) != 0)

/*
 * Bitmap global : le vertex global gv est mappé à l'offset
 *   offset = VERTEX_OWNER(gv) * max_local_verts + VERTEX_LOCAL(gv)
 */
#define GLOBAL_OFFSET(gv) \
    (VERTEX_OWNER(gv) * max_local_verts + VERTEX_LOCAL(gv))

#define TEST_BIT_GLOBAL(bm, gv) \
    (((bm)[GLOBAL_OFFSET(gv) >> 6] & (1ULL << (GLOBAL_OFFSET(gv) & 63))) != 0)

#define SET_BIT_GLOBAL(bm, gv) do { \
    int64_t _off = GLOBAL_OFFSET(gv); \
    (bm)[_off >> 6] |= (1ULL << (_off & 63)); \
} while(0)

/* Paramètres de switching Beamer et al. */
#define ALPHA 15.0
#define BETA  24.0

/* ------------------------------------------------------------------ */
/* Variables globales (non-static : attendues par csr_reference.c,    */
/* validate.c via extern)                                              */
/* ------------------------------------------------------------------ */
unsigned long  *visited;
unsigned long  *frontier_bitmap;
unsigned long  *next_frontier_bitmap;
int64_t         visited_size;

int64_t        *pred_glob;
oned_csr_graph  g;         /* attendu par csr_reference.c */
int            *rowstarts; /* attendu par validate.c      */
extern void           *column;    /* attendu par COLUMN() macro  */

int     mpi_rank, mpi_size;
int64_t max_local_verts = 0;

/* bitmap global alloué une seule fois */
static unsigned long *global_frontier    = NULL;
static int64_t        global_bitmap_size = 0;

/* ------------------------------------------------------------------ */
/* Initialisation                                                       */
/* ------------------------------------------------------------------ */

void make_graph_data_structure(const tuple_graph* const tg) {
    MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &mpi_size);

    convert_graph_to_oned_csr(tg, &g);

    /* Exposer pour validate.c et la macro COLUMN() */
    column    = g.column;
    rowstarts = g.rowstarts;

    visited_size = (g.nlocalverts + ulong_bits - 1) / ulong_bits;

    visited              = xmalloc(visited_size * sizeof(unsigned long));
    frontier_bitmap      = xmalloc(visited_size * sizeof(unsigned long));
    next_frontier_bitmap = xmalloc(visited_size * sizeof(unsigned long));

    /* Calculer max_local_verts sur tous les PEs */
    max_local_verts = g.nlocalverts;
    MPI_Allreduce(MPI_IN_PLACE, &max_local_verts, 1,
                  MPI_INT64_T, MPI_MAX, MPI_COMM_WORLD);

    /* Allouer le bitmap global une seule fois */
    global_bitmap_size = ((max_local_verts * (int64_t)mpi_size) + ulong_bits - 1) / ulong_bits;
    global_frontier    = xmalloc(global_bitmap_size * sizeof(unsigned long));
}

/* ------------------------------------------------------------------ */
/* Utilitaires                                                          */
/* ------------------------------------------------------------------ */

static inline int64_t count_set_bits(unsigned long *bitmap, int64_t size) {
    int64_t count = 0;
    for (int64_t i = 0; i < size; i++)
        count += __builtin_popcountll(bitmap[i]);
    return count;
}

/*
 * Construit et synchronise la frontière globale via Allreduce(MPI_BOR).
 */
static void build_global_frontier(void) {
    memset(global_frontier, 0, global_bitmap_size * sizeof(unsigned long));
    for (int64_t v = 0; v < g.nlocalverts; v++) {
        if (TEST_BIT(frontier_bitmap, v))
            SET_BIT_GLOBAL(global_frontier, VERTEX_TO_GLOBAL(mpi_rank, v));
    }
    MPI_Allreduce(MPI_IN_PLACE, global_frontier, global_bitmap_size,
                  MPI_UNSIGNED_LONG, MPI_BOR, MPI_COMM_WORLD);
}




/* ------------------------------------------------------------------ */
/* TOP-DOWN                                                             */
/* ------------------------------------------------------------------ */

static void bfs_step_top_down(void) {
    // 1. Synchronisation unique de la frontière
    build_global_frontier();

    // 2. Une seule passe sur les sommets locaux non-visités
    // On cherche si l'un de leurs voisins appartient à la frontière globale
    for (int64_t v = 0; v < g.nlocalverts; v++) {
        if (TEST_BIT(visited, v)) continue;

        int64_t row_start = g.rowstarts[v];
        int64_t row_end   = g.rowstarts[v + 1];

        for (int64_t j = row_start; j < row_end; j++) {
            int64_t nb_global = COLUMN(j);

            // Si le voisin est dans la frontière (qu'il soit local ou distant)
            if (TEST_BIT_GLOBAL(global_frontier, nb_global)) {
                SET_BIT(visited, v);
                SET_BIT(next_frontier_bitmap, v);
                pred_glob[v] = nb_global; // Parent trouvé
                break; // On passe au sommet local suivant
            }
        }
    }
}
/* ------------------------------------------------------------------ */
/* BOTTOM-UP                                                            */
/* ------------------------------------------------------------------ */
/*
 * Chaque vertex local non-visité cherche UN parent dans la frontière
 * globale. Le break dès le premier parent trouvé est l'optimisation
 * clé de Beamer et al.
 */
static void bfs_step_bottom_up(void) {
    build_global_frontier();

    for (int64_t v = 0; v < g.nlocalverts; v++) {
        if (TEST_BIT(visited, v)) continue;

        for (int64_t j = g.rowstarts[v]; j < g.rowstarts[v + 1]; j++) {
            int64_t nb_global = COLUMN(j);

            if (TEST_BIT_GLOBAL(global_frontier, nb_global)) {
                SET_BIT(visited, v);
                SET_BIT(next_frontier_bitmap, v);
                pred_glob[v] = nb_global;
                break;  /* Un seul parent suffit */
            }
        }
    }
}

/* ------------------------------------------------------------------ */
/* BFS principal                                                        */
/* ------------------------------------------------------------------ */

void run_bfs(int64_t root, int64_t* pred) {
    pred_glob = pred;
    clean_pred(pred); // Utilisation de votre fonction existante

    memset(visited,              0, visited_size * sizeof(unsigned long));
    memset(frontier_bitmap,      0, visited_size * sizeof(unsigned long));
    memset(next_frontier_bitmap, 0, visited_size * sizeof(unsigned long));

    int     root_pe    = VERTEX_OWNER(root);
    int64_t root_local = VERTEX_LOCAL(root);

    if (root_pe == mpi_rank) {
        pred[root_local] = root;
        SET_BIT(visited,         root_local);
        SET_BIT(frontier_bitmap, root_local);
    }

    /* Calculer les totaux une seule fois au début */
    int64_t total_vertices = 0;
    int64_t local_verts = g.nlocalverts;
    MPI_Allreduce(&local_verts, &total_vertices, 1, MPI_INT64_T, MPI_SUM, MPI_COMM_WORLD);

    int64_t total_edges = 0;
    for (int64_t i = 0; i < g.nlocalverts; i++)
        total_edges += g.rowstarts[i + 1] - g.rowstarts[i];
    MPI_Allreduce(&total_edges, &total_edges, 1, MPI_INT64_T, MPI_SUM, MPI_COMM_WORLD);

    int64_t mu = total_edges;
    int use_bottom_up = 0;

    for (int iteration = 0; iteration < 10000; iteration++) {
        // --- OPTIMISATION : UN SEUL ALLREDUCE POUR MF ET FRONTIER_COUNT ---
        int64_t local_data[2] = {0, 0}; // [0] = count, [1] = mf
        for (int64_t v = 0; v < g.nlocalverts; v++) {
            if (TEST_BIT(frontier_bitmap, v)) {
                local_data[0]++;
                local_data[1] += (g.rowstarts[v + 1] - g.rowstarts[v]);
            }
        }
        
        int64_t global_data[2];
        MPI_Allreduce(local_data, global_data, 2, MPI_INT64_T, MPI_SUM, MPI_COMM_WORLD);
        
        int64_t frontier_count = global_data[0];
        int64_t mf             = global_data[1];

        if (frontier_count == 0) break;

        // --- DÉCISION DE SWITCHING (BEAMER ET AL.) ---
        if (!use_bottom_up && (double)mf > (double)mu / ALPHA) {
            use_bottom_up = 1;
        } else if (use_bottom_up && (double)frontier_count < (double)total_vertices / BETA) {
            use_bottom_up = 0;
        }

        // --- EXÉCUTION ---
        if (use_bottom_up) {
            bfs_step_bottom_up();
        } else {
            bfs_step_top_down();
        }

        mu -= mf;
        if (mu < 0) mu = 0;

        // Swap efficace
        unsigned long *tmp = frontier_bitmap;
        frontier_bitmap = next_frontier_bitmap;
        next_frontier_bitmap = tmp;
        memset(next_frontier_bitmap, 0, visited_size * sizeof(unsigned long));
    }
}
/* ------------------------------------------------------------------ */
/* Fonctions requises par Graph500                                      */
/* ------------------------------------------------------------------ */

void get_edge_count_for_teps(int64_t* edge_visit_count) {
    int64_t edge_count = 0;
    for (int64_t i = 0; i < g.nlocalverts; i++) {
        if (pred_glob[i] != -1) {
            for (int64_t j = g.rowstarts[i]; j < g.rowstarts[i + 1]; j++) {
                if (COLUMN(j) <= VERTEX_TO_GLOBAL(mpi_rank, i))
                    edge_count++;
            }
        }
    }
    MPI_Allreduce(MPI_IN_PLACE, &edge_count, 1,
                  MPI_INT64_T, MPI_SUM, MPI_COMM_WORLD);
    *edge_visit_count = edge_count;
}

void clean_pred(int64_t* pred) {
    for (int64_t i = 0; i < g.nlocalverts; i++)
        pred[i] = -1;
}

void free_graph_data_structure(void) {
    free_oned_csr_graph(&g);
    free(visited);
    free(frontier_bitmap);
    free(next_frontier_bitmap);
    free(global_frontier);
}

size_t get_nlocalverts_for_pred(void) {
    return g.nlocalverts;
}