// bfs_custom.c - Version Hybride Corrigée (Top-Down / Bottom-Up)
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

#define SET_BIT(bm, v) do { bm[(v) >> 6] |= (1ULL << ((v) & 63)); } while (0)
#define TEST_BIT(bm, v) ((bm[(v) >> 6] & (1ULL << ((v) & 63))) != 0)

// Macro pour tester un bit dans un bitmap GLOBAL (indexé par vertex global)
#define TEST_BIT_GLOBAL(bm, global_v) ({ \
    int64_t owner = VERTEX_OWNER(global_v); \
    int64_t local_v = VERTEX_LOCAL(global_v); \
    int64_t offset = owner * max_local_verts + local_v; \
    ((bm[(offset) >> 6] & (1ULL << ((offset) & 63))) != 0); \
})

#define SET_BIT_GLOBAL(bm, global_v) do { \
    int64_t owner = VERTEX_OWNER(global_v); \
    int64_t local_v = VERTEX_LOCAL(global_v); \
    int64_t offset = owner * max_local_verts + local_v; \
    bm[(offset) >> 6] |= (1ULL << ((offset) & 63)); \
} while(0)

// Paramètres de switching
#define ALPHA 15.0
#define BETA 24.0

// Structures globales
unsigned long *visited;
unsigned long *frontier_bitmap;
unsigned long *next_frontier_bitmap;
int64_t visited_size;

int64_t *pred_glob, *column;
int *rowstarts;
oned_csr_graph g;

int mpi_rank, mpi_size;
int64_t max_local_verts = 0;  // Maximum de vertices par PE

void make_graph_data_structure(const tuple_graph* const tg) {
    MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &mpi_size);
    
    convert_graph_to_oned_csr(tg, &g);
    
    column = g.column;
    rowstarts = g.rowstarts;
    visited_size = (g.nlocalverts + ulong_bits - 1) / ulong_bits;
    
    visited = xmalloc(visited_size * sizeof(unsigned long));
    frontier_bitmap = xmalloc(visited_size * sizeof(unsigned long));
    next_frontier_bitmap = xmalloc(visited_size * sizeof(unsigned long));
    
    // Calculer le max de vertices locaux sur tous les PEs
    max_local_verts = g.nlocalverts;
    MPI_Allreduce(MPI_IN_PLACE, &max_local_verts, 1, MPI_INT64_T, MPI_MAX, MPI_COMM_WORLD);
}

static inline int64_t count_set_bits(unsigned long *bitmap, int64_t size) {
    int64_t count = 0;
    for (int64_t i = 0; i < size; i++) {
        count += __builtin_popcountll(bitmap[i]);
    }
    return count;
}

/**
 * TOP-DOWN : Explorer depuis la frontière vers les voisins
 * On synchronise la frontière globale, puis chaque PE explore ses vertices locaux
 */
static void bfs_step_top_down(void) {
    // Taille du bitmap global = nombre total de vertices possibles
    int64_t global_bitmap_size = ((max_local_verts * mpi_size) + ulong_bits - 1) / ulong_bits;
    unsigned long *global_frontier = calloc(global_bitmap_size, sizeof(unsigned long));
    
    // Chaque PE copie sa frontière locale dans le bitmap global
    for (int64_t v = 0; v < g.nlocalverts; v++) {
        if (TEST_BIT(frontier_bitmap, v)) {
            int64_t global_v = VERTEX_TO_GLOBAL(mpi_rank, v);
            SET_BIT_GLOBAL(global_frontier, global_v);
        }
    }
    
    // Synchroniser pour avoir la frontière complète
    MPI_Allreduce(MPI_IN_PLACE, global_frontier, global_bitmap_size, 
                  MPI_UNSIGNED_LONG, MPI_BOR, MPI_COMM_WORLD);
    
    // Chaque vertex local explore ses voisins
    for (int64_t v = 0; v < g.nlocalverts; v++) {
        if (TEST_BIT(visited, v)) continue;
        
        int64_t row_start = g.rowstarts[v];
        int64_t row_end = g.rowstarts[v + 1];
        
        for (int64_t j = row_start; j < row_end; j++) {
            int64_t neighbor_global = COLUMN(j);
            
            // Vérifier si ce voisin est dans la frontière globale
            if (TEST_BIT_GLOBAL(global_frontier, neighbor_global)) {
                SET_BIT(visited, v);
                SET_BIT(next_frontier_bitmap, v);
                pred_glob[v] = neighbor_global;
                break;
            }
        }
    }
    
    free(global_frontier);
}

/**
 * BOTTOM-UP : Chaque vertex non-visité cherche un parent dans la frontière
 * Plus efficace quand la frontière est grande (beaucoup de vertices actifs)
 */
static void bfs_step_bottom_up(void) {
    // Synchroniser la frontière globale
    int64_t global_bitmap_size = ((max_local_verts * mpi_size) + ulong_bits - 1) / ulong_bits;
    unsigned long *global_frontier = calloc(global_bitmap_size, sizeof(unsigned long));
    
    for (int64_t v = 0; v < g.nlocalverts; v++) {
        if (TEST_BIT(frontier_bitmap, v)) {
            int64_t global_v = VERTEX_TO_GLOBAL(mpi_rank, v);
            SET_BIT_GLOBAL(global_frontier, global_v);
        }
    }
    
    MPI_Allreduce(MPI_IN_PLACE, global_frontier, global_bitmap_size, 
                  MPI_UNSIGNED_LONG, MPI_BOR, MPI_COMM_WORLD);
    
    // Pour chaque vertex local non-visité
    for (int64_t v = 0; v < g.nlocalverts; v++) {
        if (TEST_BIT(visited, v)) continue;
        
        int64_t row_start = g.rowstarts[v];
        int64_t row_end = g.rowstarts[v + 1];
        
        // Chercher UN parent dans la frontière
        for (int64_t j = row_start; j < row_end; j++) {
            int64_t neighbor_global = COLUMN(j);
            
            if (TEST_BIT_GLOBAL(global_frontier, neighbor_global)) {
                SET_BIT(visited, v);
                SET_BIT(next_frontier_bitmap, v);
                pred_glob[v] = neighbor_global;
                break;  // Parent trouvé !
            }
        }
    }
    
    free(global_frontier);
}

void run_bfs(int64_t root, int64_t* pred) {
    pred_glob = pred;
    
    // Initialisation
    memset(visited, 0, visited_size * sizeof(unsigned long));
    memset(frontier_bitmap, 0, visited_size * sizeof(unsigned long));
    memset(next_frontier_bitmap, 0, visited_size * sizeof(unsigned long));
    
    // Setup root
    int root_pe = VERTEX_OWNER(root);
    int64_t root_local = VERTEX_LOCAL(root);
    
    if (root_pe == mpi_rank) {
        pred[root_local] = root;
        SET_BIT(visited, root_local);
        SET_BIT(frontier_bitmap, root_local);
    }
    
    MPI_Barrier(MPI_COMM_WORLD);
    
    // Calculer le nombre total de vertices
    int64_t total_vertices = g.nlocalverts;
    MPI_Allreduce(MPI_IN_PLACE, &total_vertices, 1, MPI_INT64_T, MPI_SUM, MPI_COMM_WORLD);
    
    // Calculer les edges totaux pour le critère de switching
    int64_t total_edges = 0;
    for (int64_t i = 0; i < g.nlocalverts; i++) {
        total_edges += g.rowstarts[i + 1] - g.rowstarts[i];
    }
    MPI_Allreduce(MPI_IN_PLACE, &total_edges, 1, MPI_INT64_T, MPI_SUM, MPI_COMM_WORLD);
    
    int use_bottom_up = 0;
    int iteration = 0;
    const int MAX_ITER = 1000;
    
    while (iteration < MAX_ITER) {
        iteration++;
        
        // Compter la frontière
        int64_t local_frontier = count_set_bits(frontier_bitmap, visited_size);
        int64_t frontier_count = local_frontier;
        MPI_Allreduce(MPI_IN_PLACE, &frontier_count, 1, MPI_INT64_T, MPI_SUM, MPI_COMM_WORLD);
        
        if (frontier_count == 0) break;
        
        // DÉCISION DE SWITCHING basée sur heuristiques de Beamer et al.
        if (!use_bottom_up) {
            // Estimer edges à explorer
            int64_t scout_count = frontier_count * 16;  // Approximation (edgefactor=16)
            int64_t edges_remaining = total_edges - (iteration * frontier_count * 16);
            
            if (scout_count > edges_remaining / ALPHA) {
                use_bottom_up = 1;
            }
        } else {
            // Revenir à top-down si frontière devient petite
            if (frontier_count < total_vertices / BETA) {
                use_bottom_up = 0;
            }
        }
        
        // Exécuter l'approche choisie
        if (use_bottom_up) {
            bfs_step_bottom_up();
        } else {
            bfs_step_top_down();
        }
        
        // Synchronisation
        MPI_Barrier(MPI_COMM_WORLD);
        
        // Swap frontiers
        unsigned long *temp = frontier_bitmap;
        frontier_bitmap = next_frontier_bitmap;
        next_frontier_bitmap = temp;
        memset(next_frontier_bitmap, 0, visited_size * sizeof(unsigned long));
    }
    
    MPI_Barrier(MPI_COMM_WORLD);
}

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
    MPI_Allreduce(MPI_IN_PLACE, &edge_count, 1, MPI_INT64_T, MPI_SUM, MPI_COMM_WORLD);
    *edge_visit_count = edge_count;
}

void clean_pred(int64_t* pred) {
    memset(pred, -1, g.nlocalverts * sizeof(int64_t));
}

void free_graph_data_structure(void) {
    free_oned_csr_graph(&g);
    free(visited);
    free(frontier_bitmap);
    free(next_frontier_bitmap);
}

size_t get_nlocalverts_for_pred(void) {
    return g.nlocalverts;
}
