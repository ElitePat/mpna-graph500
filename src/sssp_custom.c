// Stub for custom SSSP implementation

#include "aml.h"
#include "common.h"
#include "csr_reference.h"
#include "bitmap_reference.h"

extern oned_csr_graph g;
extern int64_t* column,*pred_glob,visited_size;
extern unsigned long * visited;

#ifdef SSSP

//user provided function to be called several times to implement kernel 3: single source shortest path
//function has filled with -1 and 1.0 pred and dist arrays
//at exit dist should have shortest distance to root for each vertice otherwise -1
//pred array should point to next vertie in shortest path or -1 if vertex unreachable
//pred[VERTEX_LOCAL(root)] should be root and dist should be 0.0

void run_sssp(int64_t root,int64_t* pred,float *dist) {
	pred_glob=pred;
	//user code for SSSP
	// initialisation phase already done

	// Varaibles
	int64_t *frontier_pred; // frontière precedente
	int current_pred;
	int64_t *frontier_next; // frontière suivante
	int current_next = 0;
	float new_dist;
	int visited[g.max_nlocalverts] = {0}; // somets visités

	// to begin with
	pred[root] = root;
	dist[root] = 0.0;
	visited[root] = 1;
	// intialisation de la frontière
	frontier_pred[0] = root;
	current=1;

	// Boucle principale
	while(true){

		// on elargi la frontière
		for(int j=0; j<current_pred; j++){ // j -> tous les noeuds dans la frontière
			for(int i=0; i<g.max_nlocalverts; i++){ // i -> tous les noeuds du graphe
				// La frontièr est faite de noeuds non visités
				if(visited[i] != 1){ 
					if(g.column[i] == frontier_pred[j]){ // tout ceux qui sont voisins avec le noeud (de la frontière) actuel
						frontier_next[current_next] = i; // on les rajoute à le frontière
						current_next++;
					}
				}
			}
		}

		// la condition d'arret ! (une frontière vide)
		if(cuurent_next == 0){
			break;
		}

		// on découvre la nouvelle frontière
		for(int i=0; i<current_next; i++){
			// on comapre la distance depuis le noeud actuel à celle qu'il a deja
			new_dist= dist[frontier_next[i]] + g.weight[frontier_next[i]];
			// si il n'a jamais été visité.
			if(visited[frontier_next[i]] == 0){
				dist[frontier_next[i]] = new_dist;
				pred[frontier_next[i]] = frontier_next[i];
				// on marque le noeud comme visité
				visited[frontier_next[i]] = 1;
			}else if(new_dist < dist[frontier_next[i]]){ // si il a deja été visité
				dist[frontier_next[i]] = new_dist;
				pred[frontier_next[i]] = frontier_next[i];
			} // et sinon on fait rien
		}

		// une fois exploré la nouvelle frontière, elle n'est plus nouvelle elle devien l'ancienne
		for(int i=0;i<current_next;i++){
			frontier_pred[i] = frontier_next[i];
		}
		current_pred = current_next;

		// on vide la nouvelle frontière (pour nouvelle exploration)
		for(int i=0;i<current_next;i++){
			frontier_next[i] = 0;
		}
		current_next = 0;
	}

	// now time for relaxation
	// on visite tous les noeuds même ceux pas encore visités
	for(int i=0; i<g.max_nlocalverts; i++){
		for(int j=0; j<g.max_nlocalverts; j++){
			if(g.column[j] == i){// pour chacun de ses voisins
				new_dist = dist[j] + g.weight[j];
				if(new_dist < dist[i]){ // on cherche si il y pas un plus court chemin !
					// si c'est le cas on actualise
					dist[i] = new_dist;
					pred[i] = j; // le voisin en question devien son predecesseur
				} // sinon on fait rien
			}
		}
	}

}

//user provided function to prefill dist array with whatever value
void clean_shortest(float* dist) {
	int i;
	for(i=0;i<g.nlocalverts;i++) dist[i]=-1.0;
}
#endif
