#ifndef BIPARTITE_H
#define BIPARTITE_H

#include <stdlib.h>
#include <string.h>

class BipartiteGraph {
public:
	int V;
	int* match; // pair of matching
	bool* used; // used in DFS
	int K; // #max matching = #min vertex cover
	int* vclist;
	bool dfs(int v,int V,int* xadj, int* adjncy){
		used[v] = true;
		for(int i = xadj[v]; i < xadj[v+1]; i++){
			int u = adjncy[i];
			int w = match[u];
			if(w < 0 || !used[w] && dfs(w,V,xadj,adjncy)){
				match[v] = u;
				match[u] = v;
				return true;
			}
		}
		return false;
	}
	BipartiteGraph(){
		V = 0;
		K = -1;
		match = NULL;
		used = NULL;
		vclist = NULL;
	}
	~BipartiteGraph(){
		if(match) free(match);
		if(used) free(used);
		if(vclist) free(vclist);
	}
	int MaxMatching(int V,int* xadj,int* adjncy){
		this->V = V;
		// xadj & adjncy <- same as Metis
		match = (int*)malloc(V*sizeof(int));
		used = (bool*)malloc(V*sizeof(bool));
		if(match == NULL || used == NULL){
			return -1;
		}
		int res = 0;
		memset(match,-1,V*sizeof(int));
		for(int v = 0; v < V; v++){
			if(match[v] < 0){
				memset(used,0,V*sizeof(bool));
				if(dfs(v,V,xadj,adjncy)){
					res++;
				}
			}
		}
		K = res;
		return res;
	}
	void mvc_dfs(int v,int X,int Y,int* xadj,int* adjncy){
		if(used[v]) return;
		used[v] = true;
		if(v < X){ // belong to X
			for(int i = xadj[v]; i < xadj[v+1]; i++){
				int idxY = adjncy[i];
				if(idxY != match[v]){ // (->) : not matching
					mvc_dfs(idxY,X,Y,xadj,adjncy);
				}
			}
		}
		else{ // belong to Y
			int idxX = match[v]; // (<-) : matching
			if(idxX < 0) return;
			mvc_dfs(idxX,X,Y,xadj,adjncy);
		}
	}
	int MinVertexCover(int X,int Y,int* xadj,int* adjncy){
		// get max matching
		if(this->V != X + Y || K == -1){
			K = MaxMatching(X+Y,xadj,adjncy);
		}
		vclist = (int*)malloc(K*sizeof(int*));
		memset(used,0,V*sizeof(bool)); // <- this array used for marking
		for(int i = 0; i < X; i++){
			if(match[i] < 0){
				// coloring vertices start from vertices in X which are not end vertices of matching
				// coloring vertices who can be reached from i
				// no match : (X->Y) direction, match : (X<-Y) direction
				mvc_dfs(i,X,Y,xadj,adjncy);
			}
		}

		int idx = 0;
		// select not marked
		for(int i = 0; i < X; i++){
			if(!used[i]){
				vclist[idx++] = i;
			}
		}
		// select marked
		for(int i = X; i < V; i++){
			if(used[i]){
				vclist[idx++] = i;
			}
		}
		if(idx != K){
			return -1; // something wrong	
		}
		return K;
	}
};

#endif
