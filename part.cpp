#include "part.h"

#include <stdio.h>
#include <stdlib.h>
#include <limits.h>
#include <set>
#include <map>
#include <queue>
#include <list>
#include <iterator>
#include <algorithm>
#include <assert.h>

/*
typedef struct{
	double ufactor; // allowed imbalance (1+x)/1000
} ndOptions;

class PartGraph {
private:
	int nvtxs;
	int nedges;
	int* xadj; // index to adjncy
	int* adjncy; // edge, size = xadj[nvtxs]
	int* vwgt; // size = nvtxs
	int* cewgt; // the wgt of edges that have been contacted to create v, size = nvtxs
	int* adjwgt; // size = xadj[nvtxs]
}
*/
void SetDefaultOptions(ndOptions* options){
	options->ufactor = 5; // 0.5% + 0.5% = 1%
	options->coarsenThreshold = 200;
}


/////////////////
// util functions
/////////////////

static void Iswap(int* x,int* y){
	int tmp = *x;
	*x = *y;
	*y = tmp;
}

static int GenerateRand(int n){
	return rand() % n;
}

static int RandomPermuteVector(int n,int* A){
	for(int i = 0; i < n; i++){
		A[i] = i;
	}
	for(int i = 0; i < n; i++){
		int idx = i + GenerateRand(n-i);
		Iswap(A+i,A+idx);
	}
	return 0;
}

static void InitMatch(int n, int* match, int** grpv, int** lrpv){
	// initialize
	for(int i = 0; i < n; i++){
		match[i] = i;
	}
	*grpv = (int*)malloc(n*sizeof(int));
	*lrpv = (int*)malloc(n*sizeof(int));
	RandomPermuteVector(n,*grpv);
}
static void finMatch(int* grpv,int* lrpv){
	free(grpv);
	free(lrpv);
}

static void SetVal(int n,int val, int* array){
	for(int i = 0; i < n; i++){
		array[i] = val;
	}
}
static void ICopy(int n,int* src,int* dst){
	for(int i = 0; i < n; i++){
		dst[i] = src[i];
	}
}

/////////////////
// sub functions
/////////////////
int PartGraph::RandomMatching(int* match){
	int* grpv, * lrpv;
	InitMatch(nvtxs,match,&grpv,&lrpv);

	for(int i = 0; i < nvtxs; i++){
		int vidx = grpv[i];
		if(match[vidx] == vidx){ // not matched yet
			int ne = xadj[vidx+1] - xadj[vidx];
			RandomPermuteVector(ne,lrpv);
			for(int j = 0; j < ne; j++){
				int eidx = adjncy[xadj[vidx] + lrpv[j]];
				if(match[eidx] == eidx){ // adj vert not matched yet -> matching
					match[eidx] = vidx;
					match[vidx] = eidx;
					break;
				}
			}
		}
	}

	finMatch(grpv,lrpv);
	return 0;
}
int PartGraph::HeavyEdgeMatching(int* match){
	int* grpv, * lrpv;
	InitMatch(nvtxs,match,&grpv,&lrpv);
	
	for(int i = 0; i < nvtxs; i++){
		int vidx = grpv[i];
		if(match[vidx] == vidx){ // not matched yet
			int mxwgt = 0;
			int mxidx = -1;
			for(int j = xadj[vidx]; j < xadj[vidx+1]; j++){
				int eidx = adjncy[j];
				if(ewgt[j] > mxwgt && match[eidx] == eidx){
					mxwgt = ewgt[j];
					mxidx = eidx;
				}
			}
			if(mxidx != -1){
				match[vidx] = mxidx;
				match[mxidx] = vidx;
			}
		}
	}

	finMatch(grpv,lrpv);
	return 0;
}

// if i < j -> map[i] <= map[j]
int PartGraph::Mapping(const int* match, int* map){
	// init
	for(int i = 0; i < nvtxs; i++){
		map[i] = -1;
	}

	int counter = 0;
	for(int i = 0; i < nvtxs; i++){
		if(map[i] == -1){
			map[i] = map[match[i]] = counter++;
		}
	}
	return counter;
}
int PartGraph::GenerateCoarserGraph(int newSize, const int* match, const int* map, PartGraph* newGraph){
	newGraph->nvtxs = newSize;
	newGraph->xadj = (int*)malloc((newSize+1)*sizeof(int));
	newGraph->vwgt = (int*)malloc(newSize*sizeof(int));
	newGraph->cewgt = (int*)malloc(newSize*sizeof(int));
	newGraph->adjwgt = (int*)malloc(newSize*sizeof(int));

	newGraph->totalvwgt = totalvwgt;
	newGraph->tolerance = tolerance;

	int* tAdjncy = (int*)malloc(nedges*sizeof(int));
	int* tEwgt = (int*)malloc(nedges*sizeof(int));

	int ecount = 0;
	for(int i = 0; i < nvtxs; i++){
		if(match[i] >= i){ // not processed yet
			std::map<int,int> st = std::map<int,int>();
			if(match[i] == i){ // only me
				// find adj
				for(int j = xadj[i]; j < xadj[i+1]; j++){
					int u = map[adjncy[j]];
					std::map<int,int>::iterator itr = st.find(u);
					if(itr == st.end()){
						st.insert(std::make_pair(u, ewgt[j]));
					}
					else{
						itr->second += ewgt[j];
					}
				}
				// set values
				newGraph->xadj[map[i]]   = ecount;
				newGraph->vwgt[map[i]]   = vwgt[i];
				newGraph->cewgt[map[i]]  = cewgt[i];
				newGraph->adjwgt[map[i]] = adjwgt[i];
			}
			else{ // matched
				//find adj
				for(int j = xadj[i]; j < xadj[i+1]; j++){
					int u = map[adjncy[j]];
					std::map<int,int>::iterator itr = st.find(u);
					if(itr == st.end()){
						st.insert(std::make_pair(u, ewgt[j]));
					}
					else{
						itr->second += ewgt[j];
					}
				}
				for(int j = xadj[match[i]]; j < xadj[match[i]+1]; j++){
					int u = map[adjncy[j]];
					std::map<int,int>::iterator itr = st.find(u);
					if(itr == st.end()){
						st.insert(std::make_pair(u, ewgt[j]));
					}
					else{
						itr->second += ewgt[j];
					}
				}
				st.erase(map[i]);
				// set values
				int wv1v2 = 0;
				for(int j = xadj[i]; j < xadj[i+1]; j++){
					if(adjncy[j] == match[i]){
						wv1v2 = ewgt[j];
						break;
					}
				}
				newGraph->xadj[map[i]]   = ecount;
				newGraph->vwgt[map[i]]   = vwgt[i] + vwgt[match[i]];
				newGraph->cewgt[map[i]]  = cewgt[i] + cewgt[match[i]] + wv1v2;
				newGraph->adjwgt[map[i]] = adjwgt[i] + adjwgt[match[i]] - 2*wv1v2;
			}
			// tAdjncy, tEwgt, ecount++
			for(std::map<int,int>::iterator im = st.begin(); im != st.end(); im++){
				tAdjncy[ecount] = im->first;
				tEwgt[ecount] = im->second;
				ecount++;
			}
		}
	}
	newGraph->xadj[newSize] = ecount;

	// set edge data
	newGraph->nedges = ecount;
	newGraph->adjncy = (int*)malloc(newGraph->nedges * sizeof(int));
	newGraph->ewgt = (int*)malloc(newGraph->nedges * sizeof(int));
	for(int i = 0; i < newGraph->nedges; i++){ // to save mem size , copy array
		newGraph->adjncy[i] = tAdjncy[i];
		newGraph->ewgt[i] = tEwgt[i];
	}

	free(tEwgt);
	free(tAdjncy);
	return 0;
}

int PartGraph::GGPartitioningEdge(double ratioX, int* partition){
	const int GGP_TIMES = 100;
	double ratio; // if ratio is close to 1, finding an unvisited value may be difficult in some cases
	int val;
	if(ratioX <= 0.5){ // consider X
		ratio = ratioX;
		val = 0;
	}
	else{ // consider Y
		ratio = 1 - ratioX;
		val = 1;
	}
	int nP = (int)(totalvwgt*ratio);

	int* part = (int*)malloc(nvtxs*sizeof(int));
	int bestScore = INT_MAX;

	for(int i = 0; i < GGP_TIMES; i++){
		SetVal(nvtxs,1 - val,part);
		std::queue<int> que;
		int counter = 0;
		while(counter < nP){
			while(1){
				int v = GenerateRand(nvtxs);
				if(part[v] == 1-val){
					que.push(v);
					break;
				}
			}
			while(!que.empty() && counter < nP){
				int v = que.front(); que.pop();
				if(part[v] == 1-val){
					part[v] = val;
					counter+= vwgt[v];
				}
				for(int j = xadj[v]; j < xadj[v+1]; j++){
					int u = adjncy[j];
					if(part[u] == 1-val){
						que.push(u);
					}
				}
			}
		}
		RefineEdge(ratioX,part);
		int newScore = GetEdgecut(part);
//		fprintf(stderr,"%d\n",newScore);
		if(bestScore > newScore){
			ICopy(nvtxs,part,partition);
			bestScore = newScore;
		}
	}
	free(part);
	edgecut = bestScore;
	currWgtX = 0;
	for(int i = 0; i < nvtxs; i++){
		if(partition[i] == 0) currWgtX += vwgt[i];
	}
	return bestScore;
}

int PartGraph::GetLargeGainVertexFromBoundary(std::vector<int> &list,int* partition,int val){
	// find large gain vertex
	int v = -1;
	if(list.empty()){
		while(1){
			v = GenerateRand(nvtxs);
			if(partition[v] != val){
				partition[v] = val;
				break;
			}
		}
		for(int i = xadj[v]; i < xadj[v+1]; i++){
			int u = adjncy[i];
			list.push_back(u);
		}
	}
	else{
		int i_max = -1;
		int g_max = INT_MIN;
		for(int i = 0; i < (int)list.size(); i++){
			int g = GetEdgeGain(list[i],partition);
			if(g > g_max){
				i_max = i;
				g_max = g;
			}
		}
		v = list[i_max];
		partition[v] = val;
		list.erase(list.begin()+i_max);
		
		// update list
		for(int i = xadj[v]; i < xadj[v+1]; i++){
			int u = adjncy[i];
			if(partition[u] != val && std::find(list.begin(),list.end(),u) == list.end()){
				list.push_back(u);
			}
		}
	}
	return v;
}

int PartGraph::GGGPartitioningEdge(double ratioX, int* partition){
	const int GGGP_TIMES = 100;
	double ratio;
	int val;
	if(ratioX <= 0.5){ // consider X
		ratio = ratioX;
		val = 0;
	}
	else{ // consider Y
		ratio = 1 - ratioX;
		val = 1;
	}
	int nP = (int)(totalvwgt*ratio);

	int* part = (int*)malloc(nvtxs*sizeof(int));
	int bestScore = INT_MAX;

	for(int i = 0; i < GGGP_TIMES; i++){
		SetVal(nvtxs,1 - val,part);
		std::vector<int> list;
		int counter = 0;
		while(counter < nP){
			int v = GetLargeGainVertexFromBoundary(list,part,val);
			counter+= vwgt[v];
		}
		RefineEdge(ratioX,part);
		int newScore = GetEdgecut(part);
//		fprintf(stderr,"%d\n",newScore);
		if(bestScore > newScore){
			ICopy(nvtxs,part,partition);
			bestScore = newScore;
		}
	}
	free(part);
	edgecut = bestScore;
	currWgtX = 0;
	for(int i = 0; i < nvtxs; i++){
		if(partition[i] == 0) currWgtX += vwgt[i];
	}

	return bestScore;
}

int PartGraph::GetEdgeGain(int v,int* partition){
	int g = 0;
	int pv = partition[v];
	for(int i = xadj[v]; i < xadj[v+1]; i++){
		int pu = partition[adjncy[i]];
		if(pu == pv){
			g -= ewgt[i];
		}
		else{
			g += ewgt[i];
		}
	}
	return g;
}

int PartGraph::GetEdgecut(int* partition){
	int score = 0;
	for(int i = 0; i < nvtxs; i++){
		int pi = partition[i];
		for(int j = xadj[i]; j < xadj[i+1]; j++){
			int pj = partition[adjncy[j]];
			if(pi != pj){
				score += ewgt[j];
			}
		}
	}
	return score / 2;
}

int PartGraph::RefineEdge(double ratioX, int* partition){
	FMDATA fm(this,partition);
	int totalWgt = totalvwgt;
	int minWgtX = (int)(totalWgt*ratioX - tolerance);
	int maxWgtX = (int)(totalWgt*ratioX + tolerance);
	currWgtX = fm.RefineEdge(currWgtX,minWgtX,maxWgtX,this,partition);
	fm.PrintSt();
	edgecut = fm.Edgecut();
	return edgecut;
}

/////////////////
// main functions
/////////////////

int PartGraph::Partition2(ndOptions* options, double ratioX, int* partition){
 // partition[nvtxs] // 0->X, 1->Y, 2->S
//	double ratioY = 1 - ratioX;
	SetTolerance(options);

	int* match = (int*)malloc(nvtxs*sizeof(int));
	int* map = (int*)malloc(nvtxs*sizeof(int));

	if(nvtxs <= options->coarsenThreshold || Esize() == 0){
		InitPartitioning(ratioX,partition);
	}
	else{
		int newSize = Coarsening(match,map);
		fprintf(stderr,"size = %d\n",newSize);

		int* coarserPart = (int*)malloc(newSize*sizeof(int));
		PartGraph newGraph;
		GenerateCoarserGraph(newSize,match,map,&newGraph);
		newGraph.Partition2(options,ratioX,coarserPart);		
	
		
		currWgtX = newGraph.currWgtX;
		edgecut = newGraph.edgecut;
		Uncoarsening(ratioX,map,coarserPart,partition);

		newGraph.DeleteGraph();
		free(coarserPart);
	}
		
	free(map);
	free(match);
	fprintf(stderr,"Xwgt %d (%f), Ecut %d (size %d)\n",currWgtX,1.0*currWgtX/totalvwgt,edgecut,nvtxs);

		int dbg = 0;
		for(int i = 0; i < nvtxs; i++){
			if(partition[i] == 0) dbg += vwgt[i];
		}
		fprintf(stderr,"wgtx true: %d\n",dbg);

	return 0; 
}


int PartGraph::Coarsening(int* match,int* map){
	// matching and map
//	RandomMatching(match);
	HeavyEdgeMatching(match);

	int newSize = Mapping(match,map);
	return newSize;
}

int PartGraph::Uncoarsening(double ratioX, int* map,int* coarserPart,int* partition){
	// project
	for(int i = 0; i < nvtxs; i++){
		partition[i] = coarserPart[map[i]];
	}

	// refinement
	int new_ecut = RefineEdge(ratioX,partition);
	assert(edgecut == GetEdgecut(partition));

	return new_ecut;
}

int PartGraph::UpperEdgeGain(){
	int ret = INT_MIN;
	for(int i = 0; i < nvtxs; i++){
		int g = 0;
		for(int j = xadj[i]; j < xadj[i+1]; j++){
			g += ewgt[j];
		}
		if(g > ret) ret = g;
	}
	return ret;
}
int PartGraph::MaxVertexWeight(){
	int ret = 0;
	for(int i = 0; i < nvtxs; i++){
		if(vwgt[i] > ret) ret = vwgt[i];
	}
	return ret;
}
void PartGraph::SetTolerance(ndOptions* options){
	int smw = 0;
	int mxw = 0;
	for(int i = 0; i < nvtxs; i++){
		 smw += vwgt[i];
		 if(vwgt[i] > mxw) mxw = vwgt[i];
	}
	tolerance = (int)(smw / 1000.0 * options->ufactor);
	if(tolerance < mxw){
		 tolerance = mxw;
	}
}

void PartGraph::Show(){
	printf("%d %d\n",nvtxs,nedges);
	for(int i = 0; i < nvtxs; i++){
		printf("%d %d\n",i,vwgt[i]);
	}
	for(int i = 0; i < nvtxs; i++){
		for(int j = xadj[i]; j < xadj[i+1]; j++){
			printf("%d %d %d\n",i,adjncy[j],ewgt[j]);
		}
	}
}
void PartGraph::Show(int* partition){
	printf("%d %d\n",nvtxs,nedges);
	for(int i = 0; i < nvtxs; i++){
		printf("%d %d %d\n",i,vwgt[i],partition[i]);
	}
	for(int i = 0; i < nvtxs; i++){
		for(int j = xadj[i]; j < xadj[i+1]; j++){
			printf("%d %d %d\n",i,adjncy[j],ewgt[j]);
		}
	}
}

int PartGraph::InitPartitioning(double ratioX, int* partition){
	return GGPartitioningEdge(ratioX,partition);
//	return GGGPartitioningEdge(ratioX,partition);
}
