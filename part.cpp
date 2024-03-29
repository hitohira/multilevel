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
#include <functional>
#include <assert.h>

#include "bipartite.h"

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
	options->matchingScheme = MATCHING_RM;
	options->initPartScheme = INITPART_GGGP;
	options->refineScheme = REFINE_EFM;
}


/////////////////
// util functions
/////////////////
static void dbg(int n){
	printf("dbg %d\n",n);
}

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
int PartGraph::LiteEdgeMatching(int* match){
	int* grpv, * lrpv;
	InitMatch(nvtxs,match,&grpv,&lrpv);
	
	for(int i = 0; i < nvtxs; i++){
		int vidx = grpv[i];
		if(match[vidx] == vidx){ // not matched yet
			int mnwgt = INT_MAX;
			int mnidx = -1;
			for(int j = xadj[vidx]; j < xadj[vidx+1]; j++){
				int eidx = adjncy[j];
				if(ewgt[j] < mnwgt && match[eidx] == eidx){
					mnwgt = ewgt[j];
					mnidx = eidx;
				}
			}
			if(mnidx != -1){
				match[vidx] = mnidx;
				match[mnidx] = vidx;
			}
		}
	}

	finMatch(grpv,lrpv);
	return 0;
}


int PartGraph::ClusteringMatchingAndMapping(MatchData &md, int* map){
	int counter_match = 0;
	int map_len = 0;
	int threshold_stop = nvtxs / 2;
	int threshold_weigt = totalvwgt / 20;

	int* grpv = (int*)malloc(nvtxs*sizeof(int));
	RandomPermuteVector(nvtxs,grpv);

	for(int i = 0; i < nvtxs; i++){
		map[i] = -1;
	}
	for(int i = 0; i < nvtxs; i++){
		int vidx = grpv[i];

		if(map[vidx] != -1) continue; // already matched, so you shold do nothing

		// v is not matched yet bellow
		
		if(counter_match >= threshold_stop){ // so many vertices has already been matched
			map[vidx] = map_len;
			md.Add(map_len,vidx,this);
			map_len++;
			continue;
		}
		if(vwgt[vidx] >= threshold_weigt){ // v is too heavy
			map[vidx] = map_len;
			md.Add(map_len,vidx,this);
			map_len++;
			continue;
		}

		// try match
		// w(e) max -> w(u) min
		int mxewgt = INT_MIN;
		int mnuwgt = INT_MAX;
		int idx = -1;
		for(int j = xadj[vidx]; j < xadj[vidx+1]; j++){
			int uidx = adjncy[j];
			// limit cluster size
			if(vwgt[uidx] + vwgt[vidx] >= threshold_weigt){
				continue;
			}
			if(map[uidx] != -1 && md.Weight(map[uidx]) + vwgt[vidx] >= threshold_weigt){
				continue;
			}

			// get best edge
			if(ewgt[j] > mxewgt){
				idx = j;
				mxewgt = ewgt[j];
				mnuwgt = vwgt[uidx];
			}
			else if(ewgt[j] == mxewgt && vwgt[uidx] < mnuwgt){
				idx = j;
				mxewgt = ewgt[j];
				mnuwgt = vwgt[uidx];
			}
		}

		if(idx != -1){ // match OK
			counter_match++;

			int uidx = adjncy[idx];
			assert(uidx < nvtxs);
			if(map[uidx] == -1){ // u has not been matched yet
				map[vidx] = map[uidx] = map_len;
				md.Add(map_len,vidx,this);
				md.Add(map_len,uidx,this);
				map_len++;
			}
			else{ // u has already matched
				map[vidx] = map[uidx];
				md.Add(map[uidx],vidx,this);
			}
		}
		else { // not matched (because all u are too heavy)
			map[vidx] = map_len;
			md.Add(map_len,vidx,this);
			map_len++;
		}
	}
	free(grpv);
	return map_len;
}

int PartGraph::GenerateClusteringCoarserGraph(MatchData &md, int* map, PartGraph* newGraph){
	int newSize = md.GetMappedSize();
/*
	if(nvtxs == newSize){
		printf("dump\n");
		for(int i = 0; i < nvtxs; i++){
			printf("%d %d %d\n",i,vwgt[i],xadj[i+1] - xadj[i]);
		}
		exit(0);
	}
*/
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
	for(int i = 0; i < newSize; i++){
		newGraph->xadj[i] = ecount;
		newGraph->vwgt[i] = md.Weight(i);

		int sum_cewgt = 0;
		int sum_adjwgt = 0;
		int len = md.MatchSize(i);
		std::map<int,int> st = std::map<int,int>();
		for(int j = 0; j < len; j++){
			int v = md.Match(i,j);
			sum_cewgt += cewgt[v];
			sum_adjwgt += adjwgt[v];
			for(int k = xadj[v]; k < xadj[v+1]; k++){
				int uidx = adjncy[k];
				int vuwgt = ewgt[k];
				int umap = map[uidx];
				if(umap == map[v]) continue; // same coarsened vertex
				std::map<int,int>::iterator itr = st.find(umap);
				if(itr == st.end()){
					st.insert(std::make_pair(umap,vuwgt));
				}
				else{
					itr->second += vuwgt;
				}
			}
		}
		for(std::map<int,int>::iterator im = st.begin(); im != st.end(); im++){
			tAdjncy[ecount] = im->first;
			tEwgt[ecount] = im->second;
			ecount++;
		}

		int wv1v2Expand = md.Cewgt(i);
//		newGraph->cewgt[i] = sum_cewgt + wv1v2Expand; // for original
		newGraph->cewgt[i] = sum_cewgt; // for #vert version
		newGraph->adjwgt[i] = sum_adjwgt - 2*wv1v2Expand; 
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
			//	newGraph->cewgt[map[i]]  = cewgt[i] + cewgt[match[i]] + wv1v2;
				newGraph->cewgt[map[i]]  = cewgt[i] + cewgt[match[i]];
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

int PartGraph::GGPartitioningEdge(ndOptions* options,double ratioX, int* partition){
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
		currWgtX = GetCurrWgt(part);
		RefineEdge(options,ratioX,part);
		int newScore = GetEdgecut(part);
/* vert ver
		edgecut = GetEdgecut(part);
		VertSepFromEdgeSep(part);
		SetWgtInfo(this,part,ratioX,tolerance,&wgtInfo);
		RefineVert(part);
		int newScore = GetVertSepSize(part);
	*/
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
	fprintf(stderr,"wgt %d (%f)\n",currWgtX,1.0*currWgtX/totalvwgt);
	return bestScore;
}

int PartGraph::GetLargeGainVertexFromBoundary(std::vector<int> &list,int* partition,int val){
	// find large gain vertex
	int v = -1;
	if(list.empty()){ // for some connected part
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

int PartGraph::GGGPartitioningEdge(ndOptions* options, double ratioX, int* partition){
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
		int tmp2 = 0;
		for(int k = 0; k < nvtxs; k++){
			if(part[k] == 0) tmp2 += vwgt[k];
		}
		currWgtX = GetCurrWgt(part);
		RefineEdge(options,ratioX,part);
		int tmp = 0;
		for(int k = 0; k < nvtxs; k++){
			if(part[k] == 0) tmp += vwgt[k];
		}
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
	fprintf(stderr,"wgt %d (%f)\n",currWgtX,1.0*currWgtX/totalvwgt);
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

int PartGraph::GetVertGain(int v,int to, int* partition){
	// partition[v] == 2
	// to = 0 or 1
	int g = vwgt[v];
	for(int i = xadj[v]; i < xadj[v+1]; i++){
		int u = adjncy[i];
		if(partition[u] == (1-to)){
			g -= vwgt[u];
		}
	}
	return g;
}
int PartGraph::GetVertGainCewgt(int v,int to, int* partition){
	// partition[v] == 2
	// to = 0 or 1
	int g = cewgt[v];
	for(int i = xadj[v]; i < xadj[v+1]; i++){
		int u = adjncy[i];
		if(partition[u] == (1-to)){
			g -= cewgt[u];
		}
	}
	return g;
}
int PartGraph::GetVertSepSize(int* partition){
	int ret = 0;
	for(int i = 0; i < nvtxs; i++){
		if(partition[i] == 2){
			ret += vwgt[i];
		}
	}
	return ret;
}
int PartGraph::GetVertSepSizeCewgt(int* partition){
	int ret = 0;
	for(int i = 0; i < nvtxs; i++){
		if(partition[i] == 2){
			ret += cewgt[i];
		}
	}
	return ret;
}


int PartGraph::RefineEdge(ndOptions* options, double ratioX, int* partition){
	FMDATA fm(this,partition);
	int totalWgt = totalvwgt;
	int minWgtX = (int)(totalWgt*ratioX - tolerance);
	int maxWgtX = (int)(totalWgt*ratioX + tolerance);
	currWgtX = fm.RefineEdge(currWgtX,minWgtX,maxWgtX,this,partition);
	fm.PrintSt();
	edgecut = fm.Edgecut();
	return edgecut;
}

int PartGraph::RefineVert(ndOptions* options, int* partition){
	int vert_sep_size;
	if(options->refineScheme == REFINE_VFM2){
		FMDATAvert2 fm(this,partition);
		vert_sep_size = fm.RefineVert(&wgtInfo,this,partition);
		fm.PrintSt();
	}
	else{
		FMDATAvert fm(this,partition);
		vert_sep_size = fm.RefineVert(&wgtInfo,this,partition);
		fm.PrintSt();
	}
	return vert_sep_size;
}

bool PartGraph::VertSepIsOK(int* partition){
	for(int i = 0; i < nvtxs; i++){
		if(partition[i] == 2) continue; // separator is skipped
		int opo = 1 - partition[i];
		for(int j = xadj[i]; j < xadj[i+1]; j++){
			int u = adjncy[j];
			if(partition[u] == opo) return false;
		}
	}
	return true;
}

int PartGraph::GetBoundary(std::vector<int>& bvs,std::map<int, int>& ibvs,int* partition){
	bvs = std::vector<int>();
	std::vector<int> bvs2;

	for(int v = 0; v < nvtxs; v++){
		int pv = partition[v];
		for(int i = xadj[v]; i < xadj[v+1]; i++){
			int u = adjncy[i];
			int pu = partition[u];
			if(pv != pu){
				if(pv == 0){
					bvs.push_back(v);
				}
				else{
					bvs2.push_back(v);
				}
				break;
			}
		}
	}
	int X = bvs.size();
	bvs.insert(bvs.end(), bvs2.begin(),bvs2.end());
	for(int i = 0; i < (int)bvs.size(); i++){
		ibvs[bvs[i]] = i;
	}
	return X;
}
int PartGraph::VertSepFromEdgeSep(int* partition){
	std::vector<int> bvs;
	std::map<int, int> ibvs;
	int X = GetBoundary(bvs,ibvs,partition);
	int Y = bvs.size() - X;
	int nv = bvs.size();
	int ne = edgecut*2;
	int* bxadj = (int*)malloc((nv+1)*sizeof(int));
	int* badjncy = (int*)malloc(ne*sizeof(int));

	bxadj[0] = 0;
	int cntr = 0;
	for(int i = 0; i < nv; i++){
		int v = bvs[i];
		for(int j = xadj[v]; j < xadj[v+1]; j++){
			int u = adjncy[j];
			if(partition[v] != partition[u]){
				int iu = ibvs[u];
				assert(iu < bvs.size()); ///////////////// 
				badjncy[cntr++] = iu;
			}
		}
		bxadj[i+1] = cntr;
	}
	BipartiteGraph bg;
	bg.MinVertexCover(X,Y,bxadj,badjncy);
	// bg.vclist is MVC
		
	for(int i = 0; i < bg.K; i++){
		int idx = bvs[bg.vclist[i]];
		partition[idx] = 2;
	}

	free(bxadj);
	free(badjncy);	
	return bg.K; // vertex separator size
}


/////////////////
// main functions
/////////////////

int PartGraph::NestedDissection(ndOptions* options, Etree& etree, int epos, int* partition, int* perm){
	fprintf(stderr,"nd node %d\n", epos);
	Enode enode = etree.Get(epos);
	
	if(enode.left == -1){
		for(int i = 0; i < nvtxs; i++) partition[i] = epos;
		fprintf(stderr,"CM start %d\n", epos);
		CuthillMcKee(perm);
		fprintf(stderr,"CM end %d\n", epos);
		return 0;
	}
	Enode eleft = etree.Get(enode.left);
	//Enode eright = etree.Get(enode.right);

	double ratioX = 1.0 * eleft.weight / enode.weight;
	Partition3(options, ratioX, partition);

	// divide Graph
	PartGraph leftGraph, rightGraph;
	int* leftMap = (int*)malloc(nvtxs*sizeof(int));
	int* rightMap = (int*)malloc(nvtxs*sizeof(int));
	DivideGraphByPartition(partition, leftGraph, leftMap, rightGraph, rightMap);

	int* leftPartition = (int*)malloc(leftGraph.Vsize()*sizeof(int));
	int* rightPartition = (int*)malloc(rightGraph.Vsize()*sizeof(int));
	int* leftPerm = (int*)malloc(leftGraph.Vsize()*sizeof(int));
	int* rightPerm = (int*)malloc(rightGraph.Vsize()*sizeof(int));
	leftGraph.NestedDissection(options, etree, enode.left, leftPartition, leftPerm);
	rightGraph.NestedDissection(options, etree, enode.right, rightPartition, rightPerm);

	// merge result
	ndMergeInfo leftMi;
	leftMi.size = leftGraph.Vsize();
	leftMi.perm = leftPerm;
	leftMi.partition = leftPartition;
	leftMi.map = leftMap;
	ndMergeInfo rightMi;
	rightMi.size = rightGraph.Vsize();
	rightMi.perm = rightPerm;
	rightMi.partition = rightPartition;
	rightMi.map = rightMap;
	MergeInfo(epos, leftMi, rightMi, partition, perm);

	free(leftMap);
	free(rightMap);
	free(leftPartition);
	free(rightPartition);
	free(leftPerm);
	free(rightPerm);
	return 0;
}

int PartGraph::DivideGraphByPartition(int* partition, PartGraph& left, int* leftMap, PartGraph& right, int* rightMap){
	left.DeleteGraph();
	right.DeleteGraph();

	// count and mapping
	int lnv = 0, rnv = 0;
	int lnnz = 0, rnnz = 0;
	int* mapper = (int*)malloc(nvtxs*sizeof(int));
	for(int i = 0; i < nvtxs; i++){
		int nnz = 0;
		for(int j = xadj[i]; j < xadj[i+1]; j++){
			int col = adjncy[j];
			if(partition[col] == partition[i]) nnz++;
		}

		if(partition[i] == 0){
			mapper[i] = lnv;
			leftMap[lnv] = i;
			lnv++;
			lnnz += nnz;
		}
		else if(partition[i] == 1){
			mapper[i] = rnv;
			rightMap[rnv] = i;
			rnv++;
			rnnz += nnz;
		}
		else{
			mapper[i] = -1;
		}
	}
	// alloc 
	left.Allocate(lnv, lnnz);
	right.Allocate(rnv, rnnz);

	// set data 
	int lvpos = 0, rvpos = 0;
	int lnzpos = 0, rnzpos = 0;
	left.xadj[lnv] = lnnz;
	right.xadj[rnv] = rnnz;
	left.totalvwgt = 0;
	right.totalvwgt = 0;
	for(int i = 0; i < nvtxs; i++){
		if(partition[i] == 0){
			left.xadj[lvpos] = lnzpos;
			left.vwgt[lvpos] = vwgt[i];
			left.totalvwgt += vwgt[i];
			left.cewgt[lvpos] = cewgt[i];

			left.adjwgt[lvpos] = 0;
			for(int j = xadj[i]; j < xadj[i+1]; j++){
				int u = adjncy[j];
				// if adj vert is in the same partition
				if(partition[u] == 0){
					left.adjwgt[lvpos] += ewgt[j];
					left.adjncy[lnzpos] = mapper[u];
					left.ewgt[lnzpos] = ewgt[j];
					lnzpos++;
				}
			}
			lvpos++;
		}
		else if(partition[i] == 1){
			right.xadj[rvpos] = rnzpos;
			right.vwgt[rvpos] = vwgt[i];
			right.totalvwgt += vwgt[i];
			right.cewgt[rvpos] = cewgt[i];

			right.adjwgt[rvpos] = 0;
			for(int j = xadj[i]; j < xadj[i+1]; j++){
				int u = adjncy[j];
				// if adj vert is in the same partition
				if(partition[u] == 1){
					right.adjwgt[rvpos] += ewgt[j];
					right.adjncy[rnzpos] = mapper[u];
					right.ewgt[rnzpos] = ewgt[j];
					rnzpos++;
				}
			}
			rvpos++;
		}
	}
	assert(lnv == lvpos);
	assert(rnv == rvpos);
	assert(lnnz == lnzpos);
	assert(rnnz == rnzpos);

	fprintf(stderr,"Nv  %d %d %d\n", nvtxs - lnv - rnv, lnv, rnv);
	fprintf(stderr,"Nnz %d %d %d\n", xadj[nvtxs] - lnnz - rnnz, lnnz, rnnz);

	free(mapper);
	return 0;
}

// map : small -> large
int PartGraph::MergeInfo(int epos, ndMergeInfo& left, ndMergeInfo& right, int* partition, int* perm){
	int sepidx = 0;
	for(int i = 0; i < nvtxs; i++){
		if(partition[i] == 2){
			partition[i] = epos;
			perm[i] = sepidx++;
		}
	}

	for(int i = 0; i < left.size; i++){
		perm[left.map[i]] = left.perm[i];
		partition[left.map[i]] = left.partition[i];
	}
	for(int i = 0; i < right.size; i++){
		perm[right.map[i]] = right.perm[i];
		partition[right.map[i]] = right.partition[i];
	}
	return 0;
}

typedef std::pair<int, int> Pair;

static Pair FindMinDegVertex(std::set<Pair> noVisit){
		return *(noVisit.begin());
}

int PartGraph::CuthillMcKee(int* perm){
	for(int i = 0; i < nvtxs; i++){
		perm[i] = -1;
	}
	int inc = 0;
	std::set<Pair> noVisit;
	std::queue<Pair> que;
	std::vector<Pair> sortV;

	for(int i = 0; i < nvtxs; i++){
		int deg = xadj[i+1] - xadj[i];
		if(deg == 0){ // this is very important for speed up. why???
			perm[i] = inc++;
		}
		else{
			noVisit.insert(std::make_pair(deg,i));
		}
	}

	long long perf_outer = 0;
	long long perf_inner = 0;

	while(!noVisit.empty()){
		perf_outer++;
		Pair firstPair = FindMinDegVertex(noVisit);
		que.push(firstPair);
		perm[firstPair.second] = -2;// -2 : added, -1 not added, 0-(m-1): processed
		noVisit.erase(noVisit.begin());

		while(!que.empty()){
			perf_inner++;
			Pair p = que.front(); que.pop();
			int idx = p.second;
			perm[idx] = inc++;
			std::vector<Pair>().swap(sortV);
			for(int i = xadj[idx]; i < xadj[idx+1]; i++){
				int col = adjncy[i];
				if(perm[col] == -1){ // if not added to que yet
					int deg = xadj[col+1] - xadj[col];
					sortV.push_back(std::make_pair(deg,col));
					perm[col] = -2;
					noVisit.erase(std::make_pair(deg,col));
				}
			}
			std::sort(sortV.begin(),sortV.end());
			for(int i = 0; i < (int)sortV.size(); i++){
				que.push(sortV[i]);
			}
		}
	}
	fprintf(stderr,"outer : %lld, inner : %lld\n",perf_outer, perf_inner);
	return 0;
}

int PartGraph::Partition2(ndOptions* options, double ratioX, int* partition){
	Partition2Inner(options,ratioX,partition);
	RefineEdge(options,ratioX,partition);
	return 0;
}
int PartGraph::Partition2Inner(ndOptions* options, double ratioX, int* partition){
 // partition[nvtxs] // 0->X, 1->Y
//	double ratioY = 1 - ratioX;
	SetTolerance(options);

	int* map = (int*)malloc(nvtxs*sizeof(int));

	if(nvtxs <= options->coarsenThreshold || Esize() == 0){
		InitPartitioningEdge(options,ratioX,partition);
	}
	else{
		int newSize;
		PartGraph newGraph;
		if(options->matchingScheme == MATCHING_CLUSTER){ 
			MatchData md(nvtxs);
			newSize = ClusteringMatchingAndMapping(md,map);
			GenerateClusteringCoarserGraph(md,map,&newGraph);
		}
		else{
			int* match = (int*)malloc(nvtxs*sizeof(int));
			newSize = Coarsening(options,match,map);
			GenerateCoarserGraph(newSize,match,map,&newGraph);
			free(match);
		}

		fprintf(stderr,"size = %d\n",newSize);
		
		if(newSize != nvtxs){
			int* coarserPart = (int*)malloc(newSize*sizeof(int));
			newGraph.Partition2(options,ratioX,coarserPart);		
			
			currWgtX = newGraph.currWgtX;
			edgecut = newGraph.edgecut;
			UncoarseningEdge(options,ratioX,map,coarserPart,partition);

			newGraph.DeleteGraph();
			free(coarserPart);
		}
		else{
			InitPartitioningEdge(options,ratioX,partition);
		}
	}
		
	free(map);
	fprintf(stderr,"Xwgt %d (%f), Ecut %d (size %d)\n",currWgtX,1.0*currWgtX/totalvwgt,edgecut,nvtxs);
/*
		int dbg = 0;
		for(int i = 0; i < nvtxs; i++){
			if(partition[i] == 0) dbg += vwgt[i];
		}
		fprintf(stderr,"wgtx true: %d\n",dbg);
*/
	return 0; 
}

int PartGraph::Partition3Inner(ndOptions* options, double ratioX, int* partition){
 // partition[nvtxs] // 0->X, 1->Y, 2->S
//	double ratioY = 1 - ratioX;
	SetTolerance(options);

	int* map = (int*)malloc(nvtxs*sizeof(int));

	if(nvtxs <= options->coarsenThreshold || Esize() == 0){
		InitPartitioningVert(options,ratioX,partition);
		SetWgtInfo(this,partition,ratioX,tolerance,&wgtInfo);
	}
	else{
		int newSize;	
		PartGraph newGraph;
		if(options->matchingScheme == MATCHING_CLUSTER){ 
			MatchData md(nvtxs);
			newSize = ClusteringMatchingAndMapping(md,map);
			GenerateClusteringCoarserGraph(md,map,&newGraph);
		}
		else{
			int* match = (int*)malloc(nvtxs*sizeof(int));
			newSize = Coarsening(options,match,map);
			GenerateCoarserGraph(newSize,match,map,&newGraph);
			free(match);
		}
		fprintf(stderr,"size = %d\n",newSize);
		
		if(nvtxs != newSize){
			int* coarserPart = (int*)malloc(newSize*sizeof(int));
			newGraph.Partition3Inner(options,ratioX,coarserPart);		
	
			wgtInfo = newGraph.wgtInfo;
			wgtInfo.tol = 1.0*tolerance / totalvwgt;
			UncoarseningVert(options,ratioX,map,coarserPart,partition);

			newGraph.DeleteGraph();
			free(coarserPart);
		}
		else{
			InitPartitioningVert(options,ratioX,partition);
			SetWgtInfo(this,partition,ratioX,tolerance,&wgtInfo);
		}
	}
		
	free(map);
	fprintf(stderr,"Xwgt %d (%f), Ssize %d (%f) (size %d)\n",wgtInfo.wgt[0],1.0*wgtInfo.wgt[0]/totalvwgt,wgtInfo.wgt[2],1.0*wgtInfo.wgt[2]/totalvwgt,nvtxs);
	return 0; 
}
int PartGraph::Partition3(ndOptions* options, double ratioX, int* partition){
	Partition3Inner(options,ratioX,partition);
	RefineVert(options,partition);
	return 0;
}


int PartGraph::Coarsening(ndOptions* options, int* match,int* map){
	// matching and map
	if(options->matchingScheme == MATCHING_RM){
		RandomMatching(match);
	}
	else if(options->matchingScheme == MATCHING_HEM){
		HeavyEdgeMatching(match);
	}
	else if(options->matchingScheme == MATCHING_LEM){
		LiteEdgeMatching(match);
	}

	int newSize = Mapping(match,map);
	return newSize;
}

int PartGraph::UncoarseningEdge(ndOptions* options,double ratioX, int* map,int* coarserPart,int* partition){
	// project
	for(int i = 0; i < nvtxs; i++){
		partition[i] = coarserPart[map[i]];
	}

	// refinement
	int new_ecut = RefineEdge(options,ratioX,partition);
//	assert(edgecut == GetEdgecut(partition));

	return new_ecut;
}

int PartGraph::UncoarseningVert(ndOptions* options, double ratioX, int* map,int* coarserPart,int* partition){
	// project
	for(int i = 0; i < nvtxs; i++){
		partition[i] = coarserPart[map[i]];
	}

	// refinement
	int vert_sep_size = RefineVert(options,partition);
//	assert(vert_sep_size == GetVertSepSize(partition));

	return vert_sep_size;
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
int PartGraph::UpperVertGain(){
	int ret = INT_MIN;
	for(int i = 0; i < nvtxs; i++){
		int g = vwgt[i];
		for(int j = xadj[i]; j < xadj[i+1]; j++){
			int u = adjncy[j];
			g += vwgt[u];
		}
		if(g > ret) ret = g;
	}
	return ret;
}
int PartGraph::UpperVertGainCewgt(){
	int ret = INT_MIN;
	for(int i = 0; i < nvtxs; i++){
		int g = cewgt[i];
		for(int j = xadj[i]; j < xadj[i+1]; j++){
			int u = adjncy[j];
			g += cewgt[u];
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
//	if(tolerance < mxw){
//		 tolerance = mxw;
//	}
}

void PartGraph::Allocate(int nvtxs, int nedges){
	this->nvtxs = nvtxs;
	this->nedges = nedges;
	this->xadj = (int*)malloc((nvtxs+1)*sizeof(int));
	this->adjncy = (int*)malloc(nedges*sizeof(int));
	this->vwgt = (int*)malloc(nvtxs*sizeof(int));
	this->ewgt = (int*)malloc(nedges*sizeof(int));
	this->cewgt = (int*)malloc(nvtxs*sizeof(int));
	this->adjwgt = (int*)malloc(nvtxs*sizeof(int));
	tolerance = totalvwgt = currWgtX = edgecut = 0;
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

int PartGraph::GetCurrWgt(int* partition){
	int ret = 0;
	for(int i = 0; i < nvtxs; i++){
		if(partition[i] == 0){
			ret += vwgt[i];
		}
	}
	return ret;
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

int PartGraph::InitPartitioningEdge(ndOptions* options, double ratioX, int* partition){
	if(options->initPartScheme == INITPART_GGP){
		return GGPartitioningEdge(options,ratioX,partition);
	}
	return GGGPartitioningEdge(options,ratioX,partition);
}

int PartGraph::InitPartitioningVert(ndOptions* options, double ratioX, int* partition){
//	GGPartitioningEdge(options,ratioX,partition);
//	VertSepFromEdgeSep(partition);
	for(int i = 0; i < nvtxs; i++) partition[i] = 2;

	SetWgtInfo(this,partition,ratioX,tolerance,&wgtInfo);
	int vert_sep_size = RefineVert(options,partition);
//	assert(vert_sep_size == GetVertSepSize(partition));
	fprintf(stderr,"Ssize %d\n",vert_sep_size);
	return vert_sep_size;
}
