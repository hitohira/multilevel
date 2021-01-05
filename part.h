#ifndef MY_PART_H
#define MY_PART_H

#include <stdlib.h>
#include <limits.h>
#include <vector>
#include <map>

#include "wgt.h"
#include "FM.h"

typedef struct ndOptions{
	int ufactor; // allowed imbalance x/1000
	int coarsenThreshold;
} ndOptions;

void SetDefaultOptions(ndOptions* options);

class PartGraph {
private:
	int nvtxs;
	int nedges;
	int* xadj; // index to adjncy
	int* adjncy; // edge, size = xadj[nvtxs]
	int* vwgt; // size = nvtxs
	int* ewgt; // size = xadj[nvtxs]
//	int* cewgt; // not used // the wgt of edges that have been contacted to create v, size = nvtxs
	int* cewgt; // #vertices contacted to create v, size = nvtxs, load.cpp is also changed
	int* adjwgt; // not used // the sum of wgt of the edges adjcent to v,  size = nvtxs

	int tolerance;
	int totalvwgt;
	int currWgtX;
	int edgecut;

	WgtInfo wgtInfo;
public:
	PartGraph(){
		nvtxs = nedges = 0;
		xadj = NULL;
		adjncy = NULL;
		vwgt = NULL;
		ewgt = NULL;
		cewgt = NULL;
		adjwgt = NULL;

		tolerance = 0;
		totalvwgt = 0;
		currWgtX = 0;
		edgecut = 0;
	};
	PartGraph(int nvtxs,int* xadj,int* adjncy,int* vwgt,int* ewgt,int* cewgt,int* adjwgt){
		this->nvtxs = nvtxs;
		this->xadj = xadj;
		this->adjncy = adjncy;
		this->vwgt = vwgt;
		this->ewgt = ewgt;
		this->cewgt = cewgt;
		this->adjwgt = adjwgt;
		nedges = xadj[nvtxs];

		tolerance = 0;
		totalvwgt = 0;
		for(int i = 0; i < nvtxs; i++) totalvwgt += vwgt[i];
		currWgtX = 0;
		edgecut = 0;
	};
	~PartGraph(){ };
	void DeleteGraph(){
		free(adjwgt);
		free(cewgt);
		free(ewgt);
		free(vwgt);
		free(adjncy);
		free(xadj);
	};
	int Vsize(){
		return nvtxs;
	};
	int Esize(){
		return nedges;
	}
	int Vwgt(int v){
		return vwgt[v];
	}
	int Ewgt(int i){
		return ewgt[i];
	}
	int Xadj(int v){
		return xadj[v];
	}
	int Adjncy(int i){
		return adjncy[i];
	}
	int Tolerance(){
		return tolerance;
	}
	int UpperEdgeGain();
	int UpperVertGain();
	int UpperVertGainCewgt();
	int MaxVertexWeight();
	void SetTolerance(ndOptions* options);

	void Show();
	void Show(int* partition);

	// partition[nvtxs] // 0->X, 1->Y // ratioY = 1.0 - ratioX
	int Partition2(ndOptions* options, double ratioX, int* partition); 

	// partition[nvtxs] // 0->X, 1->Y, 2->S // ratioY = 1.0 - ratioX
	int Partition3(ndOptions* options, double ratioX, int* partition);

	int GetCurrWgt(int* partition);

	int VertSepFromEdgeSep(int* partition);
	int GetBoundary(std::vector<int>& bvs,std::map<int,int>& ibvs,int* partition);
	bool VertSepIsOK(int* partition);

	int GetEdgeGain(int v, int* partition);
	int GetEdgecut(int* partition);

	int GetVertGain(int v, int to, int* partition); // gain when separator vertex "v" is moved to partition "to"
	int GetVertSepSize(int* partition);
	int GetVertGainCewgt(int v, int to, int* partition); 
	int GetVertSepSizeCewgt(int* partition);

	int Coarsening(int* match, int* map); // map[nvtxs], match[nvtxs], ret val = nvtxs of coarser graph
	int RandomMatching(int* match);
	int HeavyEdgeMatching(int* match);
	int LiteEdgeMatching(int* match);
	int Mapping(const int* match, int* map); // ret val = nvtxs of coarser graph
	int GenerateCoarserGraph(int newSize, const int* match, const int* map, PartGraph* newGraph);
	
	int InitPartitioningEdge(double ratioX, int* partition); // divide G to X=0 and Y=1. ret val = |Y|
	int GGPartitioningEdge(double ratioX, int* partition); // Graph Growing Algorithm
	int GGGPartitioningEdge(double ratioX, int* partition); // Greedy Graph Growing Algorithm
	int GetLargeGainVertexFromBoundary(std::vector<int> &list,int* partition,int val);

	int InitPartitioningVert(double ratioX, int* partition);

	int UncoarseningEdge(double ratioX, int* map,int* coarserPart,int* partition);
	int RefineEdge(double ratioX, int* partition);

	int UncoarseningVert(double ratioX, int* map, int* coarserPart, int* partition);
	int RefineVert(int* partition); // using this->wgtInfo internally
};

#endif
