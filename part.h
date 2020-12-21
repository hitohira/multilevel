#ifndef MY_PART_H
#define MY_PART_H

#include <stdlib.h>
#include <limits.h>
#include <vector>

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
	int* cewgt; // the wgt of edges that have been contacted to create v, size = nvtxs
	int* adjwgt; // the sum of wgt of the edges adjcent to v,  size = nvtxs

	int tolerance;
public:
	PartGraph(){
		nvtxs = nedges = 0;
		xadj = NULL;
		adjncy = NULL;
		vwgt = NULL;
		ewgt = NULL;
		cewgt = NULL;
		adjwgt = NULL;

		tolerance = INT_MAX;
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

		tolerance = INT_MAX;
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
	int MaxEdgeGain();
	int MaxVertexWeight();
	void SetTolerance(ndOptions* options);

	void Show();
	void Show(int* partition);

	int Partition2(ndOptions* options, double ratioX, int* partition); 
	// partition[nvtxs] // 0->X, 1->Y, 2->S // ratioY = 1.0 - ratioX

	int GetEdgeGain(int v, int* partition);
	int GetEdgecut(int* partition);

	int Coarsening(int* match, int* map); // map[nvtxs], match[nvtxs], ret val = nvtxs of coarser graph
	int RandomMatching(int* match);
	int HeavyEdgeMatching(int* match);
	int Mapping(const int* match, int* map); // ret val = nvtxs of coarser graph
	int GenerateCoarserGraph(int newSize, const int* match, const int* map, PartGraph* newGraph);
	
	int InitPartitioning(double ratioX, int* partition); // divide G to X=0 and Y=1. ret val = |Y|
	int GGPartitioningEdge(double ratioX, int* partition); // Graph Growing Algorithm
	int GGGPartitioningEdge(double ratioX, int* partition); // Greedy Graph Growing Algorithm
	int GetLargeGainVertexFromBoundary(std::vector<int> &list,int* partition,int val);

	int Uncoarsening(int* map,int* coarserPart,int* partition);
};

#endif
