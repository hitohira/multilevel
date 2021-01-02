#ifndef FM_STRUCTURE_H
#define FM_STRUCTURE_H

#include <stdio.h>
#include <stdlib.h>
#include <vector>

#include "part.h"
class PartGraph;

class DoubleLinkedListItem {
public:
	DoubleLinkedListItem* pre;
	DoubleLinkedListItem* next;
	int item;

	DoubleLinkedListItem();
	~DoubleLinkedListItem();
	DoubleLinkedListItem(int item);
	DoubleLinkedListItem* Insert(int iitem);
	DoubleLinkedListItem* InsNext(int iitem);
	DoubleLinkedListItem* GetFront();
};

class DoubleLinkedList {
private:
	DoubleLinkedListItem* dummy;
public:
	DoubleLinkedList();
	~DoubleLinkedList();
	void SetState(int state_item);
	DoubleLinkedListItem* PushFront(int iitem);
	DoubleLinkedListItem* Insert(DoubleLinkedListItem* pitem, int iitem);
	DoubleLinkedListItem* Erase(DoubleLinkedListItem* pitem);
	DoubleLinkedListItem* Begin();
	void EraseAll();
};

typedef struct VGPair{
	int v, gain;
} VGPair;

class FMDATA {
private:
	DoubleLinkedList* bucket[2];
	DoubleLinkedListItem** cell;
	int* cell_gain;
	int nvtxs;
	int pmax;
	int max_gain[2];
	std::vector<int> free_cell;
	int end_update_idx;
	int min_edgecut;
	int edgecut;

	double st_create;
	double st_pop;
	double st_update;
	double st_backtrack;
	double st_option;
public:
	FMDATA(PartGraph* pg, int* partition);
	~FMDATA();
	void Construct(PartGraph* pg, int* partition);
	void Reconstruct(PartGraph* pg, int* partition); // reset free_cell -> reconstruct buckets
	void ResetSt();
	int Edgecut(){
		return edgecut;
	}

	int RefineEdge(int currWgtX,int minWgtX, int maxWgtX, PartGraph* pg, int* partition); // ret currWgtX
	int RefineEdgeInner(int currWgtX,int minWgtX, int maxWgtX, PartGraph* pg, int* partition);

	VGPair PopMovedVertexEdge(int currWgtX, int minWgtX,int maxWgtX, PartGraph* pg);
	int UpdatePartEdgecutAndGainEdge(VGPair vg, int currWgtX, PartGraph* pg, int* partition);// ret newWgtX
	int BackTrackPartitionEdge(int currWgtX, PartGraph* pg, int* partition); // ret newWgtX

	int IndexBucket(int gain); // -pmax <= gain <= pmax -> 0 <= ret <= 2*pmax
	void PrintSt();
};
/*+
class FMDATAvert {
private:
	DoubleLinkedList* bucket;
	DoubleLinkedListItem** cell;
	int* cell_gain;
	int nvtxs;
	int pmax;
	int max_gain;
	std::vector<int> free_cell;
	int end_update_idx;
	int min_score;
	int score;

	double st_create;
	double st_pop;
	double st_update;
	double st_backtrack;
	double st_option;
public:
	FMDATAvert(PartGraph* pg, int* partition);
	~FMDATAvert();
	void Construct(PartGraph* pg, int* partition);
	void Reconstruct(PartGraph* pg, int* partition); // reset free_cell -> reconstruct buckets
	void ResetSt();
	int Score(){
		return score;
	}

	int RefineVert(int currWgtX,int minWgtX, int maxWgtX, PartGraph* pg, int* partition); // ret currWgtX
	int RefineVertInner(int currWgtX,int minWgtX, int maxWgtX, PartGraph* pg, int* partition);

	VGPair PopMovedVertexVert(int currWgtX, int minWgtX,int maxWgtX, PartGraph* pg);
	int UpdatePartScoreAndGainVert(VGPair vg, int currWgtX, PartGraph* pg, int* partition);// ret newWgtX
	int BackTrackPartitionVert(int currWgtX, PartGraph* pg, int* partition); // ret newWgtX

	int IndexBucket(int gain); // -pmax <= gain <= pmax -> 0 <= ret <= 2*pmax
	void PrintSt();
};
*/
	
#endif
