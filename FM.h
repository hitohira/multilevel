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
};

class DoubleLinkedList {
private:
	DoubleLinkedListItem* dummy;
public:
	DoubleLinkedList();
	~DoubleLinkedList();
	DoubleLinkedListItem* PushFront(int iitem);
	DoubleLinkedListItem* Insert(DoubleLinkedListItem* pitem, int iitem);
	DoubleLinkedListItem* Erase(DoubleLinkedListItem* pitem);
	DoubleLinkedListItem* Begin();
};

class FMDATA {
private:
	DoubleLinkedList* bucket[2];
	DoubleLinkedListItem** cell;
	int nvtxs;
	int pmax;
	int max_gain[2];
	std::vector<int> free_cell;
	int end_update_idx;
public:
	FMDATA(PartGraph* pg, int* partition);
	~FMDATA();
	void Construct(PartGraph* pg, int* partition);
	void Reconstruct(PartGraph* pg, int* partition); // reset free_cell -> reconstruct buckets

	int IndexDLList(int gain); // -pmax <= gain <= pmax -> 0 <= ret <= 2*pmax
};

#endif
