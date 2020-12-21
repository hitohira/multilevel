#include <stdio.h>
#include <stdlib.h>

#include "FM.h"

DoubleLinkedListItem::DoubleLinkedListItem(){
	pre = NULL;
	next = NULL;
	item = -1;
}
DoubleLinkedListItem::~DoubleLinkedListItem(){
	if(next) next->pre = pre;
	if(pre) pre->next = next;
}
DoubleLinkedListItem::DoubleLinkedListItem(int item){
	DoubleLinkedListItem();
	this->item = item;
}
DoubleLinkedListItem* DoubleLinkedListItem::Insert(int iitem){
	DoubleLinkedListItem* obj = new DoubleLinkedListItem();
	if(obj == NULL) return NULL;
	obj->next = this;
	obj->pre = this->pre;
	obj->item = iitem;
	if(this->pre) this->pre->next = obj;
	this->pre = obj;
	return obj;
}
DoubleLinkedListItem* DoubleLinkedListItem::InsNext(int iitem){
	DoubleLinkedListItem* obj = new DoubleLinkedListItem();
	obj->next = this->next;
	obj->pre = this;
	obj->item = iitem;
	if(this->next) this->next->pre = obj;
	this->next = obj;
	return obj;
}



DoubleLinkedList::DoubleLinkedList(){
	dummy = new DoubleLinkedListItem();
}
DoubleLinkedList::~DoubleLinkedList(){
	DoubleLinkedListItem* itr = dummy;
	while(itr != NULL){
		itr = Erase(itr);
	}
}
DoubleLinkedListItem* DoubleLinkedList::PushFront(int iitem){
	return dummy->InsNext(iitem);
}
DoubleLinkedListItem* DoubleLinkedList::Insert(DoubleLinkedListItem* pitem, int iitem){
	if(pitem == NULL) return NULL;
	return pitem->Insert(iitem);
}
DoubleLinkedListItem* DoubleLinkedList::Erase(DoubleLinkedListItem* pitem){
	if(pitem == NULL) return NULL;
	DoubleLinkedListItem* ret = pitem->next;
	delete pitem;
	return ret;
}
DoubleLinkedListItem* DoubleLinkedList::Begin(){
	return dummy->next;
}



FMDATA::FMDATA(PartGraph* pg, int* partition){
	Construct(pg,partition);
}
void FMDATA::Construct(PartGraph* pg,int* partition){
	free_cell = std::vector<int>();
	end_update_idx = 0;
	pmax = pg->MaxEdgeGain();
	nvtxs = pg->Vsize();
	max_gain[0] = max_gain[1] = -pmax;
	bucket[0] = new DoubleLinkedList[2*pmax+1];
	bucket[1] = new DoubleLinkedList[2*pmax+1];
	cell = new DoubleLinkedListItem*[nvtxs];

	for(int i = 0; i < nvtxs; i++){
		int c = partition[i];
		int gain = pg->GetEdgeGain(i,partition);
		int idx = IndexDLList(gain);
		cell[i] = bucket[c][idx].PushFront(gain);
		if(gain > max_gain[c]){
			max_gain[c] = gain;
		}
	}
}
FMDATA::~FMDATA(){
	delete[] bucket[0];
	delete[] bucket[1];
	delete[] cell;
}
void FMDATA::Reconstruct(PartGraph* pg, int* partition){
	for(int i = 0; i < (int)free_cell.size(); i++){
		int v = free_cell[i];
		int c = partition[v];
		int gain = pg->GetEdgeGain(v,partition);
		int idx = IndexDLList(gain);
		cell[v] = bucket[c][idx].PushFront(gain);
		if(gain > max_gain[c]){
			max_gain[c] = gain;
		}
	}
	std::vector<int>().swap(free_cell);
	end_update_idx = 0;
}


int FMDATA::IndexDLList(int gain){
	return gain + pmax;
}
