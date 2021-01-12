#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <sys/time.h>

#include "FM.h"

static double GetTime(){
#ifdef FM_DEBUG
	struct timeval tv;
	gettimeofday(&tv,NULL);
	return tv.tv_sec + tv.tv_usec*1e-6;
#endif
	return 0.0;
}

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
DoubleLinkedListItem* DoubleLinkedListItem::GetFront(){
	if(this->pre == NULL) return this;
	return this->pre->GetFront();
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
void DoubleLinkedList::SetState(int state_item){
	dummy->item = state_item;
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
void DoubleLinkedList::EraseAll(){
	DoubleLinkedListItem* ptr = dummy->next;
	while(ptr != NULL){
		ptr = Erase(ptr);
	}
}



FMDATA::FMDATA(PartGraph* pg, int* partition){
	ResetSt();
	double t1 = GetTime();
	Construct(pg,partition);
	double t2 = GetTime();
	st_create += t2-t1;
}
void FMDATA::Construct(PartGraph* pg,int* partition){
	free_cell = std::vector<int>();
	end_update_idx = 0;
	edgecut = pg->GetEdgecut(partition);
	min_edgecut = edgecut;
	pmax = pg->UpperEdgeGain();
	nvtxs = pg->Vsize();
	max_gain[0] = max_gain[1] = -pmax;
	bucket[0] = new DoubleLinkedList[2*pmax+1];
	bucket[1] = new DoubleLinkedList[2*pmax+1];
	cell = new DoubleLinkedListItem*[nvtxs];
	cell_gain = new int[nvtxs];

	for(int i = -pmax; i <= pmax; i++){
		int idx = IndexBucket(i);
		bucket[0][idx].SetState(i);
		bucket[1][idx].SetState(i);
	}
	for(int i = 0; i < nvtxs; i++){
		int c = partition[i];
		int gain = pg->GetEdgeGain(i,partition);
		int idx = IndexBucket(gain);
		cell[i] = bucket[c][idx].PushFront(i);
		cell_gain[i] = gain;
		if(gain > max_gain[c]){
			max_gain[c] = gain;
		}
	}
}
FMDATA::~FMDATA(){
	delete[] bucket[0];
	delete[] bucket[1];
	delete[] cell;
	delete[] cell_gain;
}
void FMDATA::Reconstruct(PartGraph* pg, int* partition){
	std::vector<int>().swap(free_cell);
	end_update_idx = 0;
	max_gain[0] = max_gain[1] = -pmax;

	for(int i = -pmax; i <= pmax; i++){
		int idx = IndexBucket(i);
		bucket[0][idx].EraseAll();
		bucket[1][idx].EraseAll();
	}
	for(int i = 0; i < nvtxs; i++){
		int c = partition[i];
		int gain = pg->GetEdgeGain(i,partition);
		int idx = IndexBucket(gain);
		cell[i] = bucket[c][idx].PushFront(i);
		cell_gain[i] = gain;
		if(gain > max_gain[c]){
			max_gain[c] = gain;
		}
	}
}
void FMDATA::ResetSt(){
	st_create = st_pop = st_update = st_backtrack = st_option = 0.0;
}
void FMDATA::PrintSt(){
#ifdef FM_DEBUG
	fprintf(stderr,"FMDATA status\ncreate: %f\npop: %f\nupdate: %f\nbacktrack: %f\noption: %f\n",st_create,st_pop,st_update,st_backtrack,st_option);
#endif
}
static bool ConditionWgt(int newWgtX,int minWgtX, int maxWgtX,int oldWgtX, int part){
	return (minWgtX <= newWgtX && newWgtX <= maxWgtX) ||
	(newWgtX < minWgtX && oldWgtX < minWgtX && newWgtX > oldWgtX) ||
	(newWgtX > maxWgtX && oldWgtX > maxWgtX && newWgtX < oldWgtX);
/*
	static int tms = 0;
	tms++;
	if(minWgtX <= newWgtX && newWgtX <= maxWgtX) return true;
	if(newWgtX < minWgtX && oldWgtX < minWgtX && newWgtX > oldWgtX){
		if(tms < 20) printf("%d %d %d\n",oldWgtX,newWgtX,minWgtX);
		return false;
	}
	if(newWgtX > maxWgtX && oldWgtX > maxWgtX && newWgtX < oldWgtX) return true;
	return false;
*/
}
static bool IsBalancedWgt(int newWgtX,int minWgtX, int maxWgtX){
	return (minWgtX <= newWgtX && newWgtX <= maxWgtX);
}
VGPair FMDATA::PopMovedVertexEdge(int currWgtX, int minWgtX,int maxWgtX, PartGraph* pg){
	int lgb = max_gain[0] > max_gain[1] ? 0 : 1;
	int smb = 1 - lgb;
	int threshold = -pmax;
	for(int i = max_gain[lgb]; i >= threshold; i--){
		int idx = IndexBucket(i);
		if(max_gain[smb] < i){ // not consider bucket[smb]
			DoubleLinkedListItem* ptr = bucket[lgb][idx].Begin();
			while(ptr != NULL){
				int v = ptr->item;
				int w = pg->Vwgt(v);
				int newWgtX = lgb == 0 ? currWgtX - w : currWgtX + w;
				if(ConditionWgt(newWgtX,minWgtX,maxWgtX,currWgtX,lgb)){
					bucket[lgb][idx].Erase(ptr);
					cell[v] = NULL;
					free_cell.push_back(v);
					VGPair pr; pr.v = v; pr.gain = i;		
//					fprintf(stderr,"vgm %d %d",pr.v,pr.gain);
					return pr;
				}
				ptr = ptr->next;
			}
		}
		else{ // consider both buckets
			int fib,seb;
			DoubleLinkedListItem* fi, * se;
			if(abs(currWgtX-minWgtX) < abs(currWgtX-maxWgtX)){ // 1 -> 0 prior
				fib = 1;
				seb = 0;
			}
			else{ // 0 -> 1 prior
				fib = 0;
				seb = 1;
			}
			fi = bucket[fib][idx].Begin();
			se = bucket[seb][idx].Begin();
			while(fi != NULL){
				int v = fi->item;
				int w = pg->Vwgt(v);
				int newWgtX = fib == 0 ? currWgtX - w : currWgtX + w;
				if(ConditionWgt(newWgtX,minWgtX,maxWgtX,currWgtX,fib)){
					bucket[fib][idx].Erase(fi);
					cell[v] = NULL;
					free_cell.push_back(v);
					VGPair pr; pr.v = v; pr.gain = i;		
//					fprintf(stderr,"vgf %d %d",pr.v,pr.gain);
					return pr;
				}
				fi = fi->next;
			}
			while(se != NULL){
				int v = se->item;
				int w = pg->Vwgt(v);
				int newWgtX = seb == 0 ? currWgtX - w : currWgtX + w;
				if(ConditionWgt(newWgtX,minWgtX,maxWgtX,currWgtX,seb)){
					bucket[seb][idx].Erase(se);
					cell[v] = NULL;
					free_cell.push_back(v);
					VGPair pr; pr.v = v; pr.gain = i;		
//					fprintf(stderr,"vgs %d %d",pr.v,pr.gain);
					return pr;
				}
				se = se->next;
			}
		}
	}
	VGPair vgpair; vgpair.v = -1; vgpair.gain = 0;
	return vgpair;
}

int FMDATA::UpdatePartEdgecutAndGainEdge(VGPair vg, int currWgtX, int minWgtX, int maxWgtX, PartGraph* pg, int* partition){
//	fprintf(stderr,"enter\n");
	int movedV = vg.v;
	if(movedV == -1) return currWgtX;
//	fprintf(stderr,"\t%d gain  %d == %d\n",movedV,vg.gain,pg->GetEdgeGain(movedV,partition));
	int sgn = partition[movedV] == 0 ? -1 : 1;
	int newWgtX = currWgtX + sgn*pg->Vwgt(movedV);

//	fprintf(stderr,"adj gain\n");
	for(int i = pg->Xadj(movedV); i < pg->Xadj(movedV+1); i++){
		int u = pg->Adjncy(i);
		int ew = pg->Ewgt(i);
		DoubleLinkedListItem* ptr = cell[u];
		if(ptr == NULL) continue;
double t1 = GetTime();
//		int gain = ptr->GetFront()->item;
		int gain = cell_gain[u];
double t2 = GetTime();
st_option += t2-t1;
		delete ptr;
		if(partition[u] == partition[movedV]){
			gain += 2*ew;
		}
		else{
			gain -= 2*ew;
		}
		int partu = partition[u];
		int idx = IndexBucket(gain);
		cell[u] = bucket[partu][idx].PushFront(u);
		cell_gain[u] = gain;
		if(gain > max_gain[partu]) max_gain[partu] = gain;
	}
	
//	fprintf(stderr,"max_gain\n");
	for(int i = max_gain[0]; i > -pmax; i--){
		int idx = IndexBucket(i);
		if(bucket[0][idx].Begin() != NULL){
			max_gain[0] = i;
			break;
		}
//		fprintf(stderr,"idx %d / %d\n",i,-pmax);
	}
	for(int i = max_gain[1]; i > -pmax; i--){
		int idx = IndexBucket(i);
		if(bucket[1][idx].Begin() != NULL){
			max_gain[1] = i;
			break;
		}
	}
	
	edgecut -= vg.gain;
	if(!IsBalancedWgt(newWgtX,minWgtX,maxWgtX)){ // TODO consider when nvtxs is close to original nvtxs
		min_edgecut = edgecut;
		end_update_idx = free_cell.size();
	}
	if(edgecut < min_edgecut){
		min_edgecut = edgecut;
		end_update_idx = free_cell.size();
	}
//	fprintf(stderr,"partition\n");
	partition[movedV] = 1 - partition[movedV];
//	fprintf(stderr,"exit\n");
	return newWgtX;
}

int FMDATA::BackTrackPartitionEdge(int currWgtX, PartGraph* pg, int* partition){
	for(int i = end_update_idx; i < (int)free_cell.size(); i++){
		int v = free_cell[i];
		partition[v] = 1 - partition[v];
		int sgn = partition[v] == 0 ? 1 : -1;
		currWgtX += sgn * pg->Vwgt(v);
	}
	end_update_idx = free_cell.size();
	edgecut = min_edgecut;
	return currWgtX;
}

int FMDATA::IndexBucket(int gain){
	return gain + pmax;
}


int FMDATA::RefineEdge(int currWgtX,int minWgtX, int maxWgtX, PartGraph* pg, int* partition){
	int old_edgecut;
	while(1){
//		fprintf(stderr,"edgecut : wgt = %d : %d\n",edgecut,currWgtX);
//		fprintf(stderr,"min max %d  %d\n",minWgtX,maxWgtX);
		old_edgecut = edgecut;
		currWgtX = RefineEdgeInner(currWgtX,minWgtX,maxWgtX,pg,partition);
		if(old_edgecut == edgecut) break;
//		if(old_edgecut < edgecut){ // it happen when imbalance
//			fprintf(stderr,"old_edgecut < edgecut!!!\n");
//			return -1;
//		}
//		fprintf(stderr,"recunstruct\n");
		double t1 = GetTime();
		Reconstruct(pg,partition);
		double t2 = GetTime();
		st_create += t2-t1;
	}
	return currWgtX;
}

int FMDATA::RefineEdgeInner(int currWgtX,int minWgtX, int maxWgtX, PartGraph* pg, int* partition){
	for(int rep=0;;rep++){
//		fprintf(stderr,"score %d == %d, min %d\n",edgecut,pg->GetEdgecut(partition),min_edgecut);
//		assert(edgecut == pg->GetEdgecut(partition));
//		fprintf(stderr,"%dpop",rep);
		double t1 = GetTime();
		VGPair vg = PopMovedVertexEdge(currWgtX,minWgtX,maxWgtX,pg);
		double t2 = GetTime();
		st_pop += t2-t1;
		if(vg.v == -1) break;
//		fprintf(stderr," update\n");
		currWgtX = UpdatePartEdgecutAndGainEdge(vg,currWgtX,minWgtX,maxWgtX,pg,partition);
//		fprintf(stderr,"%d %d %d\n",currWgtX, pg->GetCurrWgt(partition), pg->GetEdgecut(partition));
		double t3 = GetTime();
		st_update += t3-t2;
	}
//	fprintf(stderr," backtrack\n");
	double t4 = GetTime();
	currWgtX = BackTrackPartitionEdge(currWgtX,pg,partition);
	double t5 = GetTime();
	st_backtrack += t5-t4;
//	fprintf(stderr,"score end inner %d == %d\n",edgecut,pg->GetEdgecut(partition));
	return currWgtX;
}
