#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <sys/time.h>
#include <set>

#include "FM.h"

static double GetTime(){
#ifdef FM_DEBUG
	struct timeval tv;
	gettimeofday(&tv,NULL);
	return tv.tv_sec + tv.tv_usec*1e-6;
#endif
	return 0.0;
}


FMDATAvert2::FMDATAvert2(PartGraph* pg, int* partition){
	ResetSt();
	double t1 = GetTime();
	Construct(pg,partition);
	double t2 = GetTime();
	st_create += t2-t1;
}
void FMDATAvert2::Construct(PartGraph* pg,int* partition){
	free_cell = std::vector<Trace>();
	end_update_idx = 0;
	score = pg->GetVertSepSizeCewgt(partition);
	min_score = score;
	pmax = pg->UpperVertGainCewgt();
	nvtxs = pg->Vsize();
	max_gain[0] = max_gain[1] = -pmax;
	bucket[0] = new DoubleLinkedList[2*pmax+1];
	bucket[1] = new DoubleLinkedList[2*pmax+1];
	cell[0] = new DoubleLinkedListItem*[nvtxs];
	cell[1] = new DoubleLinkedListItem*[nvtxs];
	cell_gain[0] = new int[nvtxs];
	cell_gain[1] = new int[nvtxs];

	for(int i = -pmax; i <= pmax; i++){
		int idx = IndexBucket(i);
		bucket[0][idx].SetState(i);
		bucket[1][idx].SetState(i);
	}
	for(int i = 0; i < nvtxs; i++){
		if(partition[i] != 2) continue;
		for(int j = 0; j < 2; j++){
			int gain = pg->GetVertGainCewgt(i,j,partition);
			int idx = IndexBucket(gain);
			cell[j][i] = bucket[j][idx].PushFront(i);
			cell_gain[j][i] = gain;
			if(gain > max_gain[j]){
				max_gain[j] = gain;
			}
		}
	}
}
FMDATAvert2::~FMDATAvert2(){
	delete[] bucket[0];
	delete[] bucket[1];
	delete[] cell[0];
	delete[] cell[1];
	delete[] cell_gain[0];
	delete[] cell_gain[1];
}
void FMDATAvert2::Reconstruct(PartGraph* pg, int* partition){
	std::vector<Trace>().swap(free_cell);
	end_update_idx = 0;
	max_gain[0] = max_gain[1] = -pmax;

	for(int i = -pmax; i <= pmax; i++){
		int idx = IndexBucket(i);
		bucket[0][idx].EraseAll();
		bucket[1][idx].EraseAll();
	}
	for(int i = 0; i < nvtxs; i++){
		if(partition[i] != 2) continue;
		for(int j = 0; j < 2; j++){
			int gain = pg->GetVertGainCewgt(i,j,partition);
			int idx = IndexBucket(gain);
			cell[j][i] = bucket[j][idx].PushFront(i);
			cell_gain[j][i] = gain;
			if(gain > max_gain[j]){
				max_gain[j] = gain;
			}
		}
	}
}
void FMDATAvert2::ResetSt(){
	st_create = st_pop = st_update = st_backtrack = st_option = 0.0;
}
void FMDATAvert2::PrintSt(){
#ifdef FM_DEBUG
	fprintf(stderr,"FMDATAvert2 status\ncreate: %f\npop: %f\nupdate: %f\nbacktrack: %f\noption: %f\n",st_create,st_pop,st_update,st_backtrack,st_option);
#endif
}
static double CalcRatio(WgtInfo* wgtInfo){
	return 1.0*(wgtInfo->wgt[0] + wgtInfo->wgt[2]) / (wgtInfo->wgt[0] + wgtInfo->wgt[1] + 2*wgtInfo->wgt[2]);
}
//TODO consider vertices that come to separator
static double CalcRatioNew(int part, int wgt, WgtInfo* wgtInfo){
	double sep = wgtInfo->wgt[2];
	if(part == 0){
		return 1.0*(wgtInfo->wgt[0]+wgt + sep) / (wgtInfo->wgt[0]+wgt + wgtInfo->wgt[1] + 2*sep);
	}
	else {
		return 1.0*(wgtInfo->wgt[0] + sep) / (wgtInfo->wgt[0] + wgtInfo->wgt[1]+wgt + 2*sep);
	}
}
static bool ConditionWgt(double ratioNow, double ratioNew, WgtInfo* wgtInfo){
	double ratioX = wgtInfo->ratioX;
	double tol = wgtInfo->tol;
	double ratioL = ratioX - tol;
	double ratioH = ratioX + tol;
	return (ratioL <= ratioNew && ratioNew <= ratioH) ||
		(ratioNew < ratioL && ratioNow < ratioNew) ||
		(ratioNew > ratioH && ratioNow > ratioNew);
}
VGTri FMDATAvert2::PopMovedVertexVert(WgtInfo* wgtInfo, PartGraph* pg){
	int lgb = max_gain[0] > max_gain[1] ? 0 : 1;
	int smb = 1 - lgb;
	int threshold = -pmax;
	double ratioNow = CalcRatio(wgtInfo);
	for(int i = max_gain[lgb]; i >= threshold; i--){
		int idx = IndexBucket(i);
		if(max_gain[smb] < i){ // not consider bucket[smb]
			DoubleLinkedListItem* ptr = bucket[lgb][idx].Begin();
			while(ptr != NULL){
				int v = ptr->item;
				int w = pg->Vwgt(v);
				double ratioNew = CalcRatioNew(lgb,w,wgtInfo);
				if(ConditionWgt(ratioNow,ratioNew,wgtInfo)){
					bucket[lgb][idx].Erase(ptr);
					bucket[smb][idx].Erase(cell[smb][v]);
					cell[0][v] = &moved;
					cell[1][v] = &moved;
					VGTri pr; pr.v = v; pr.gain = i;		
					pr.to = lgb;
//					fprintf(stderr,"vgm %d %d",pr.v,pr.gain);
					return pr;
				}
				ptr = ptr->next;
			}
		}
		else{ // consider both buckets
			int fib,seb;
			DoubleLinkedListItem* fi, * se;
			if(ratioNow < wgtInfo->ratioX){ // 1 -> 0 prior
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
				double ratioNew = CalcRatioNew(fib,w,wgtInfo);
				if(ConditionWgt(ratioNow,ratioNew,wgtInfo)){
					bucket[fib][idx].Erase(fi);
					bucket[seb][idx].Erase(cell[seb][v]);
					cell[0][v] = &moved;
					cell[1][v] = &moved;
					VGTri pr; pr.v = v; pr.gain = i;		
					pr.to = fib;
//					fprintf(stderr,"vgf %d %d",pr.v,pr.gain);
					return pr;
				}
				fi = fi->next;
			}
			while(se != NULL){
				int v = se->item;
				int w = pg->Vwgt(v);
				double ratioNew = CalcRatioNew(seb,w,wgtInfo);
				if(ConditionWgt(ratioNow,ratioNew,wgtInfo)){
					bucket[seb][idx].Erase(se);
					bucket[fib][idx].Erase(cell[fib][v]);
					cell[0][v] = &moved;
					cell[1][v] = &moved;
					VGTri pr; pr.v = v; pr.gain = i;		
					pr.to = seb;
//					fprintf(stderr,"vgs %d %d",pr.v,pr.gain);
					return pr;
				}
				se = se->next;
			}
		}
	}
	VGTri vgt; vgt.v = -1; vgt.to = -1; vgt.gain = 0;
	return vgt;
}

int FMDATAvert2::UpdatePartScoreAndGainVert(VGTri vg, WgtInfo* wgtInfo, PartGraph* pg, int* partition){
//	fprintf(stderr,"enter\n");
	int movedV = vg.v;
	if(movedV == -1) return score;

	assert(pg->GetVertGainCewgt(movedV,0,partition) == cell_gain[0][movedV]);
	assert(pg->GetVertGainCewgt(movedV,1,partition) == cell_gain[1][movedV]);
	

	Trace trace; trace.v = movedV; trace.from = 2;
	free_cell.push_back(trace);
	partition[movedV] = vg.to;
	wgtInfo->wgt[vg.to] += pg->Vwgt(movedV);
	wgtInfo->wgt[2] -= pg->Vwgt(movedV);
	
	double t1 = GetTime();
//	fprintf(stderr,"update list\n");
	std::set<int> updateList;
	for(int i = pg->Xadj(movedV); i < pg->Xadj(movedV+1); i++){
		int u = pg->Adjncy(i);
		if(partition[u] != vg.to){
			int ppu = partition[u];
			if(partition[u] != 2){
				Trace trace; trace.v = u; trace.from = partition[u];
				free_cell.push_back(trace);
				partition[u] = 2;
				wgtInfo->wgt[2] += pg->Vwgt(u);
				wgtInfo->wgt[ppu] -= pg->Vwgt(u);
				if(cell[0][u] != &moved){
					cell[0][u] = NULL;
					cell[1][u] = NULL;
				}
				for(int k = pg->Xadj(u); k < pg->Xadj(u+1); k++){
					int u2 = pg->Adjncy(k); // 2 hop adjacent of movedV
					if(partition[u2] == 2 && cell[0][u2] != &moved){ // if u2 is in S, u2's gain is also changed
						updateList.insert(u2);
					}
				}
			}
			if(cell[0][u] != &moved){ // because this has not already moved and shold add bucket
				updateList.insert(u);
			}
		}
	}
//	fprintf(stderr,"update gain\n");
	for(std::set<int>::iterator itr = updateList.begin(); itr != updateList.end(); itr++){
		int u = *itr;
		if(cell[0][u] != NULL){ // vertex in S before move -> delete bucket item
			delete cell[0][u];
			delete cell[1][u];
		}
		for(int j = 0; j < 2; j++){ // re-calcuration of gain and register to bucket
			int gain = pg->GetVertGainCewgt(u,j,partition);
			int idx = IndexBucket(gain);
			cell[j][u] = bucket[j][idx].PushFront(u);
			cell_gain[j][u] = gain;
			if(gain > max_gain[j]) max_gain[j] = gain;
		}
	}
	double t2 = GetTime();
	st_option += t2-t1;
	
//	fprintf(stderr,"update score\n");
	score -= vg.gain;
	if(score < min_score){
		min_score = score;
		end_update_idx = free_cell.size();
	}
//	assert(score == pg->GetVertSepSize(partition));
	
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
	
//	fprintf(stderr,"exit\n");
	return score;
}

int FMDATAvert2::BackTrackPartitionVert(WgtInfo* wgtInfo, PartGraph* pg, int* partition){
	for(int i = (int)free_cell.size()-1 ; i >= end_update_idx; i--){
		int v = free_cell[i].v;
		int ppv = partition[v];
		partition[v] = free_cell[i].from;
		wgtInfo->wgt[ppv] -= pg->Vwgt(v);
		wgtInfo->wgt[partition[v]] += pg->Vwgt(v);
	}
	end_update_idx = free_cell.size();
	score = min_score;
	return score;
}

int FMDATAvert2::IndexBucket(int gain){
	return gain + pmax;
}


int FMDATAvert2::RefineVert(WgtInfo* wgtInfo, PartGraph* pg, int* partition){
	int old_score;
	while(1){
//		fprintf(stderr,"edgecut : wgt = %d : %d\n",edgecut,currWgtX);
//		fprintf(stderr,"min max %d  %d\n",minWgtX,maxWgtX);
		old_score = score;
		score = RefineVertInner(wgtInfo,pg,partition);
		if(old_score == score) break;
		if(old_score < score){
//			fprintf(stderr,"old_edgecut < edgecut!!!\n");
			return -1;
		}
//		fprintf(stderr,"recunstruct\n");
		double t1 = GetTime();
		Reconstruct(pg,partition);
		double t2 = GetTime();
		st_create += t2-t1;
	}
	return score;
}

int FMDATAvert2::RefineVertInner(WgtInfo* wgtInfo, PartGraph* pg, int* partition){
	for(int rep=0;;rep++){
//		fprintf(stderr,"score %d == %d, min %d\n",edgecut,pg->GetEdgecut(partition),min_edgecut);
//		assert(edgecut == pg->GetEdgecut(partition));
//		fprintf(stderr,"%dpop",rep);
		double t1 = GetTime();
		VGTri vg = PopMovedVertexVert(wgtInfo,pg);
		double t2 = GetTime();
		st_pop += t2-t1;
		if(vg.v == -1) break;
//	fprintf(stderr," begin update\n");
		UpdatePartScoreAndGainVert(vg,wgtInfo,pg,partition);
//	fprintf(stderr," end update\n");
//		fprintf(stderr,"%d %d %d\n",currWgtX, pg->GetCurrWgt(partition), pg->GetEdgecut(partition));
		double t3 = GetTime();
		st_update += t3-t2;
	}
//	fprintf(stderr," backtrack\n");
	double t4 = GetTime();
	BackTrackPartitionVert(wgtInfo,pg,partition);
	double t5 = GetTime();
	st_backtrack += t5-t4;
//	fprintf(stderr,"score end inner %d == %d\n",edgecut,pg->GetEdgecut(partition));
	return score;
}
