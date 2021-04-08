#include "ordering.h"

#include <algorithm>
#include <queue>
#include <map>
#include <set>
#include <functional>
#include <time.h>
#include <assert.h>


typedef std::pair<int, int> Pair;


Pair FindMinDegVertex(std::set<Pair> noVisit){
	return *(noVisit.begin());
}

int CuthillMcKee(SparseMatrix& csr, IVec& perm){
	// perf
	long long nfmin = 0;
	long long npop = 0;
//	double t1 = elasped();
	//perf
	int m = csr.m;
	perm.Alloc(m);
	
	for(int i = 0; i < m; i++){
		perm.val[i] = -1;
	}
	int inc = 0;
	std::set<Pair> noVisit;
	std::queue<Pair> que;
	std::vector<Pair> sortV;

	for(int i = 0; i < m; i++){
		int deg = csr.rowptr[i+1] - csr.rowptr[i];
		noVisit.insert(std::make_pair(deg,i));
	}

	while(!noVisit.empty()){
		Pair firstPair = FindMinDegVertex(noVisit);
		nfmin++; // perf nfmin
//		fprintf(stderr,"f %d %d\n",firstPair.first,firstPair.second);
		que.push(firstPair);
		perm.val[firstPair.second] = -2; // -2 : added, -1 not added, 0-(m-1): processed
		noVisit.erase(noVisit.begin());

		while(!que.empty()){
			Pair p = que.front(); que.pop();
			npop++; // perf
			int idx = p.second;
			perm.val[idx] = inc++;
			std::vector<Pair>().swap(sortV);
			for(int i = csr.rowptr[idx]; i < csr.rowptr[idx+1]; i++){
				int col = csr.colind[i];
				if(perm.val[col] == -1){ // if not added to que yet
					int deg = csr.rowptr[col+1] - csr.rowptr[col];
					sortV.push_back(std::make_pair(deg,col));
					perm.val[col] = -2;
					noVisit.erase(std::make_pair(deg,col));
				}
			}
			std::sort(sortV.begin(),sortV.end());
			for(int i = 0; i < (int)sortV.size(); i++){
				que.push(sortV[i]);
			}
		}
	}
	// perf
//	double t2 = elasped();
//	printf("%d %d %lld %lld %f\n",
//		csr.m,csr.rowptr[csr.m],nfmin,npop,(t2-t1)*1e6);
//	printf("nvert: %d, nedge: %d nfmin: %d, npop: %d, time: %f usec\n",
//		csr.m,csr.c_rowptr[csr.m],nfmin,npop,(t2-t1)*1e6);
	// perf
	return 0;
}

int ReverseCuthillMcKee(SparseMatrix& csr, IVec& perm){
	int res = CuthillMcKee(csr,perm);
	if(res != 0) return res;
	std::reverse(perm.val,perm.val+perm.m);
	return 0;
}


int Rearrange(SparseMatrix& oldcsr, IVec& perm, SparseMatrix& newcsr){
	int m = oldcsr.m;
	int n = oldcsr.n;
	int nnz = oldcsr.rowptr[m];

	double* n_val = (double*)malloc(nnz*sizeof(double));
	int* n_colind = (int*)malloc(nnz*sizeof(int));
	int* n_rowptr = (int*)malloc((m+1)*sizeof(int));
		
	n_rowptr[m] = nnz;

	IVec nzpr(m);
	for(int i = 0; i < m; i++){
		int row = perm.val[i];
		nzpr.val[row] = oldcsr.rowptr[i+1] - oldcsr.rowptr[i];
	}
	n_rowptr[0] = 0;
	for(int i = 1; i < m; i++){
		n_rowptr[i] = n_rowptr[i-1] + nzpr.val[i-1];
	}
	for(int i = 0; i < m; i++){
		int row = perm.val[i];
		int idx = n_rowptr[row];
		for(int j = oldcsr.rowptr[i]; j < oldcsr.rowptr[i+1]; j++){
			n_val[idx] = oldcsr.val[j];
			n_colind[idx] = perm.val[oldcsr.colind[j]];
			idx++;
		}
	}

	newcsr.Copy(m,n,n_rowptr,n_colind,n_val);
	free(n_val);
	free(n_colind);
	free(n_rowptr);
	return 0;
}

int SetPermFromPart(int nvtxs,int nparts, int* part, IVec& perm, IVec& partitionIdx){
	perm.Alloc(nvtxs);
	partitionIdx.Alloc(nparts+1);
	
	for(int i = 0; i < nparts; i++){
		partitionIdx.val[i] = 0;
	}
	for(int i = 0; i < nvtxs; i++){
		partitionIdx.val[part[i]+1]++;
	}
	for(int i = 1; i <= nparts; i++){
		partitionIdx.val[i] += partitionIdx.val[i-1];
	}
	IVec curIdx;
	curIdx.Copy(partitionIdx);
	for(int i = 0; i < nvtxs; i++){
		int idx = curIdx.val[part[i]]++;
		 assert(idx < nvtxs);
		 perm.val[i] = idx;
	}
	return 0;
}

int ReorderForCPUGPU(SparseMatrix& csr, IVec& partitionIdx, IVec& perm, IVec& sepIdx){
	fprintf(stderr,"ReorderForCPUGPU: this function is used for '1 CPU'-'1 GPU' hetero\n");

	int m = csr.m;
	perm.Alloc(m);
	
	// 0=CPU, 1=CPU2GPU, 2=GPU2CPU, 3=GPU
	int sep = partitionIdx.val[1];
	int count[4] = {0,0,0,0};
	count[0] = sep;
	for(int i = 0; i < sep; i++){
		for(int j = csr.rowptr[i]; j < csr.rowptr[i+1]; j++){
			if(csr.colind[j] >= sep){
				count[1]++;
				break;
			}
		}
	}
	count[0] -= count[1];
	count[3] = m - sep;
	for(int i = sep; i < m; i++){
		for(int j = csr.rowptr[i]; j < csr.rowptr[i+1]; j++){
			if(csr.colind[j] < sep){
				count[2]++;
				break;
			}
		}
	}
	count[3] -= count[2];

	int idx[4];
	idx[0] = 0;
	idx[1] = count[0];
	idx[2] = count[0] + count[1];
	idx[3] = count[0] + count[1] + count[2];
	sepIdx.Alloc(4);
	for(int i = 0; i < 4; i++) 
		sepIdx.val[i] = idx[i];
	
	for(int i = 0; i < sep; i++){
		bool stored = false;
		for(int j = csr.rowptr[i]; j < csr.rowptr[i+1]; j++){
			if(csr.colind[j] >= sep){
				perm.val[i] = idx[1]++;
				stored = true;
				break;
			}
		}
		if(!stored){
			perm.val[i] = idx[0]++;
		}
	}
	for(int i = sep; i < m; i++){
		bool stored = false;
		for(int j = csr.rowptr[i]; j < csr.rowptr[i+1]; j++){
			if(csr.colind[j] < sep){
				perm.val[i] = idx[2]++;
				stored = true;
				break;
			}
		}
		if(!stored){
			perm.val[i] = idx[3]++;
		}
	}
	return 0;
}

int PartitionSimple(SparseMatrix& csr, CVec& divRatio, IVec& partitionIdx){
	assert(divRatio.m == 2);
	int m = csr.m;
	int nparts = divRatio.m;
	int nnz = csr.rowptr[m];

	partitionIdx.Alloc(nparts+1);
	
	int sepnz = (int)(nnz * divRatio.c_val[0]);
	int sepIdx = (int)((long long)(std::lower_bound(csr.rowptr,csr.rowptr+m,sepnz)) - (long long)csr.rowptr)
	                 / sizeof(csr.rowptr[0]);
	fprintf(stderr,"partition row : %d\n",sepIdx);
	partitionIdx.val[0] = 0;
	partitionIdx.val[1] = sepIdx;
	partitionIdx.val[2] = m;
	return 0;
}

