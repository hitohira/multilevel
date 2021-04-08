#ifndef ORDERING_H
#define ORDERING_H

#include <stdio.h>
#include <stdlib.h>
#include <vector>
#include <map>
#include <algorithm>
#include <omp.h>
#include <math.h>

#include "load.h"

class IVec{
public:
	int m;
	int* val;
	IVec(){
		m = 0;
		val = NULL;
	}
	IVec(int m){
		val = NULL;
		Alloc(m);
	}
	~IVec(){
		if(val) free(val);
	}
	void Alloc(int m){
		this->m = m;
		if(val) free(val);
		val = (int*)malloc(m*sizeof(int));
	}
	void Copy(IVec& vec){
		Alloc(vec.m);
		for(int i = 0; i < m; i++){
		val[i] = vec.val[i];
		}
	}
	void Dump(){
		printf("Dump begin %d\n",m);
		for(int i = 0; i < m; i++){
			printf("%d ",val[i]);
		}
		puts("");
		printf("Dump end\n");
	}
};

class CVec{
public:
	int m;
	double* c_val;
	CVec(){
		m = 0;
		c_val = NULL;
	}
	CVec(int m){
		c_val = NULL;
		Alloc(m);
	}
	CVec(int m,double v){
		c_val = NULL;
		Alloc(m);
		for(int i = 0; i < m; i++){
			c_val[i] = v;
		}
	}
	~CVec(){
		free(c_val);
	}
	void Alloc(int m){
		free(c_val);
		this->m = m;
		c_val = (double*)malloc(m*sizeof(double));
	}
	void SetZero(){
		for(int i = 0; i < m; i++){
			c_val[i] = 0.0;
		}
	}
	void SetRand(){
		for(int i = 0; i < m; i++){
			c_val[i] = rand() / 10000.0;
		}
	}
	void Dump(){
		printf("Dump begin %d\n",m);
		for(int i = 0; i < m; i++){
			printf("%f ",c_val[i]);
		}
		puts("");
		printf("Dump end\n");
	}
};


int CuthillMcKee(SparseMatrix& csr, IVec& perm);
int ReverseCuthillMcKee(SparseMatrix& csr, IVec& perm);
int Rearrange(SparseMatrix& oldcsr, IVec& perm, SparseMatrix& newcsr);
int SetPermFromPart(int nvtxs,int nparts, int* part, IVec& perm, IVec& partitionIdx);
int ReorderForCPUGPU(SparseMatrix& csr, IVec& partitionIdx, IVec& perm, IVec& sepIdx);
int PartitionSimple(SparseMatrix& csr, CVec& divRatio, IVec& partitionIdx);

#endif
