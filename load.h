#ifndef LOAD_H
#define LOAD_H

#include "etree.h"

typedef struct {
	int nvtxs;
	int* xadj;
	int* adjncy;
	int* vwgt;
	int* ewgt;
	int* cewgt;
	int* adjwgt;
} GraphData;

class SparseMatrix{
public:
	int m,n;
	int* rowptr; // size = m+1
	int* colind; // size = rowptr[m]
	double* val; // size = rowptr[m]

	SparseMatrix();
	SparseMatrix(int m, int n, int* rowptr, int* colind, double* val);
	SparseMatrix(const char* filename); // read MM format
	~SparseMatrix();

	void Reset();
	void Copy(SparseMatrix& src);
	void Copy(int m, int n, int* rowptr, int* colind, double* val);

	void Dump();
	void GenerateBitmap(const char* filename);

	void GenerateFile(const char* filename);

	void RemoveDiagonal();

	void SubMatInfo(Etree& etree);

	int* xadj(){
		return rowptr;
	}
	int* adjncy(){
		return colind;
	}
	int Vsize(){
		return m;
	}
	
	// each memory should be alloced
	void Get(GraphData* gd,bool vwgt_is_deg);

	// subroutine
	int tCSRLoadFromMM(const char* filename);
	int Expand();
};


#endif
