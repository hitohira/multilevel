#ifndef LOAD_H
#define LOAD_H

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
private:
	int m,n;
	int* rowptr; // size = m+1
	int* colind; // size = rowptr[m]
	double* val; // size = rowptr[m]
public:
	SparseMatrix(const char* filename); // read MM format

	void Dump();
	void GenerateBitmap(const char* filename);

	void RemoveDiagonal();

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
	void Get(GraphData* gd);

	// subroutine
	int tCSRLoadFromMM(const char* filename);
	int Expand();
};


#endif
