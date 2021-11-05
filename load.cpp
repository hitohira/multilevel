#include "load.h"
#include <stdio.h>
#include <stdlib.h>
#include <vector>
#include <algorithm>
#include <map>
#include <assert.h>

#include "bitmap.h"

SparseMatrix::SparseMatrix(){
	m = n = 0;
	rowptr = NULL;
	colind = NULL;
	val = NULL;
}
SparseMatrix::SparseMatrix(int m, int n, int* rowptr, int* colind, double* val){
	this->rowptr = NULL;
	this->colind = NULL;
	this->val = NULL;
	Copy(m,n,rowptr,colind,val);
}

void SparseMatrix::Copy(SparseMatrix& src){
	Copy(src.m,src.n,src.rowptr,src.colind,src.val);
}
void SparseMatrix::Copy(int m, int n, int* rowptr, int* colind, double* val){
	Reset();
	this->m = m;
	this->n = n;
	int nnz = rowptr[m];
	this->rowptr = (int*)malloc((m+1)*sizeof(int));
	this->colind = (int*)malloc(nnz*sizeof(int));
	this->val = (double*)malloc(nnz*sizeof(double));

	for(int i = 0; i <= m; i++){
		this->rowptr[i] = rowptr[i];
	}
	for(int i = 0; i < nnz; i++){
		this->colind[i] = colind[i];
		this->val[i] = val[i];
	}
}

// read MM format
SparseMatrix::SparseMatrix(const char* filename){
	tCSRLoadFromMM(filename);
	Expand();
}

SparseMatrix::~SparseMatrix(){
	Reset();
}

void SparseMatrix::Reset(){
	if(rowptr){
		free(rowptr);
		rowptr = NULL;
	}
	if(colind){
		free(colind);
		colind = NULL;
	}
	if(val){
		free(val);
		val = NULL;
	}
}

void SparseMatrix::Dump(){
	printf("m = %d, n = %d, nnz = %d\n",m,n,rowptr[m]);
}

void SparseMatrix::GenerateBitmap(const char* filename){
	if(m < 2000 || n < 2000){
		fprintf(stderr,"Generate Bmp : m = %d, nnz = %d, h = 1, w = 1\n",m,rowptr[m]);
		Bitmap bitmap(m,n);
		for(int i = 0; i < m; i++){
			for(int j = rowptr[i]; j < rowptr[i+1]; j++){
				unsigned x = i;
				unsigned y = colind[j];
				bitmap.set1(x,y);
			}
		}
		bitmap.write(filename);
	}
	else{
		int h = m / 1000;
		int w = n / 1000;
		fprintf(stderr,"Generate Bmp : m = %d, nnz = %d, h = %d, w = %d\n",m,rowptr[m],h,w);
		Bitmap bitmap(1000,1000);
		for(int i = 0; i < m; i++){
			for(int j = rowptr[i]; j < rowptr[i+1]; j++){
				unsigned x = i / h;
				unsigned y = colind[j] / w;
				bitmap.set1(x,y);
			}
		}
		bitmap.write(filename);
	}
}

void SparseMatrix::GenerateFile(const char* filename){
	FILE* fp = fopen(filename,"w");
	
	fprintf(fp,"%d %d\n", m, rowptr[m]);
	
	for(int i = 0; i < m; i++){
		fprintf(fp,"%d ", rowptr[i]);
	}
	fprintf(fp,"\n");

	for(int i = 0; i < rowptr[m]; i++){
		fprintf(fp,"%d ", colind[i]);
	}
	fprintf(fp,"\n");

	for(int i = 0; i < rowptr[m]; i++){
		fprintf(fp,"%lf ", val[i]);
	}
	fprintf(fp,"\n");

	fclose(fp);
}

void SparseMatrix::RemoveDiagonal(){
	int nnz = rowptr[m];

	int diagcntr = 0;
	for(int i = 0; i < m; i++){
		for(int j = rowptr[i]; j < rowptr[i+1]; j++){
			int v = colind[j];
			if(i == v){
				diagcntr++;
				break;
			}
		}
	}
	int* n_rowptr = (int*)malloc((m+1)*sizeof(int));
	int* n_colind = (int*)malloc((nnz-diagcntr)*sizeof(int));
	double* n_val = (double*)malloc((nnz-diagcntr)*sizeof(double));

	int counter = 0;
	for(int i = 0; i < m; i++){
		n_rowptr[i] = counter;
		for(int j = rowptr[i]; j < rowptr[i+1]; j++){
			int v = colind[j];
			if(i != v){
				n_colind[counter] = v;
				n_val[counter] = val[j];
				counter++;
			}
		}
	}
	n_rowptr[m] = counter;
	assert(counter == nnz - diagcntr);
	
	free(rowptr);
	free(colind);
	free(val);
	rowptr = n_rowptr;
	colind = n_colind;
	val = n_val;
}

void SparseMatrix::SubMatInfo(Etree& etree){
	std::vector<int> preOrder = etree.PreOrder();
	std::vector<int> ofs(preOrder.size());

	for(int i = 0; i < (int)preOrder.size(); i++){
		ofs[i] = etree.Get(preOrder[i]).ofs;
	}
	ofs.push_back(m);

	std::vector<std::vector<int> > res(preOrder.size());
	for(int i = 0; i < (int)preOrder.size(); i++) res[i] = std::vector<int>(preOrder.size(), 0);
	
	//count
	int r_bidx = 0;
	int r_next_ofs = ofs[r_bidx+1];
	for(int i = 0; i < m; i++){
		while(i >= r_next_ofs){
			r_bidx++;
			r_next_ofs = ofs[r_bidx+1];
		}

		for(int j = rowptr[i]; j < rowptr[i+1]; j++){
			int c = colind[j];
			std::vector<int>::iterator itr = std::lower_bound(ofs.begin(), ofs.end(), c);
			int c_bidx = std::distance(ofs.begin(), itr) - (*itr == c ? 0 : 1);
			res[r_bidx][c_bidx]++;
		}
	}

	//print
	printf("---submat info---\n");
	for(int i = 0; i < (int)preOrder.size(); i++){
		for(int j = 0; j < (int)preOrder.size(); j++){
			printf("(%8d, %8d, %8d)\t", ofs[i+1]-ofs[i], ofs[j+1]-ofs[j], res[i][j]);
		}
		printf("\n");
	}
	printf("---submat info---\n");

}


void SparseMatrix::Get(GraphData* gd, bool vwgt_is_deg){
	gd->nvtxs = m;

	for(int i = 0; i < m; i++){
		if(vwgt_is_deg){
			gd->vwgt[i] = rowptr[i+1] - rowptr[i] + 1;
		}
		else{
			gd->vwgt[i] = 1;
		}
//		gd->cewgt[i] = 0;
		gd->cewgt[i] = 1;
		gd->xadj[i] = rowptr[i];
	}
	gd->xadj[m] = rowptr[m];

	for(int i = 0; i < gd->xadj[m]; i++){
		gd->ewgt[i] = 1;
		gd->adjncy[i] = colind[i];
	}

	for(int i = 0; i < m; i++){
		int acc = 0;
		for(int j = gd->xadj[i]; j < gd->xadj[i+1]; j++){
			acc += gd->ewgt[j];
		}
		gd->adjwgt[i] = acc;
	}
}


int SparseMatrix::tCSRLoadFromMM(const char* filename){
	static char buf[4096];
	FILE* fp;
	fp = fopen(filename,"r");
	if(fp == NULL){
		fprintf(stderr,"fail at opening file\n");
		return -1;
	}
	// skip comment line
	while(1){
		if(fgets(buf,4096,fp) == NULL){
			fclose(fp);
			return -1;
		}
		if(buf[0] != '%'){
			break;
		}
	}
	int nz;
	sscanf(buf,"%d %d %d",&n,&m,&nz);
	val = (double*)malloc(nz*sizeof(double));
	colind = (int*)malloc(nz*sizeof(int));
	rowptr = (int*)malloc((m+1)*sizeof(int));
	rowptr[m] = nz;
	int tmpR = -1;
	for(int i = 0; i < nz; i++){
		if(i % 1000000 == 0) fprintf(stderr,"%d / %d\n",i,nz);
		int x,y;
		double f;
		if(fgets(buf,4096,fp) == NULL){
			fclose(fp);
			return -1;
		}
		int fn = sscanf(buf,"%d %d %lf",&y,&x,&f);
		if(!(fn == 2 || fn == 3)) {
			fclose(fp);
			m = n = 0;
			return -1;
		}
		if(fn == 3){
			val[i] = f;
		}
		else{
			val[i] = 1.0;
		}
		colind[i] = y-1;
		if(x-1 != tmpR){
			while(tmpR != x-1){
				tmpR++;
				if(tmpR <= 0){
					rowptr[tmpR] = 0;
				}
				else{
					rowptr[tmpR] = rowptr[tmpR-1];
				}
			}
			rowptr[x-1] = i;
			tmpR = x-1;
		}
	}
	for(int i = tmpR+1; i < m; i++){ // if row( > tmpR) has no nnz value, this loop is needed.
		rowptr[i] = rowptr[i-1];
	}
	fclose(fp);
	return 0;
}

int SparseMatrix::Expand(){
	std::vector<std::vector<std::pair<int,double> > > vlist(n);
	for(int i = 0; i < m; i++){
		for(int j = rowptr[i]; j < rowptr[i+1]; j++){
			vlist[colind[j]].push_back(std::make_pair(i,val[j]));
			if(i != colind[j]){
				vlist[i].push_back(std::make_pair(colind[j],val[j]));
			}
		}
	}
	int nnz = 0;
	for(int i = 0; i < m; i++){
		std::sort(vlist[i].begin(),vlist[i].end());
		nnz += vlist[i].size();
	}
	bool both = false;
	for(int i = 0; i < m; i++){
		for(int j = 1; j < (int)vlist[i].size(); j++){
			if(vlist[i][j-1] == vlist[i][j]){
				both = true;
				break;
			}
		}
		if(both) break;
	}
	if(both){ // upper and lower part exist
		fprintf(stderr,"both data exist\n");
		return 0;
	}
	// else upper or lower only
	int newN = m;
	int newM = n;
	free(rowptr);
	free(colind);
	free(val);
	rowptr = (int*)malloc((newM+1)*sizeof(int));
	colind = (int*)malloc(nnz*sizeof(int));
	val = (double*)malloc(nnz*sizeof(double));
	if(rowptr == NULL) {
		fprintf(stderr,"malloc error\n");
		return -1;
	}
	int cntr = 0;
	rowptr[newM] = nnz;
	for(int i = 0; i < newM; i++){
		rowptr[i] = cntr;
		for(int j = 0; j < (int)vlist[i].size(); j++){
			colind[cntr] = vlist[i][j].first;
			val[cntr] = vlist[i][j].second;
			cntr++;
		}
	}
	if(cntr != nnz) {
		fprintf(stderr,"not match nnz\n");
		return -1;
	}
	m = newM;
	n = newN;
	return 0;
}


