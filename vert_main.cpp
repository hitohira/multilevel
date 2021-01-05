#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include "part.h"
#include "load.h"

//const char* filename = "../centrality/matrix/bcsstk17/bcsstk17.mtx";
//const char* filename = "/mnt/d/DATA/Documents/IS/M1/Krylov/matrices/G3_circuit/G3_circuit.mtx";
//const char* filename = "/mnt/d/DATA/Documents/IS/M1/Krylov/matrices/posdef2/crystm02/crystm02.mtx";
//const char* filename = "/mnt/d/DATA/Documents/IS/M1/Krylov/matrices/ecology2/ecology2.mtx";
//const char* filename = "/mnt/d/DATA/Documents/IS/M1/Krylov/matrices/posdef/bundle1/bundle1.mtx";
const char* filename = "/mnt/d/DATA/Documents/IS/M1/Krylov/matrices/posdef2/parabolic_fem/parabolic_fem.mtx";
//const char* filename = "/mnt/d/DATA/Documents/IS/M1/Krylov/matrices/posdef2/nd12k/nd12k.mtx";
//const char* filename = "/mnt/d/DATA/Documents/IS/M1/Krylov/matrices/posdef2/crankseg_1/crankseg_1.mtx";
//const char* filename = "/mnt/d/DATA/Documents/IS/M1/Krylov/matrices/superlu/ldoor/ldoor.mtx";
//const char* filename = "/mnt/d/DATA/Documents/IS/M1/Krylov/matrices/superlu/Serena/Serena.mtx";

const double ratioX = 0.5;

int main(){
	unsigned seed = (unsigned)time(NULL);
	srand(seed);
	fprintf(stderr,"seed %u\n",seed);


	SparseMatrix csr(filename);	
	fprintf(stderr,"%s\n",filename);
//	csr.Dump();
	csr.RemoveDiagonal();
//	csr.Dump();
//	csr.GenerateBitmap("mat.bmp");
/*
	int nvtxs      = 5;
//	int nedges     = 12;
	int xadj[6]    = {0,     3,   5,     8,   10, 12};
	int adjncy[12] = {1,2,4, 0,2, 0,1,3, 2,4, 0,3};
	int vwgt[5]    = {1,1,1,1,1};
	int ewgt[12]   = {1,1,1, 1,1, 1,1,1, 1,1, 1,1};
	int cewgt[5]   = {0,0,0,0,0};
	int adjwgt[5]  = {3,2,3,2,2};
	int match[5],map[5];
*/
	int nvtxs = csr.Vsize();
	GraphData gd;
	gd.nvtxs = nvtxs;
	gd.xadj = (int*)malloc((nvtxs+1)*sizeof(int));
	gd.adjncy = (int*)malloc((csr.xadj()[nvtxs])*sizeof(int));
	gd.vwgt = (int*)malloc((nvtxs)*sizeof(int));
	gd.ewgt = (int*)malloc((csr.xadj()[nvtxs])*sizeof(int));
	gd.cewgt = (int*)malloc((nvtxs)*sizeof(int));
	gd.adjwgt = (int*)malloc((nvtxs)*sizeof(int));
	csr.Get(&gd,false);
//	csr.Get(&gd,true);

	int* match = (int*)malloc((nvtxs)*sizeof(int)); 
	int* map = (int*)malloc((nvtxs)*sizeof(int)); 
	
	int* partition = (int*)malloc((nvtxs)*sizeof(int));

	PartGraph pg(gd.nvtxs,gd.xadj,gd.adjncy,gd.vwgt,gd.ewgt,gd.cewgt,gd.adjwgt);

	ndOptions options;
	SetDefaultOptions(&options);
	
	for(int i = 0; i < 100; i++)
	{
		pg.Partition3(&options,ratioX,partition);
		if(pg.VertSepIsOK(partition))
//			printf("vert sep = %d\n",pg.GetVertSepSize(partition));
			printf("%d\n",pg.GetVertSepSizeCewgt(partition));
		else
			printf("vert sep wrong\n");
	}

	free(partition);

	free(gd.xadj);
	free(gd.adjncy);
	free(gd.vwgt);
	free(gd.ewgt);
	free(gd.cewgt);
	free(gd.adjwgt);
	free(match);
	free(map);
	return 0;
}
