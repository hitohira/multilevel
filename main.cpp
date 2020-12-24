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
//const char* filename = "/mnt/d/DATA/Documents/IS/M1/Krylov/matrices/posdef2/parabolic_fem/parabolic_fem.mtx";
//const char* filename = "/mnt/d/DATA/Documents/IS/M1/Krylov/matrices/posdef2/nd12k/nd12k.mtx";
//const char* filename = "/mnt/d/DATA/Documents/IS/M1/Krylov/matrices/posdef2/crankseg_1/crankseg_1.mtx";
//const char* filename = "/mnt/d/DATA/Documents/IS/M1/Krylov/matrices/superlu/ldoor/ldoor.mtx";
//const char* filename = "/mnt/d/DATA/Documents/IS/M1/Krylov/matrices/superlu/Serena/Serena.mtx";

const double ratioX = 0.2;

int rec(PartGraph* pg,int* match,int* map){
	fprintf(stderr,"%d\t%d\n",pg->Vsize(), pg->Esize());
	if(pg->Vsize() < 200 || pg->Esize() < 200){
		int* partition = (int*)malloc(pg->Vsize()*sizeof(int));
		if(partition == NULL) fprintf(stderr,"mem X\n");
		ndOptions options;
		SetDefaultOptions(&options);
		pg->SetTolerance(&options);
		pg->InitPartitioning(ratioX,partition);
		/*
		FMDATA fm(pg,partition);
		int currWgtX = 0;
		int totalWgt = 0;
		for(int i = 0; i < pg->Vsize(); i++){
			if(partition[i] == 0) currWgtX += pg->Vwgt(i);
			totalWgt += pg->Vwgt(i);
		}
		int minWgtX = totalWgt*ratioX - pg->Tolerance();
		int maxWgtX = totalWgt*ratioX + pg->Tolerance();
		fprintf(stderr,"start refine %d %d %d\n",currWgtX,minWgtX,maxWgtX);
		currWgtX = fm.RefineEdge(currWgtX,minWgtX,maxWgtX,pg,partition);
		fprintf(stderr,"edge cut %d -> %d\n",old_ecut,fm.Edgecut());
		*/
		pg->Show(partition);
		free(partition);
		return 0;
	}
	int newSize = pg->Coarsening(match,map);
	PartGraph coarser;
	pg->GenerateCoarserGraph(newSize,match,map,&coarser);
	rec(&coarser,match,map);
	coarser.DeleteGraph();
	return 0;
}

int main(){
	srand((unsigned)time(NULL));
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
	csr.Get(&gd);

	int* match = (int*)malloc((nvtxs)*sizeof(int)); 
	int* map = (int*)malloc((nvtxs)*sizeof(int)); 
	
	int* partition = (int*)malloc((nvtxs)*sizeof(int));

	PartGraph pg(gd.nvtxs,gd.xadj,gd.adjncy,gd.vwgt,gd.ewgt,gd.cewgt,gd.adjwgt);

//	rec(&pg,match,map);
//	return 0;

	ndOptions options;
	SetDefaultOptions(&options);

	pg.Partition2(&options,ratioX,partition);

	fprintf(stderr,"edgecut %d\n",pg.GetEdgecut(partition));
	int wgtx = 0;
	for(int i = 0; i < nvtxs; i++){
		if(partition[i] == 0) wgtx += pg.Vwgt(i);
	}
	fprintf(stderr,"wgtX %d\n",wgtx);
//	pg.Show(partition);
	
	int boundary = 0;
	for(int i = 0; i < nvtxs; i++){
		for(int j = gd.xadj[i]; j < gd.xadj[i+1]; j++){
			if(partition[i] != partition[gd.adjncy[j]]){
				boundary++;
				break;
			}
		}
	}
	printf("boundary %d (%f)\n",boundary,1.0*boundary/nvtxs);
	
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
