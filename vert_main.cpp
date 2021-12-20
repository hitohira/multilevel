#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include "part.h"
#include "load.h"
#include "ordering.h"
#include <math.h>

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

//const char* filename = "/mnt/d/DATA/Documents/IS/M1/Krylov/matrices/sym/engine/engine.mtx";
//const char* filename = "/mnt/d/DATA/Documents/IS/M1/Krylov/matrices/sym/F1/F1.mtx";
//const char* filename = "/mnt/d/DATA/Documents/IS/M1/Krylov/matrices/sym/TEM152078/TEM152078.mtx";

//small
//const char* filename = "/mnt/d/DATA/Documents/IS/M1/Krylov/matrices/can_62/can_62.mtx";

//const char* filename = "/mnt/d/DATA/Documents/IS/M1/Krylov/matrices/ScaleFree/com-Youtube/com-Youtube.mtx";
//const char* filename = "/mnt/d/DATA/Documents/IS/M1/Krylov/matrices/ScaleFree/com-Amazon/com-Amazon.mtx";
//const char* filename = "/mnt/d/DATA/Documents/IS/M1/Krylov/matrices/ScaleFree/com-DBLP/com-DBLP.mtx";

//const char* filename = "/mnt/d/DATA/Documents/IS/M1/Krylov/matrices/ScaleFree/com-LiveJournal/com-LiveJournal.mtx";
//const char* filename = "/mnt/d/DATA/Documents/IS/M1/Krylov/matrices/ScaleFree/hollywood-2009/hollywood-2009.mtx";

const char* filename = "/mnt/d/DATA/Documents/IS/M1/Krylov/matrices/DIMACS10/coAuthorsCiteseer/coAuthorsCiteseer.mtx";
//const char* filename = "/mnt/d/DATA/Documents/IS/M1/Krylov/matrices/DIMACS10/m14b/m14b.mtx";


const double ratioX = 0.5;


void dump_deg(PartGraph*);

int main(){
//	unsigned seed = (unsigned)time(NULL);
//	seed = 1610591357;
//	srand(seed);
//	fprintf(stderr,"seed %u\n",seed);


	SparseMatrix csr_origin(filename);	
	SparseMatrix csr;
	fprintf(stderr,"%s\n",filename);
//	csr.Dump();
	csr.Copy(csr_origin);
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
	///!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
//	csr.Get(&gd,false);
	csr.Get(&gd,true);

	int* match = (int*)malloc((nvtxs)*sizeof(int)); 
	int* map = (int*)malloc((nvtxs)*sizeof(int)); 
	
	int* partition = (int*)malloc((nvtxs)*sizeof(int));
	int* local_perm = (int*)malloc(nvtxs*sizeof(int));

	PartGraph pg(gd.nvtxs,gd.xadj,gd.adjncy,gd.vwgt,gd.ewgt,gd.cewgt,gd.adjwgt);
	
	if(0){
		dump_deg(&pg);
		return 0;
	}

	ndOptions options;
	SetDefaultOptions(&options);
	options.matchingScheme = MATCHING_CLUSTER;
	options.ufactor = 5;
	options.coarsenThreshold = 1000;
	options.refineScheme = REFINE_VFM2;


	Etree etree;
	etree.ConstructFromFile("mynd_input6th.txt");
//	etree.ConstructFromFile("mynd_input5th.txt");
//	etree.ConstructFromFile("mynd_input4th.txt");
//	etree.ConstructFromFile("mynd_input3rd.txt");
//	etree.ConstructFromFile("mynd_input2nd.txt");
//	etree.ConstructFromFile("mynd_input.txt");
	etree.Dump(0,0);


	if(1){
		pg.NestedDissection(&options, etree, 0, partition, local_perm);
//		pg.Partition3(&options,ratioX,partition);
//		for(int i = 0; i < nvtxs; i++){
//			partition[i] = (partition[i]+1) % 3;
//		}
		printf("nd fin\n");
		
		SparseMatrix newcsr;
		IVec perm;
		GeneratePermFromEtreeCM(nvtxs, partition, local_perm, etree,perm);
//		GeneratePermFromEtree(nvtxs, partition,etree,perm);
		printf("perm fin\n");
		Rearrange(csr_origin,perm,newcsr);
		newcsr.GenerateBitmap("mynd.bmp");

		newcsr.GenerateFile("mynd_csr.txt");

		etree.GenerateBlkInfo("mynd_blk.txt");

		newcsr.SubMatInfo(etree);
		return 0;
	}

	
	for(int i = 0; i < 100; i++)
	{
		pg.Partition3(&options,ratioX,partition);
		if(pg.VertSepIsOK(partition))
//			printf("vert sep = %d\n",pg.GetVertSepSize(partition));
			printf("%d\n",pg.GetVertSepSizeCewgt(partition));
		else
			printf("vert sep wrong\n");
	}

	free(local_perm);
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

void dump_deg(PartGraph* pg){
	double acc = 0.0;
	double ave = 0.0;
	for(int i = 0; i < pg->Vsize(); i++){
		int deg = pg->Xadj(i+1) - pg->Xadj(i) + 1;
		ave += deg;
		acc += deg*deg;
		printf("%d\n",deg);
	}
	ave /= pg->Vsize();
	ave = ave*ave;
	acc /= pg->Vsize();
	fprintf(stderr,"%f\n",sqrt(acc - ave));
}
