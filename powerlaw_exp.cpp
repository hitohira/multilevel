#include "BAnetwork.h"
#include <time.h>
#include <stdio.h>
#include <stdlib.h>
#include "part.h"
#include "ordering.h"
#include "load.h"

int type_part_nnz = 0;
int type_sep_nnz = 0;
int type_ba = 0;
int type_powerlaw = 1;
int N = 20000;
double ratioX = 0.5;

void SetMyOptions(ndOptions* options){
	SetDefaultOptions(options);
	options->ufactor = 5;
	options->coarsenThreshold = 1000;

	if(type_ba || type_powerlaw){	
		options->matchingScheme = MATCHING_CLUSTER;
	}
	else{
		options->matchingScheme = MATCHING_RM;
	}

	if(type_part_nnz != type_sep_nnz){
		options->refineScheme = REFINE_VFM2;
	}
	else{
		options->refineScheme = REFINE_VFM;
	}
}

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

//const char* filename = "/mnt/d/DATA/Documents/IS/M1/Krylov/matrices/ScaleFree/com-Amazon/com-Amazon.mtx";
const char* filename = "/mnt/d/DATA/Documents/IS/M1/Krylov/matrices/ScaleFree/com-DBLP/com-DBLP.mtx";



void PrintData(int n, int* partition,int* array){
	int sm[3] = {0,0,0};
	for(int i = 0; i < n; i++){
		sm[partition[i]] += array[i];
	}
	int total = sm[0] + sm[1] + sm[2];
	printf("%d %d %d %f %f\n",sm[0],sm[1],sm[2],1.0*(sm[0]+sm[2])/total,1.0*(sm[1]+sm[2])/total);
}

int main(int argc,char** argv){
	PartGraph pg;
	int* vwgt = NULL;
	int* ewgt = NULL;
	int* cewgt = NULL;
	int* adjwgt = NULL;
	int* partition = NULL;
	int* dbg = NULL;

	if(type_ba){
		unsigned seed = (unsigned)time(NULL);
		fprintf(stderr,"%d\n",seed);
		srand(seed);
		BAnetwork ba(N);

		int nnz = ba.xadj[N];
	
		vwgt = new int[N];
		ewgt = new int[nnz];
		cewgt = new int[N];
		adjwgt = new int[N];

		partition = new int[N];
		dbg = new int[N];

		for(int i = 0; i < N; i++){
			if(type_part_nnz){
				vwgt[i] = ba.xadj[i+1] - ba.xadj[i] + 1;
			}
			else{
				vwgt[i] = 1;
			}
			cewgt[i] = 1;
			adjwgt[i] = ba.xadj[i+1] - ba.xadj[i];
			dbg[i] = ba.xadj[i+1] - ba.xadj[i] + 1;
		}
		for(int i = 0; i < nnz; i++){
			ewgt[i] = 1;
		}
		pg.SetGraph(N,ba.xadj,ba.adjncy,vwgt,ewgt,cewgt,adjwgt);
	}
	else{
		SparseMatrix csr(filename);
		fprintf(stderr,"%s\n",filename);
		csr.RemoveDiagonal();
		int N = csr.Vsize();
		int nnz = csr.xadj()[N];
		
		vwgt = new int[N];
		ewgt = new int[nnz];
		cewgt = new int[N];
		adjwgt = new int[N];

		partition = new int[N];
		dbg = new int[N];
		for(int i = 0; i < N; i++){
			dbg[i] = csr.xadj()[i+1] - csr.xadj()[i] + 1;
		}

		GraphData gd;
		gd.nvtxs = N;
		gd.xadj = new int[N+1];
		gd.adjncy = new int[nnz];
		gd.vwgt = vwgt;
		gd.ewgt = ewgt;
		gd.cewgt = cewgt;
		gd.adjwgt = adjwgt;

		csr.Get(&gd,type_part_nnz);
		pg.SetGraph(gd.nvtxs,gd.xadj,gd.adjncy,gd.vwgt,gd.ewgt,gd.cewgt,gd.adjwgt);

		delete[] gd.xadj;
		delete[] gd.adjncy;
	}

	ndOptions options;
	SetMyOptions(&options);

	if(0){
		// TODO impl clustering version for partition2	
		pg.Partition2(&options,ratioX,partition);	
		int mvc = pg.VertSepFromEdgeSep(partition);
		if(pg.VertSepIsOK(partition)){
			printf("%d\n",mvc);
		}
		else printf("vert sep is wrong\n");
	}
	if(1){
		pg.Partition3(&options,ratioX,partition);
		PrintData(pg.Vsize(),partition,dbg);
		PrintData(pg.Vsize(),partition,cewgt);
		if(pg.VertSepIsOK(partition))
			printf("%d\n",pg.GetVertSepSizeCewgt(partition));
		else
			printf("vert sep wrong\n");
	}

	delete[] vwgt;
	delete[] ewgt;
	delete[] cewgt;
	delete[] adjwgt;
	delete[] partition;
	delete[] dbg;
	return 0;
}


