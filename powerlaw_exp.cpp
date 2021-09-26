#include "BAnetwork.h"
#include <time.h>
#include <stdio.h>
#include <stdlib.h>
#include "part.h"
#include "ordering.h"
#include "load.h"

int type_part_nnz = 1;
int type_sep_nnz = 1;
int type_ba = 1;
int type_powerlaw = 1;
int N = 20000;
double ratioX = 0.1;
int rep_times = 10;

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
//const char* filename = "/mnt/d/DATA/Documents/IS/M1/Krylov/matrices/ScaleFree/com-DBLP/com-DBLP.mtx";

const char* filename = "/mnt/d/DATA/Documents/IS/M1/Krylov/matrices/DIMACS10/coAuthorsCiteseer/coAuthorsCiteseer.mtx";
//const char* filename = "/mnt/d/DATA/Documents/IS/M1/Krylov/matrices/DIMACS10/m14b/m14b.mtx";

void PrintData(int n, int* partition,int* array){
	int sm[3] = {0,0,0};
	for(int i = 0; i < n; i++){
		sm[partition[i]] += array[i];
	}
	int total = sm[0] + sm[1] + 2*sm[2];
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

//	unsigned seed = (unsigned)time(NULL);
//	fprintf(stderr,"%d\n",seed);
//	srand(seed);

	if(type_ba){
//		BAnetwork ba(N);
//		BAnetwork ba(N,10,3); // data2
//		BAnetwork ba(N,10,10); // data3
		BAnetwork ba(N,10,6); // data4

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

	int* data_nnz = new int[3*rep_times];
	int* data_nvtxs = new int[3*rep_times];

	// partitioning
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
		for(int t = 0; t < rep_times; t++){
			pg.Partition3(&options,ratioX,partition);
//			PrintData(pg.Vsize(),partition,dbg);
//			PrintData(pg.Vsize(),partition,cewgt);
			int sm_nnz[3] = {0,0,0};
			int sm_nvtxs[3] = {0,0,0};
			for(int i = 0; i < pg.Vsize(); i++){
				sm_nnz[partition[i]] += dbg[i];
				sm_nvtxs[partition[i]] += cewgt[i];
			}
			data_nnz[t*3  ]   = sm_nnz[0];
			data_nnz[t*3+1]   = sm_nnz[1];
			data_nnz[t*3+2]   = sm_nnz[2];
			data_nvtxs[t*3  ] = sm_nvtxs[0];
			data_nvtxs[t*3+1] = sm_nvtxs[1];
			data_nvtxs[t*3+2] = sm_nvtxs[2];

			if(pg.VertSepIsOK(partition)){
				fprintf(stderr,"%d\n",pg.GetVertSepSizeCewgt(partition));
			}
			else{
				printf("vert sep wrong\n");
				break;
			}
		}
	}

	// print data
	printf("type_part_nnz,%d\n",type_part_nnz);
	printf("type_sep_nnz,%d\n",type_sep_nnz);
	printf("is_ba,%d\n",type_ba);
	printf("is_powerlaw,%d\n",type_powerlaw);
	printf("ratio,%f\n",ratioX);
	if(type_ba){
		printf("size,%d\n",N);
	}
	else{
		printf("name,%s\n",filename);
	}
	printf("nnz[0],nnz[1],nnz[2],nvtxs[0],nvtxs[1],nvtxs[2]\n");
	for(int t = 0; t < rep_times; t++){
		int k = t*3;
		printf("%d,%d,%d,%d,%d,%d\n",
			data_nnz[k],data_nnz[k+1],data_nnz[k+2],data_nvtxs[k],data_nvtxs[k+1],data_nvtxs[k+2]);
	}

	delete[] data_nvtxs;
	delete[] data_nnz;

	delete[] vwgt;
	delete[] ewgt;
	delete[] cewgt;
	delete[] adjwgt;
	delete[] partition;
	delete[] dbg;
	return 0;
}


