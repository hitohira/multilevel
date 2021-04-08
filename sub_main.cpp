#include "BAnetwork.h"
#include <time.h>
#include "part.h"
#include "ordering.h"

void PrintData(int n, int* partition,int* array){
	int sm[3] = {0,0,0};
	for(int i = 0; i < n; i++){
		sm[partition[i]] += array[i];
	}
	int total = sm[0] + sm[1] + sm[2];
	printf("%d %d %d %f %f\n",sm[0],sm[1],sm[2],1.0*sm[0]/total,1.0*sm[1]/total);
}

int main(int argc,char** argv){
	int N = 20000;
	if(argc == 2){
		N = atoi(argv[1]);
	}
	unsigned seed = (unsigned)time(NULL);
//	unsigned seed = 1610422192;
	fprintf(stderr,"%d\n",seed);
	srand(seed);
	BAnetwork ba(N);

//	ba.DumpDegrees();
//	ba.DumpGraph();
	int nnz = ba.xadj[N];

	int* vwgt = new int[N];
	int* ewgt = new int[nnz];
	int* cewgt = new int[N];
	int* adjwgt = new int[N];

	int* partition = new int[N];

	int* dbg = new int[N];

	for(int i = 0; i < N; i++){
//		vwgt[i] = 1;
		vwgt[i] = ba.xadj[i+1] - ba.xadj[i] + 1;
		cewgt[i] = 1;
		adjwgt[i] = ba.xadj[i+1] - ba.xadj[i];
		dbg[i] = ba.xadj[i+1] - ba.xadj[i] + 1;
	}
	for(int i = 0; i < nnz; i++){
		ewgt[i] = 1;
	}
	
	// CM
	/*
	{
		double* val = (double*)malloc(ba.xadj[N]*sizeof(double));
		for(int i = 0; i < ba.xadj[N]; i++) val[i] = 1.0;
		SparseMatrix mat(N, N, ba.xadj, ba.adjncy, val);
		SparseMatrix cmmat;
		IVec perm;
		CuthillMcKee(mat,perm);
		Rearrange(mat,perm,cmmat);
		mat.GenerateBitmap("mat.bmp");
		cmmat.GenerateBitmap("cm.bmp");
		free(val);
	}
	*/
	

	PartGraph pg(N,ba.xadj,ba.adjncy,vwgt,ewgt,cewgt,adjwgt);
	
	ndOptions options;
	SetDefaultOptions(&options);
	options.matchingScheme = MATCHING_CLUSTER;
//	options.ufactor = 50;
//	options.coarsenThreshold = 1000;
//	options.refineScheme = REFINE_VFM2;

	double ratioX = 0.5;
	if(0){	
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
