#include "BAnetwork.h"
#include <time.h>
#include "part.h"

int main(int argc,char** argv){
	int N = 3000;
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

	for(int i = 0; i < N; i++){
		vwgt[i] = 1;
//		vwgt[i] = ba.xadj[i+1] - ba.xadj[i] + 1;
		cewgt[i] = 1;
		adjwgt[i] = ba.xadj[i+1] - ba.xadj[i];
	}
	for(int i = 0; i < nnz; i++){
		ewgt[i] = 1;
	}

	PartGraph pg(N,ba.xadj,ba.adjncy,vwgt,ewgt,cewgt,adjwgt);
	
	ndOptions options;
	SetDefaultOptions(&options);

	double ratioX = 0.5;
	if(1){	
		pg.Partition2(&options,ratioX,partition);	
		int mvc = pg.VertSepFromEdgeSep(partition);
		if(pg.VertSepIsOK(partition)){
			printf("%d\n",mvc);
		}
		else printf("vert sep is wrong\n");
	}
	if(1){
		pg.Partition3(&options,ratioX,partition);
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
	return 0;
}
