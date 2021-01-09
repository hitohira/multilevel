#include "BAnetwork.h"
#include <time.h>

int main(int argc,char** argv){
	int N = 1000;
	if(argc == 2){
		N = atoi(argv[1]);
	}
	srand((unsigned)time(NULL));
	BAnetwork ba(N);
	ba.DumpDegrees();
//	ba.DumpGraph();
	return 0;
}
