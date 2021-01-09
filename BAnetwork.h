#ifndef BANETWORK_H
#define BANETWORK_H

#include <vector>
#include <stdio.h>
#include <stdlib.h>

class BAnetwork {
private:
	int N;
	int ncores;
	int* xadj;
	int* adjncy;
public:
	BAnetwork(int N, int ncores=2){ // N = #vertices, start from $(ncores)-clique
		this->N = N;
		this->ncores = ncores;
		if(ncores < 2){
			fprintf(stderr,"ncores shold >= 2\n");
			exit(0);
		}

		xadj = new int[N+1];
		xadj[0] = 0;

		std::vector<std::vector<int> > edges(N);
		int nlink = 0;
		
		for(int i = 0; i < ncores; i++){
			for(int j = 0; j < ncores; j++){
				if(i != j){
					edges[i].push_back(j);
					nlink++;
				}
			}
		}
		nlink /= 2;

		for(int i = ncores; i < N; i++){
			for(int j = 0; j < i; j++){
				double p = 1.0 * edges[j].size() / nlink;
				if(drand() < p){
					edges[j].push_back(i);
					edges[i].push_back(j);
					nlink++;
				}
			}
		}
		adjncy = new int[nlink*2];
		int cntr = 0;
		for(int i = 0; i < N; i++){
			for(int j = 0; j < edges[i].size(); j++){
				adjncy[cntr++] = edges[i][j];
			}
			xadj[i+1] = xadj[i] + edges[i].size();
		}
		fprintf(stderr,"N = %d, nlink = %d\n",N,nlink);
	}
	~BAnetwork(){
		delete[] xadj;
		delete[] adjncy;
	}
	void DumpDegrees(){
		for(int i = 0; i < N; i++){
			printf("%d\n",xadj[i+1] - xadj[i]);
		}
	}
	void DumpGraph(){ // for show.py
		printf("%d %d\n",N,xadj[N]);
		for(int i = 0; i < N; i++){
			printf("%d %d 0\n",i,xadj[i+1] - xadj[i]);
		}
		for(int i = 0; i < N; i++){
			for(int j = xadj[i]; j < xadj[i+1]; j++){
				printf("%d %d 1\n",i,adjncy[j]);
			}
		}
	}

	double drand(){
		return 1.0 * rand() / RAND_MAX;
	}
};

#endif
