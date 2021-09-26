#ifndef BANETWORK_H
#define BANETWORK_H

#include <vector>
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>

class BAnetwork {
public:
	int N;
	int m0;
	int m;
	int* xadj;
	int* adjncy;


	BAnetwork(int N, int m0, int m){ // N = #vertices, m0 = #init_vertices, m = #edged added at each step
		this->N = N;
		this->m0 = m0;
		this->m = m;
		if(N < m0 || m0 < m){
			fprintf(stderr,"wrong args\n");
			exit(0);
		}

		xadj = new int[N+1];
		xadj[0] = 0;

		std::vector<std::vector<int> > edges(N);
		std::vector<double> prob(N+1);
		std::vector<bool> is_connected(N+1);
		int nlink = 0;

		fprintf(stderr,"N %d, m0 %d, m %d\n",N,m0,m);

		//set initial vertices. connected, tree
		for(int i = 1; i < m0; i++){
			double r = drand();
			int x = (int)(r*i);
			fprintf(stderr,"%d %d %f\n",x,i,r);
			assert(x < i);
			edges[x].push_back(i);
			edges[i].push_back(x);
		}
		nlink += m0-1;

		//add vertices
		for(int i = m0; i < N; i++){
			if(i % 10000 == 0) fprintf(stderr,"%d / %d\n",i,N);
			// init prob, is_connected
			prob[0] = 0.0;
			is_connected[0] = false;
			for(int j = 1; j <= i; j++){
				prob[j] = prob[j-1] + 1.0*edges[j-1].size() / nlink;
				is_connected[j] = false;
			}

			// add m edges at each step
			for(int j = 0; j < m; j++){
				while(1){
					double r = drand();
					int li = 0;
					int hi = i;
					int mi = (li+hi)/2;
					int ai = -1;
					// search
					while(1){
						assert(li < hi);
						if(prob[mi] > r){
							hi = mi;
							mi = (li + hi) / 2;
						}
						else{
							if(prob[mi+1] < r){
								li = mi+1;
								mi = (li + hi) / 2;
							}
							else{ // prob[mi] <= r <= prob[mi+1]
								ai = mi;
								break;
							}
						}
					}
//					fprintf(stderr,"r %d %f  %d\n",ai,r,i);
					if(is_connected[ai]) exit(0); // continue;
					edges[ai].push_back(i);
					is_connected[ai] = true;
					double pai = prob[ai+1] - prob[ai];
					for(int k = ai+1; k <= i; k++){
						prob[k] -= pai;
					}
					double maxp = prob[i];
					for(int k = 0; k <= i; k++){
						prob[k] /= maxp;
//						printf("%d %f\n",k,prob[k]);
					}
					break;
				}
			}
			for(int j = 0; j < i; j++){
				if(is_connected[j]) edges[i].push_back(j);
			}
			nlink += m;
		}

		adjncy = new int[nlink*2];
		int cntr = 0;
		for(int i = 0; i < N; i++){
			for(int j = 0; j < (int)edges[i].size(); j++){
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
		return 1.0 * (rand()-1) / RAND_MAX;
	}
};

#endif
