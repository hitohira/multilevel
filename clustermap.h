#ifndef CLUSTERMAP_H
#define CLUSTERMAP_H

#include <stdlib.h>
#include <vector>

class PartGraph;

class Cmatch {
public:
	std::vector<int> list;
	int wgt;
	int cewgt;
	Cmatch(){
		wgt = 0;
		cewgt = 0;	
	}
	~Cmatch(){}
};

class MatchData {
private:
	int n; // 余剰込みのmatch確保サイズ
	int used_n; // 実際に使われているmatchのサイズ
	Cmatch* match;

public:
	MatchData(int n){
		this->n = n;
		used_n = 0;
		match = (Cmatch*)malloc(n*sizeof(Cmatch));
	}
	~MatchData(){
		free(match);
	}
	int GetMappedSize(){
		return used_n;
	}
	int MatchSize(int i){
		return match[i].list.size();
	}
	int Match(int i,int j){
		return match[i].list[j];
	}
	int Weigt(int i){
		return match[i].wgt;
	}
	int Cewgt(int i){
		return match[i].cewgt;
	}
	
	// mapidxにvidxを追加
	void Add(int mapidx, int vidx, PartGraph* pg);
};

#endif
