#include "clustermap.h"
#include "part.h" 

void MatchData::Add(int mapidx, int vidx, PartGraph* pg){
	// update cewgt
	for(int i = 0; i < (int)match[mapidx].list.size(); i++){
		int uidx = match[mapidx].list[i];
		for(int j = pg->Xadj(uidx); j < pg->Xadj(uidx+1); j++){
			int eidx = pg->Adjncy(j);
			if(eidx == vidx){
				match[mapidx].cewgt += pg->Ewgt(j);
				break;
			}
		}
	}

	match[mapidx].list.push_back(vidx);
	match[mapidx].wgt += pg->Vwgt(vidx);

	if(mapidx >= used_n) used_n = mapidx+1;
}

