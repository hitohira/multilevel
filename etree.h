#ifndef ETREE_H
#define ETREE_H

#include <cstdio>
#include <vector>

class Enode {
public:
	// set by this->Part()
	int parent;
	int left, right;

	// set when generating with proc info file
	int proc_id, unit_id;
	int weight; // sum of descendants

	// set when ordering
	int ofs, len;
	int order;

	Enode();
	Enode(const Enode& src){
		parent = src.parent;
		left = src.left;
		right = src.right;
		proc_id = src.proc_id;
		unit_id = src.unit_id;
		weight = src.weight;
		ofs = src.ofs;
		len = src.len;
		order = src.order;
	}
	Enode& operator=(const Enode& src){
		parent = src.parent;
		left = src.left;
		right = src.right;
		proc_id = src.proc_id;
		unit_id = src.unit_id;
		weight = src.weight;
		ofs = src.ofs;
		len = src.len;
		order = src.order;
		return (*this);
	}
};

typedef struct {
	int nodeid, procid, ncpu, ngpu;
} Pinfo;

class Etree {
private:
	std::vector<Enode> nodes;

	int PreOrderSub(std::vector<int>& res, int pos, int num);
	int SetOfsLenSub(std::vector<int>& blkptr, int step, int pos);
	int GenerateBlkInfoSub(FILE* fp, int pos);
	int setWeight(int wcpu, int wgpu, Pinfo* pinfo, int pos);
public:
	Etree() {
		nodes.push_back(Enode());
	}
	~Etree() {}
	
	int NumNode(){
		return nodes.size();
	}
	Enode Get(int pos){
		return nodes[pos];
	}

	// generate etree from proc info file
	int ConstructFromFile(const char* file_name);

	// etree have proc, unit, ofs, len info
	int SetOfsLen(std::vector<int>& blkptr);
	
	// divide leaf x, add leaf (NumNode) and (NumNode+1)
	int Part(int x);
	// output preorder node-num array, using for blk reordering
	std::vector<int> PreOrder();

	int GenerateBlkInfo(const char* file_name);

	void Dump(int depth, int pos);
};

#endif
