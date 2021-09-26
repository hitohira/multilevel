#ifndef ETREE_H
#define ETREE_H

#include <cstdio>
#include <vector>

class Enode {
public:
	int parent;
	int left, right;
	Enode(){
		parent = left = right = -1;
	}
};

class Etree {
private:
	std::vector<Enode> nodes;

	int PreOrderSub(std::vector<int>& res, int pos);
public:
	Etree() {
		nodes.push_back(Enode());
	}
	~Etree() {}
	
	int NumNode(){
		return nodes.size();
	}

	// divide leaf x, add leaf (NumNode) and (NumNode+1)
	int Part(int x);
	// output preorder node-num array, using for blk reordering
	std::vector<int> PreOrder();
};

#endif
