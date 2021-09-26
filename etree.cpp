#include "etree.h"

int Etree::Part(int x){
	if(x < 0 || x >= NumNode()){
		return -1;
	}
	if(nodes[x].left != -1 || nodes[x].right != -1){
		return -1;
	}
	nodes[x].left = NumNode();
	nodes[x].right = NumNode() + 1;

	nodes.push_back(Enode());
	nodes[NumNode()-1].parent = x;
	nodes.push_back(Enode());
	nodes[NumNode()-1].parent = x;

	return NumNode();
}


std::vector<int> Etree::PreOrder(){
	std::vector<int> res;
	
	PreOrderSub(res, 0);

	return res;
}

int Etree::PreOrderSub(std::vector<int>& res, int pos){
	if(pos == -1) return 0;

	res.push_back(pos);
	
	PreOrderSub(res, nodes[pos].left);
	return PreOrderSub(res, nodes[pos].right);
}

