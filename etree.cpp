#include "etree.h"
#include <assert.h>

Enode::Enode(){
	parent = left = right = -1;
	proc_id = unit_id = -1;
	weight = 0;
	ofs = len = -1;
}

/* proc info file
	n_node n_proc n_unit
	cpu_power gpu_power (int)	
	(for each node)
	nodeid_i procid_i n_inner_cpu n_inner gpu
	...

 */


int Etree::ConstructFromFile(const char* file_name){
	int n_node, n_proc, n_unit;
	int wcpu, wgpu;
	Pinfo* pinfo;


	// read file
	FILE* fp = fopen(file_name,"r");
	if(fscanf(fp,"%d %d %d\n", &n_node, &n_proc, &n_unit) != 3){
		fclose(fp); return -1;
	}
	if(fscanf(fp,"%d %d\n",&wcpu, &wgpu) != 2){
		fclose(fp); return -1;
	}
	pinfo = new Pinfo[n_proc];
	for(int i = 0; i < n_proc; i++){
		if(fscanf(fp,"%d %d %d %d\n",&pinfo[i].nodeid, &pinfo[i].procid, &pinfo[i].ncpu, &pinfo[i].ngpu) != 4){
			fclose(fp); return -1;
		}
	}
	fclose(fp);

	// set etree
	int* pleaf = new int[n_proc]; // under pleaf[i], node of units in i-th proc is constructed
	int* uofs = new int[n_proc];
	// construct proc leaves. #leaves = nproc
	int nowN = 1;
	int pos = 0;
	while(nowN != n_proc){
		int res = Part(pos);
		if(res == -1) fprintf(stderr,"pos=%d\n",pos);
		assert(res != -1);
		pos++;
		nowN++;
	}
	for(int i = 0; i < n_proc; i++){
		pleaf[i] = NumNode() - 1 - i;
	}

	// construct unit leaves. devide proc leaves
	for(int i = 0; i < n_proc; i++){
		pos = pleaf[i];
		uofs[i] = i==0 ? 0 : uofs[i-1] + pinfo[i-1].ncpu + pinfo[i-1].ngpu;
		bool first = true;

		//set ancestor's id
		while(pos != -1 && nodes[pos].proc_id == -1){
			nodes[pos].proc_id = i;
			nodes[pos].unit_id = uofs[i];
			pos = nodes[pos].parent;
		}

		pos = pleaf[i];
		nodes[pos].weight = pinfo[i].ncpu > 0 ? -1 : -2; // -1=CPU, -2=GPU
		for(int j = 0; j < pinfo[i].ncpu - 1; j++){
			int res = Part(pos);
			nodes[res-2].proc_id = i;
			nodes[res-2].unit_id = nodes[pos].unit_id;
			nodes[res-2].weight = nodes[pos].weight;
			nodes[res-1].proc_id = i;
			nodes[res-1].unit_id = uofs[i] + j + 1;
			nodes[res-1].weight = -1; // cpu
			if(first){
				first = false;
				pos = res - 2;
			}
			else{
				pos++;
			}
		}
		int shift = pinfo[i].ncpu == 0 ? 1 : 0;
		for(int j = 0; j < pinfo[i].ngpu - shift; j++){
			int res = Part(pos);
			nodes[res-2].proc_id = i;
			nodes[res-2].unit_id = nodes[pos].unit_id;
			nodes[res-2].weight = nodes[pos].weight;
			nodes[res-1].proc_id = i;
			nodes[res-1].unit_id = uofs[i] + pinfo[i].ncpu + j + shift;
			nodes[res-1].weight = -2; //gpu
			if(first){
				first = false;
				pos = res - 2;
			}
			else{
				pos++;
			}
		}
	}

	// set weight info
	setWeight(wcpu, wgpu, pinfo, 0);
	
	delete[] uofs;
	delete[] pleaf;
	delete[] pinfo;
	return 0;
}
int Etree::setWeight(int wcpu, int wgpu, Pinfo* pinfo, int pos){
	if(nodes[pos].left == -1){ // leaf
		if(nodes[pos].weight == -2){ // GPU
			nodes[pos].weight = wgpu;
		}
		else{ // CPU
			nodes[pos].weight = wcpu;
		}
		return nodes[pos].weight;
	}
	// not leaf
	int wl = setWeight(wcpu, wgpu, pinfo, nodes[pos].left);
	int wr = setWeight(wcpu, wgpu, pinfo, nodes[pos].right);
	nodes[pos].weight = wl + wr;
	return wl + wr;
}


int Etree::Part(int x){
	if(x < 0 || x >= NumNode()){
		fprintf(stderr,"wrong x: %d / %d\n",x, NumNode());
		return -1;
	}
	if(nodes[x].left != -1 || nodes[x].right != -1){
		fprintf(stderr,"not leaf: l=%d r=%d\n",nodes[x].left, nodes[x].right);
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
	
	PreOrderSub(res, 0, 0);

	return res;
}

int Etree::PreOrderSub(std::vector<int>& res, int pos, int num){
	if(pos == -1) return num;

	res.push_back(pos);
	nodes[pos].order = num++;
	
	num = PreOrderSub(res, nodes[pos].left, num);
	return PreOrderSub(res, nodes[pos].right, num);
}

int Etree::SetOfsLen(std::vector<int>& blkptr){
	SetOfsLenSub(blkptr, 0, 0);
	return 0;
}

int Etree::SetOfsLenSub(std::vector<int>& blkptr, int step, int pos){
	if(pos == -1) return step;

	nodes[pos].ofs = blkptr[step];
	nodes[pos].len = blkptr[step+1] - blkptr[step];

	step = SetOfsLenSub(blkptr, step+1, nodes[pos].left);
	return SetOfsLenSub(blkptr, step, nodes[pos].right);
}

int Etree::GenerateBlkInfo(const char* file_name){
	FILE* fp = fopen(file_name, "w");
	if(fp == NULL){
		fprintf(stderr,"fail at gen blkinfo\n");
		return -1;
	}
	
	fprintf(fp,"%d\n",NumNode());

	GenerateBlkInfoSub(fp,0);

	fclose(fp);
	return 0;
}

int Etree::GenerateBlkInfoSub(FILE* fp, int pos){
	if(pos == -1) return 0;

	Enode en = nodes[pos];
	int is_sep = en.left != -1 ? 1 : 0;

	std::vector<int> tmp = PreOrder();
	int par = nodes[en.parent].order;
	fprintf(fp, "%d %d %d %d %d %d\n",en.proc_id, en.unit_id, par, en.ofs, en.len, is_sep);

	GenerateBlkInfoSub(fp, nodes[pos].left);
	return GenerateBlkInfoSub(fp, nodes[pos].right);
}

void Etree::Dump(int depth, int pos){
	if(pos == -1) return;

	for(int i = 0; i < depth; i++){
		fprintf(stderr,"-");	
	}
	fprintf(stderr,"%d: pid=%d uid=%d w=%d\n", pos, nodes[pos].proc_id, nodes[pos].unit_id, nodes[pos].weight);

	Dump(depth+1, nodes[pos].left);
	Dump(depth+1, nodes[pos].right);
}
