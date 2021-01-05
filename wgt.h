#ifndef FM_WGT_H
#define FM_WGT_H

typedef struct VGPair{
	int v, gain;
} VGPair;
	
typedef struct VGTri{
	int v, to, gain;
} VGTri;
typedef struct Trace{
	int v,from;
} Trace;
typedef struct WgtInfo{
	int wgt[3];
	double ratioX,tol; // ratioX - tol < wgt[0]/(wgt[0]+wgt[1]) < ratio + tol
} WgtInfo;

#endif
