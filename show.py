import sys
import networkx as nx
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import numpy as np
input = sys.stdin.readline

# --- imput data ---
# N M
# v0 vw0 p0
# v1 vw0 p1
# ...
# vN vwN pN
# a0 b0 ew0
# ...
# aM bM ewM
# --- end ---
# N : #vertices
# M : #edges
# vi vwi pi : i-vertex weight partition
# ai bi ewi : e(a1,b1) ewi=w(a1,b1)

def load_graph_data(fname):
	fp = open(fname)
	G = nx.Graph()
	color = []
	N,M = map(int,fp.readline().split())
	for i in range(N):
		v,vw,p = map(int,fp.readline().split())
		G.add_node(v,weight=vw)
		c = "red"
		if p == 1:
			c = "blue"
		color.append(c)
	for i in range(M):
		a,b,ew = map(int,fp.readline().split())
		G.add_edge(a,b,weight=ew)
	fp.close()
	return G,color

def draw_graph(G,color):
	nx.draw_networkx(G,node_color=color,node_size=[10]*G.number_of_nodes(),with_labels=False);
	plt.show()

G,color = load_graph_data("p.txt")
draw_graph(G,color)
