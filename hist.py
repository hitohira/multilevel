import sys
import matplotlib.pyplot as plt
import numpy as np
input = sys.stdin.readline

def load_data(fname):
	fp = open(fname,"r")
	X = map(int,fp.read().split())
	fp.close()
	return X 


X = load_data("p.txt")
X = np.array(X)
print("max",np.max(X))
print("ratio of 0",1.0*(len(X)-np.count_nonzero(X)) / len(X))
print("> 100", np.count_nonzero(X > 100))

plt.hist(X,bins=100)
plt.show()
