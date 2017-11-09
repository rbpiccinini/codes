import matplotlib.pyplot as plt
import numpy as np
from numpy import loadtxt
from scipy.sparse import coo_matrix
from scipy.sparse.linalg import eigsh

print('Reading matrix.')
#i,j,value=loadtxt('./27cellsAniso/linsys_mimetic-1-mat.dat').T
i,j,value=loadtxt('./lindat/linsys_mimetic-1-mat.dat').T
m=coo_matrix((value,(i,j)))*1e15
#m=m/m.max()
#m = m.todense()

print('Computing eigenvalues. It may take a while.')
w = eigsh(m, k=5, return_eigenvectors=False, which='LM', sigma=1e-20)
#w = eigsh(m, k=5, return_eigenvectors=False, which='SM', maxiter=1000)
print(w)

#plt.figure()
#for i in range(m.shape[0]):
#	for j in  range(m.shape[1]):
#		if m[i,j] > 0.:
#			plt.scatter(i,j)
#
#plt.savefig('matrix.png')
# plt.show()
