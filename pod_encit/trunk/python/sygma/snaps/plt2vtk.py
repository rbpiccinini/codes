#!/usr/bin/env python
import visit_writer
from pylab import *
import glob

files=glob.glob('sygma/pods/*.plt')

for file in files:
    x=loadtxt(file, skiprows=4)
    pts=[]
    vec=[]
    for i in range(len(x)):
        if sqrt(x[i,0]**2+x[i,1]**2) <= 64:
            pts=pts+[x[i,0],x[i,1],x[i,2]]
            vec=vec+[x[i,3],x[i,4],x[i,5]]

    vars = (("data", 3, 1, vec), ("ptsvec", 3, 1, pts))
    visit_writer.WritePointMesh(file[:-4]+'.vtk', 1, pts, vars)
