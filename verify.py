import otherFunctions as of
from GOPack import GOPack
import math
import numpy as np
gop=GOPack()

gop.readObj("cow.obj") #input file
gop.riffle(20)
T=gop.loadTangency()
A,B=of.affineNormalizer(T)
vec3 = np.zeros((gop.nodeCount,3),dtype=float)
radii3 = np.zeros(gop.nodeCount,dtype=float)
with open("f://test.obj","w") as f: #output file
    for i in range(gop.centers.shape[0]):
        center=A*gop.centers[i]+B
        sph_corr,sr=of.e_to_s(center.real,center.imag,gop.radii[i]*A)
        vec = of.s_pt_to_vec(sph_corr[0],sph_corr[1])
        vec3[i,:]=vec
        radii3[i]=sr
        f.write("v {} {} {} \n".format(vec[0], vec[1], vec[2]))

for i in range(gop.centers.shape[0]):
    curr_vec3=vec3[i]
    for f in gop.flowers[i]:
        conn_vec3=vec3[f]
        dist=math.sqrt((curr_vec3[0]-conn_vec3[0])**2+(curr_vec3[1]-conn_vec3[1])**2+(curr_vec3[2]-conn_vec3[2])**2)
        sumdist=radii3[i]+radii3[f]
        print(dist,sumdist)
    break