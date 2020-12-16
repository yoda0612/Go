import otherFunctions as of
from GOPack import GOPack
gop=GOPack()
gop.readObj("cow.obj") #input file
gop.riffle(1)
T=gop.loadTangency()
A,B=of.affineNormalizer(T)
with open("f://test.obj","w") as f: #output file

    for i in  range(gop.centers.shape[0]):
        center=A*gop.centers[i]+B
        sph_corr,_=of.e_to_s(center.real,center.imag,gop.radii[i]*A)
        vec = of.s_pt_to_vec(sph_corr[0],sph_corr[1])
        f.write("v {} {} {} \n".format(vec[0], vec[1], vec[2]))