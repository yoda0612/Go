import matplotlib.pyplot as plt
from GOPack import GOPack
gop=GOPack()
gop.readpack("rawTri")
gop.riffle(10)
x=gop.centers.real
y=gop.centers.imag

plt.figure()
plt.scatter(x,y)
plt.show()

