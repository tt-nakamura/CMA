import numpy as np
import matplotlib.pyplot as plt
from CMA import PlotBoundary, WNS

mu = 5
Xmax,Ymax = 4.3, 6

plt.figure(figsize=(5, 5*Ymax/Xmax))

PlotBoundary(mu, Xmax, Ymax,
             R='b', L='r', S='g', RLPS=':')

XY = [(0.4, 0.3), # region 1
      (0.65, 0.5),# region 2
      (0.8, 0.8), # region 3
      (1.3, 0.8), # region 4
      (0.5, 3),   # region 6a
      (0.1, 4.7), # region 6b
      (1.3, 2.5), # region 7
      (2.5, 1.75),# region 8a
      (2.5, 3.5), # region 8b
      (0.8, 4.5), # region 9
      (0.8, 4.85),# region 10
      (2.5, 4.5), # region 11
      (0.5, 5.5), # region 12
      (2.5, 5.5)] # region 13
scale = [0.1, 0.04, 0.07, 0.05, 0.1, 0.07, 0.1,
         0.4, 0.2, 0.05, 0.08, 0.2, 0.2, 0.2]

for (X,Y), scale in zip(XY, scale):
    s1,s2 = WNS(X,Y,mu,scale,True)
    plt.plot(s1[:,0],s1[:,1],'k')
    plt.plot(s2[:,0],s2[:,1],'k')

plt.xlabel(r'$X$')
plt.ylabel(r'$Y$')
#plt.axis('equal')
plt.axis([0,Xmax,0,Ymax])
plt.tight_layout()
plt.savefig('fig2.eps')
plt.show()
