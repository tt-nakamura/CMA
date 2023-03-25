import numpy as np
import matplotlib.pyplot as plt
from CMA import PlotBoundary, WNS

mu = 5

plt.figure(figsize=(5, 5.3))

Xmax,Ymax = 2, 6.15
plt.subplot(2,1,1)

PlotBoundary(mu, Xmax, Ymax,
             R='b', L='r', S='g', RLPS=':')

XY = [(0.2, 6),  # region 12
      (0.5, 5.1),# region 12
      (0.5, 5.6),# region 12
      (0.8, 6),  # region 12
      (1.2, 6),  # region 13
      (1.8, 5.1),# region 13
      (1.2, 5.2),# region 13
      (1.9, 5.9)]# region 13
scale = [0.1]*len(XY)

for (X,Y), scale in zip(XY, scale):
    s1,s2 = WNS(X,Y,mu,scale,True)
    plt.plot(s1[:,0],s1[:,1],'k')
    plt.plot(s2[:,0],s2[:,1],'k')

plt.ylabel(r'$Y$')
plt.tick_params(labelbottom=False)
plt.axis([0,Xmax,mu,Ymax])

#############################################

Xmax,Ymax = 2, 1
plt.subplot(2,1,2)

PlotBoundary(mu, Xmax, Ymax,
             R='b', L='r', S='g', RLPS=':')

XY = [(0.1, 0.8), # region 1
      (0.2, 0.2), # region 1
      (0.8, 0.1), # region 1
      (0.4, 0.75),# region 2
      (0.6, 0.6), # region 2
      (0.9, 0.3), # region 2
      (0.4, 0.9), # region 3
      (0.9, 0.5), # region 3
      (0.9, 0.9), # region 3
      (1.15,0.6)] # region 4
scale = [0.05]*len(XY)

for (X,Y), scale in zip(XY, scale):
    s1,s2 = WNS(X,Y,mu,scale,True)
    plt.plot(s1[:,0],s1[:,1],'k')
    plt.plot(s2[:,0],s2[:,1],'k')

plt.xlabel(r'$X$')
plt.ylabel(r'$Y$')
plt.axis([0,Xmax,0,Ymax])

plt.tight_layout()
plt.savefig('fig3.eps')
plt.show()
