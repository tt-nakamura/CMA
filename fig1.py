import matplotlib.pyplot as plt
from CMA import WNS

def XY(beta, gamma, mu=1836):
    alpha = mu*gamma*beta**2
    X = alpha*(1 + 1/mu)
    Y = mu*beta
    return X,Y

plt.figure(figsize=(5,7.5))

XY = [(1e6, 1e6), (2e6, 2e3), (11.5, 4.6),
      (1.1, 1.4), (0.9, 1.3), (0.7, 0.7)]

scale = [1.3, 0.06, 0.6, 1.6, 4, 2.1]
format = ['%.0f', '%.0f', '%g', '%g', '%g', '%g']

for i,(X,Y) in enumerate(XY):
    plt.subplot(3,2,i+1)
    u1,u2 = WNS(X,Y)
    s = scale[i]
    plt.axis('equal')
    plt.plot(u1[:,0],u1[:,1],'b')
    plt.plot(u2[:,0],u2[:,1],'r')
    plt.axis([-s,s,-s,s])
    if i%2==0: plt.ylabel(r'$u\cos\theta$')
    if i>=4:   plt.xlabel(r'$u\sin\theta$')
    t = plt.xticks()[0]; t = t[1]-t[0]
    f = r'$X=' + format[i] + '$, $Y=' + format[i] + '$'
    plt.tick_params(labelbottom=False, labelleft=False)
    plt.text(-s*0.9, s*0.9, f%(X,Y), va='center')
    plt.text(-s*0.9,-s*0.9, r'$L=%g$'%t, va='center')

plt.tight_layout()
plt.savefig('fig1.eps')
plt.show()
