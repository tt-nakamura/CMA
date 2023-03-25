import numpy as np
import matplotlib.pyplot as plt
from CMA import CPDR

wp_wce = 0.32
k = np.linspace(0,2,200)
th = [0, np.pi/6, np.pi/2]
label = [r'0$', r'\pi/6$', r'\pi/2$']

plt.figure(figsize=(5, 7.5))

for i,th in enumerate(th):
    plt.subplot(3,1,i+1)
    w = CPDR(k, th, wp_wce)
    plt.plot(k, w)
    plt.ylabel(r'$\omega/|\omega_{\rm ce}|$')
    if i<2: plt.tick_params(labelbottom=False)
    plt.text(0, 1.9, r'$\theta=' + label[i])

plt.xlabel(r'$ck/|\omega_{\rm ce}|$')
plt.tight_layout()
plt.savefig('fig4.eps')
plt.show()
