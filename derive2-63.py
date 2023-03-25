# derivation of equation(2.63) of 
#   D.G.Swanson "Plasma Waves" 2nd edition

import sympy as sp

w = sp.symbols('omega') # frequency
wp = sp.symbols('omega_p') # plasma frequency
Wi = sp.symbols('Omega_i') # ion cyclotron frequency
We = sp.symbols('Omega_e') # electron cyclotron frequency (We<0)
k = sp.symbols('k') # wave number
th = sp.symbols('theta') # angle between k and B0
c,s = sp.cos(th), sp.sin(th)
c2,s2 = c**2, s**2
n = k/w

P = 1 - (wp/w)**2
R = 1 - wp**2/(w+Wi)/(w+We)
L = 1 - wp**2/(w-Wi)/(w-We)
S = (R+L)/2
A = S*s2 + P*c2
B = R*L*s2 + P*S*(1 + c2)
C = P*R*L
D = A*n**4 - B*n**2 + C
E = D * w**6 * (w**2 - Wi**2) * (w**2 - We**2)
E = sp.expand(E)
E = sp.simplify(E)
E = sp.collect(E,w)

WW = -We*Wi
WW2 = We**2 + Wi**2
k2,k4 = k**2, k**4
wp2 = wp**2
E8 = -(2*k2 + WW2 + 3*wp2)
E6 = k4 + (2*k2 + wp2)*(WW2 + 2*wp2) + (wp2 + WW)**2
E4 = -(k4*(WW2 + wp2) + 2*k2*(wp2 + WW)**2
       + k2*wp2*(WW2 - WW)*(1+c2) + wp2*(wp2 + WW)**2)
E2 = k4*(wp2*(WW2-WW)*c2 + WW*(wp2 + WW)) + k2*wp2*WW*(wp2 + WW)*(1+c2)
E0 = -k4*WW**2*wp2*c2

for i,Ei in enumerate([E0,E2,E4,E6,E8]):
    a = E.coeff(w,2*i) - Ei
    print(sp.simplify(a)) # must be zero
