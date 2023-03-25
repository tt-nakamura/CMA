# Clemmow-Mullaly-Allis diagram for waves in cold plasma
# references:
#   T. H. Stix "Waves in Plasmas" chapters 1-2
#   D. G. Swanson "Plasma Waves" 2nd edition, chapter 2

import numpy as np
from numpy.polynomial import Polynomial
import matplotlib.pyplot as plt

def PlotBoundary(mu, Xmax=None, Ymax=None, N=20, **kw):
    """
    mu = mass ratio of ion to electron
    Xmax = right limit of boundary
    Ymax = upper limit of boundary
    if Xmax is None, Xmax is set to mu-1+1/mu
    if Yamx is None, Ymax is set to mu*1.2
    N = number of plot points on boundary
    kw = keyward arguments
    kw['R'] = plot option of R=0 boundary
    kw['L'] = plot option of L=0 boundary
    kw['S'] = plot option of S=0 boundary
    kw['RLPS'] = plot option of RS=PS boundary
    """
    if Xmax is None: Xmax = mu - 1 + 1/mu
    if Ymax is None: Ymax = mu*1.2

    plt.hlines(1, 0, Xmax) # R\to\infty boundary
    plt.hlines(mu, 0, Xmax)# L\to\infty boundary
    plt.vlines(1, 0, Ymax) # P=0 boundary

    # R=0 boundary
    Y = np.linspace(0, 1, N)
    X = (1 + Y/mu)*(1 - Y)
    if 'R' in kw: plt.plot(X,Y,kw['R'])
    else: plt.plot(X,Y)

    # L=0 boundary
    Y = np.linspace(0, mu, N)
    X = (1 - Y/mu)*(1 + Y)
    if 'L' in kw: plt.plot(X,Y,kw['L'])
    else: plt.plot(X,Y)

    # S=0 boundary 1
    Y = np.linspace(0, 1, N)
    Y2 = Y**2
    X = (mu - Y2/mu)*(1 - Y2)/(mu - Y2)
    if 'S' in kw: plt.plot(X,Y,kw['S'])
    else: plt.plot(X,Y)

    # S=0 boundary 2
    X = np.linspace(0, Xmax, N)
    t = 1 - mu*X + mu**2
    Y2 = (t + np.sqrt(t**2 - 4*mu**2*(1-X)))/2
    Y = np.sqrt(Y2)
    if 'S' in kw: plt.plot(X,Y,kw['S'])
    else: plt.plot(X,Y)

    # RL = PS boundary
    Y = np.linspace(0, np.sqrt(mu*(mu-1) + 1), N)
    X = mu - 1 + (1-Y*Y)/mu
    if 'RLPS' in kw: plt.plot(X,Y,kw['RLPS'])
    else: plt.plot(X,Y)


def WNS(X, Y, mu=1836, scale=1, transl=False, N=200):
    """ Wave Normal Surface
    X = (plasma frequency)^2/(wave frequency)^2
    Y = (electron cyclotron freq)/(wave freq)
    mu = mass ratio of ion to electron
    scale = factor of magnification (or shrink)
    if transl, origin is set to (X,Y)
    N = number of plot points (maximum)
    return v1,v2 = two solutions for phase velocity
    each of v1,v2 has shape (n,2), n<=N
    v1[i,j] = j-th component of i-th point
    j=0,1 for x,y components, respectively
    """
    P = 1 - X
    R = 1 - X/(1 + Y/mu)/(1 - Y)
    L = 1 - X/(1 - Y/mu)/(1 + Y)
    S = (R+L)/2
    RL = R*L
    PS = P*S

    if PS<0:
        # resonance angle
        th_res = np.arctan(np.sqrt(-P/S))
        x = np.tanh(np.linspace(0,8,N//2))
        th1 = th_res * x
        th2 = np.pi/2*(1-x) + th1
        th = np.r_[th1,th2[::-1]]
    else:
        th = np.linspace(0, np.pi/2, N)
 
    c,s = np.cos(th), np.sin(th)
    c2,s2 = c*c, s*s

    A = S*s2 + P*c2
    B = RL*s2 + PS*(1 + c2)
    C = P*RL
    F = np.sqrt(B**2 - 4*A*C)

    u1 = (B+F)/2/C
    u2 = (B-F)/2/C
    w = []
    for u in [u1,u2]:
        r = scale * np.sqrt(u[u>=0])
        x,y = r*s[u>=0], r*c[u>=0]
        r = np.c_[np.r_[x, x[::-1], -x, -x[::-1]],
                  np.r_[y, -y[::-1], -y, y[::-1]]]
        if transl: r += [X,Y]
        w.append(r)
            
    return w


def GVS(X, Y, mu=1836, scale=1, transl=False, N=200):
    """ Group Velocity Surface
    reference: D.G.Swanson, section 2.3.2
    """
    Ym = Y/mu
    P = 1 - X
    R = 1 - X/(1 + Ym)/(1 - Y)
    L = 1 - X/(1 - Ym)/(1 + Y)
    S = (R+L)/2
    RL = R*L
    PS = P*S

    if PS<0:
        # resonance angle
        th_res = np.arctan(np.sqrt(-P/S))
        x = np.tanh(np.linspace(0,8,N//2))
        th1 = th_res * x
        th2 = np.pi/2*(1-x) + th1
        th = np.r_[th1,th2[::-1]]
    else:
        th = np.linspace(0, np.pi/2, N)
 
    c,s = np.cos(th), np.sin(th)
    c2,s2,sc2 = c*c, s*s, 2*s*c

    A = S*s2 + P*c2
    B = RL*s2 + PS*(1 + c2)
    C = P*RL
    F = np.sqrt(B**2 - 4*A*C)

    n1 = (B+F)/2/A
    n2 = (B-F)/2/A
    
    dA_dth = sc2*(S - P)
    dB_dth = sc2*(RL - PS)
    dF_dth = (B*dB_dth - 2*dA_dth*C)/F
    dP_dlnw = 2*X
    dR_dlnw = (1-R)*(1/(1+Ym) + 1/(1-Y))
    dL_dlnw = (1-L)*(1/(1-Ym) + 1/(1+Y))
    dS_dlnw = (dR_dlnw + dL_dlnw)/2
    dA_dlnw = dS_dlnw*s2 + dP_dlnw*c2
    dRL_dlnw = dR_dlnw*L + R*dL_dlnw
    dPS_dlnw = dP_dlnw*S + P*dS_dlnw
    dB_dlnw = dRL_dlnw*s2 + dPS_dlnw*(1+c2)
    dC_dlnw = dP_dlnw*RL + P*dRL_dlnw
    dlnA_dth = dA_dth/A
    dln1_dth = (dB_dth + dF_dth)/(B+F) - dlnA_dth
    dln2_dth = (dB_dth - dF_dth)/(B-F) - dlnA_dth
    alpha1 = -np.arctan(dln1_dth/2)
    alpha2 = -np.arctan(dln2_dth/2)

    vg = []
    for n, alpha in zip([n1,n2], [alpha1,alpha2]):
        i = (n>=0)
        n = n[i]
        v = (4*(2*A[i]*n - B[i])**2
             + (dA_dth[i]*n - dB_dth[i])**2)*n \
            /((dA_dlnw[i] - 4*A[i])*n**2
            - (dB_dlnw[i] - 2*B[i])*n + dC_dlnw)**2
        v = np.sqrt(v)
        x = v*np.sin(th[i] + alpha[i])
        y = v*np.cos(th[i] + alpha[i])
        v = np.c_[np.r_[x, x[::-1], -x, -x[::-1]],
                  np.r_[y, -y[::-1], -y, y[::-1]]]
        if transl: v += [X,Y]
        vg.append(v)

    return vg


def CPDR(k, theta, wp_wce=0.32, mu=1836):
    """ Cold Plasma Dispersion Relation
    k = (wave number)*c/(electron cyclotron freq.)
    theta = angle between k and B0 / radian
    wp_wce = (plasma freq)/(electron cyclotron freq.)
    mu = mass ratio of ion to electron
    return w = (wave freq.)/(electron cyclotron freq.)
    w has the same shape as k
    reference: D.G.Swanson, section 2.4
    """
    c2 = np.cos(theta)**2
    Wi = 1/mu # omega_ci / omega_ce
    Wi2 = Wi**2
    WW = 1 + Wi2
    ww = wp_wce**2
    wW = ww + Wi
    wW2 = wW**2
    k1,w = np.asarray(k),[]
    for k in np.atleast_1d(k1).flat:
        k2,k4 = k**2, k**4
        a8 = -(2*k2 + WW + 3*ww)
        a6 = k4 + (2*k2 + ww)*(WW + 2*ww) + wW2
        a4 = -(k4*(WW + ww) + 2*k2*wW2
             + k2*ww*(WW - Wi)*(1+c2) + ww*wW2)
        a2 = k4*(ww*(WW - Wi)*c2 + Wi*wW)\
           + k2*ww*Wi*wW*(1+c2)
        a0 = -k4*Wi2*ww*c2
        # equation(2.63) of Swanson
        a = Polynomial([a0,a2,a4,a6,a8,1])
        w.append(np.sqrt(a.roots()))

    return np.real(w).reshape(k1.shape + (-1,))
