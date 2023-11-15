from scipy.optimize import newton
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

def VdW(p, t, r, a, b):
    """ This function calculates the molar volume from VdW's eos using Newton's method

    Parameters
    ----------
    p, t, r, a, b : int, optional
        All parameters of Van der Waal's equation of state in J/mol*K
    """
    def f(V):
        return r*t/(V-b)-a/(V**2)-p
    V_guess = r * t / p
    V = newton(f, V_guess)
    return V

def phimix(p,t,r,a,b,y1):
    """ This function calculates phi_mix from phi_1,mix and phi_2,mix

    Parameter
    ---------
    p, t, r : floats
        All parameters of Van der Waal's equation of state in J/mol*K
    a, b : double lists
        Each contain two values corresponding to a_i and b_i
    y1 : float
        mol fraction of species 1 in mixture
    """
    y2 = 1-y1
    y = [y1,y2]
    # Calculate a and b for VdW's eos
    am = a[0]*y[0]+a[1]*y[1]
    bm = b[0]*y[0]+b[1]*y[1]
    v = VdW(p,t,r,am,bm)

    # Calculate phi_i,mix using equations from table 7.1 in Koretsky
    phi1m = np.exp(b[0]/(v-bm)-np.log((v-bm)*p/(r*t))
                -2*(y[0]*a[0]+y[1]*(a[0]*a[1])**(1/2))/(r*t*v))
    
    phi2m = np.exp(b[1]/(v-bm)-np.log((v-bm)*p/(r*t))
                -2*(y[1]*a[1]+y[0]*(a[0]*a[1])**(1/2))/(r*t*v))
    
    # The sum of phi_i,mix*y_i should be phi_mix
    phim = y[0]*phi1m+y[1]*phi2m
    return phim

def phi_dy(p,t):
    """ This function calculates phi among various y1 values

    Parameters
    ----------
    p,t : int, optional
        The constant values of pressure and temperature

    Returns
    -------
    y : list
        An ordered list of y-values from 0 to 1
    phi : list
        Any index i of phi will return the fugacity coeficient of the mixture
        at y[i], t, and p.
    """
    # Initialize array of y values
    dy = 0.001
    y = np.arange(0,1+dy,dy)
    # Create array of phi values corresponding to y_1 values in y
    phi = [phimix(p,t,c[2],c[0],c[1],val)for val in y]
    return [y, phi]

def phi_dp(y,t):
    """ This function calculates phi among various P values

    Parameters
    ----------
    y,t : int, optional
        The constant values of mixture composition and temperature

    Returns
    -------
    p : list
        An ordered list of y-values from 1e4 to 1e7
    phi : list
        Any index i of phi will return the fugacity coeficient of the mixture
        at p[i], t, and y.
    """
    # Initialize array of P values
    dp = 1e4
    # It pisses me off how np.arange() isn't inclusive of the upper bound
    # We don't want to divide by 0 so we start at dp
    p = np.arange(dp,1e7+dp,dp)
    # Create array of phi values corresponding to P values in p
    phi = [phimix(val,t,c[2],c[0],c[1],y)for val in p]
    return [p,phi]

def phi_dt(y,p):
    """ This function calculates phi among various T values

    Parameters
    ----------
    y,p : int, optional
        The constant values of mixture composition and pressure

    Returns
    -------
    t : list
        An ordered list of y-values from 273 to 1073
    phi : list
        Any index i of phi will return the fugacity coeficient of the mixture
        at t[i], p, and y.
    """
    # Initialize array of T values of some arbitrary range (~200-1000?)
    dt = 0.8
    t = np.arange(273,1073+dt,dt)
    # Create array of phi values corresponding to T values in t
    phi = [phimix(p,val,c[2],c[0],c[1],y)for val in t]
    return [t,phi]

def graph(p,t,y1):
    """ This function 'compiles' the data into consumable portions

    Parameters
    ----------
    p,t,yi : floats
        Initialized values for constant properties
    """
    [y, phiy] = phi_dy(p,t)
    [dt,phit] = phi_dt(y1,p)
    [dp, phip] = phi_dp(y1,t)

    if txt:
        ydf = pd.DataFrame(y)
        phiydf = pd.DataFrame(phiy)
        yphi = pd.concat([ydf,phiydf], axis=1)
        yphi.columns = ['y1','phi']
        yphi.to_csv('data/dy.csv',index=False)

        tdf = pd.DataFrame(dt)
        phitdf = pd.DataFrame(phit)
        tphi = pd.concat([tdf,phitdf], axis=1)
        tphi.columns = ['t[K]','phi']
        tphi.to_csv('data/dt.csv',index=False)

        # Converts Pa to Bar
        pbar = (val/1e5 for val in dp)
        pdf = pd.DataFrame(pbar)
        phipdf = pd.DataFrame(phip)
        pphi = pd.concat([pdf,phipdf], axis=1)
        pphi.columns = ['p[bar]','phi']
        pphi.to_csv('data/dp.csv',index=False)

    # If we don't need any model, we can end here
    if not jpg and not figure:
        pass

    # Plotting y vs phi
    fig1 = plt.figure(1)
    ax1 = fig1.add_subplot(1,1,1)
    plt.plot(y,phiy)

    min1 = float(min(phiy))
    max1 = float(max(phiy))
    dphi1 = (max1-min1)/20
    mdphi1 = dphi1/5
    min1 = min1-dphi1
    max1 = max1+dphi1

    mt1 = np.arange(0,1.05,0.1)
    ax1.set_xticks(mt1)
    nt1 = np.arange(0,1.05,0.02)
    ax1.set_xticks(nt1,minor=True)
    myt1 = np.arange(min1,max1,dphi1)
    ax1.set_yticks(myt1)
    nyt1 = np.arange(min1,max1,mdphi1)
    ax1.set_yticks(nyt1,minor=True)

    plt.grid(which='major')
    plt.grid(which='minor',alpha=0.2)
    plt.title(rf'$y_\text{species1}$ vs $\phi$ ($T={t}[K]$ and $P={int(p*10**-5)}$[bar])')
    plt.xlabel(rf'$y_\text{species1}$')
    plt.ylabel(r'$\phi$')
    plt.xlim([0,1])
    plt.ylim([min1,max1])
    if jpg:
        plt.savefig('data/dy.jpg')

    # Plotting P vs phi
    fig1 = plt.figure(2)
    ax2 = fig1.add_subplot(1,1,1)
    plt.plot(dp*1e-5,phip)

    min2 = float(min(phip))
    max2 = float(max(phip))
    dphi2 = (max2-min2)/20
    mdphi2 = dphi2/5
    min2 = min2-dphi2
    max2 = max2+dphi2

    mt2 = np.arange(0,110,10)
    ax2.set_xticks(mt2)
    nt2 = np.arange(0,105,2)
    ax2.set_xticks(nt2,minor=True)
    myt2 = np.arange(min2,max2,dphi2)
    ax2.set_yticks(myt2)
    nyt2 = np.arange(min2,max2,mdphi2)
    ax2.set_yticks(nyt2,minor=True)

    plt.grid(which='major')
    plt.grid(which='minor',alpha=0.2)
    plt.title(rf'$P$[bar] vs $\phi$ ($T={t}[K]$ and $y_\text{species1}={y1}$)')
    plt.xlabel(r'$P$[bar]')
    plt.ylabel(r'$\phi$')

    plt.xlim([0,100])
    plt.ylim([min2,max2])
    if jpg:
        plt.savefig('data/dp.jpg')

    # Plotting T vs phi
    fig3 = plt.figure(3)
    ax3 = fig3.add_subplot(1,1,1)
    plt.plot(dt,phit)

    min3 = float(min(phit))
    max3 = float(max(phit))
    dphi3 = (max3-min3)/20
    mdphi3 = dphi3/5
    min3 = min3-dphi3
    max3 = max3+dphi3

    mt3 = np.arange(273,933,60)
    ax3.set_xticks(mt3)
    nt3 = np.arange(273,933,12)
    ax3.set_xticks(nt3,minor=True)
    myt3 = np.arange(min3,max3,dphi3)
    ax3.set_yticks(myt3)
    nyt3 = np.arange(min3,max3,mdphi3)
    ax3.set_yticks(nyt3,minor=True)

    plt.grid(which='major')
    plt.grid(which='minor',alpha=0.2)
    plt.title(rf'$T$[K] vs $\phi$ ($P={int(p*10**-5)}$[bar] and $y_\text{species1}={y1}$)')
    plt.xlabel(r'$T$[K]')
    plt.ylabel(r'$\phi$')

    plt.xlim([273,873])
    plt.ylim([min3,max3])
    if jpg:
        plt.savefig('data/dt.jpg')
    
    if figure:
        plt.show()
    pass

global c, txt, jpg, figure, c0, c1, t0, p0
# Initializing constant values - Methane and Propane by default
a=[0.1303, 0.939]
b=[4.31e-5 , 9.05e-5]
r=8.314
c = [a,b,r]
species1 = '{methane}'
omega = [0.008, 0.152]
tc = [190.6, 370]
pc = [46, 42.44]
jpg = False
txt = False
figure = False
