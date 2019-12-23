#!/usr/bin/python3
# -*- coding: UTF-8 -*-
import matplotlib.pyplot as plt
import numpy as np
from astropy.io import fits
from astropy.table import Table
from astropy.wcs import WCS
from astropy import units as u
from astropy.coordinates import SkyCoord
from astropy.nddata import Cutout2D
import astropy.io.fits.hdu.hdulist
from scipy.special import gamma
from scipy.special import gammaincc
from scipy.integrate import solve_ivp, RK23
from scipy.optimize import curve_fit
import os

def main(V0=None,alf=None):
    if V0 is None:
        V0 = float(input("input V0(m/sec): "))
    if alf is None:
        alf = float(input("input alf(degrees): "))


    '''
        Re:    5.44587833 (init= 5)
        n:     3.23984320 (init= 4)
        I0:    4967.29722 (init= 1000)
        bkg:   1.49366509 (init= 10)
    '''

    pi = 3.14
    G = 6.67*10**(-11)
    M0 = 2*10**30

    print('G = '+repr(G))

    path=os.path.dirname(os.path.abspath(__file__) )
    image_hst = fits.open(os.path.join(path,'hst_14219_34_wfc3_ir_f110w_drz.fits'),memmap=True)
    hst_header = image_hst['PRIMARY'].header
    arc_pix = hst_header['D001SCAL'] 
    zeropoint = hst_header['PHOTZPT']
    phot = hst_header['PHOTFLAM']
    phot2 = hst_header['PHOTPLAM']

    zeropoint = -2.5*np.log10(phot)-21.10-5*np.log10(phot2)+18.692
    print('zeropoint = '+repr(zeropoint))
    print('arc_pix = '+repr(arc_pix))
    m = -2.5*np.log10(4967.29722)+zeropoint
    m_c = m+2.5*np.log10(arc_pix**2)
    print('m = '+repr(m))
    print('m_c ='+repr(m_c))
    m_abs = m - 5*np.log10((53.6*10**6)/10)-0.1
    print('m_abs = '+repr(m_abs))
    m_abs_m = m_abs+2.5*np.log10(arc_pix**2*0.260**2*9.5*10**38)
    print('m_abs/m2 = '+repr(m_abs_m))
    I0_c = 10**((m_abs_m-3.82)/-2.5)*(3.08*10**16)**2
    print('I0_c = '+repr(I0_c))
    I0 = 3.8*10**26*10**((m_abs_m-3.82)/-2.5)
    print('I0 = '+repr(I0))

    Re = 5.44587833 * 3.08*10**19
    Re_pc = 5445.87833
    n = 3.23984320
    print('n = '+repr(n))


    b = 2*n-1./3+4./(405*n)+46./(25525*n**2)+131./(1148175*n**3)-2194697./(30690717750*n**4)
    p = 1.0 - 0.6097/n+0.05563/n**2
    M_L = 0.36
    print('M/L = '+repr(M_L))
    q0_mpc = M_L*I0_c*b**(n*(1-p))*gamma(2*n)/(2*Re_pc*gamma(n*(3-p)))
    q0 = q0_mpc*M0/(3.08*10**16)**3
    print('Re = '+repr(Re))
    print('b = '+repr(b)+' p = '+repr(p)+' q_mpc = '+repr(q0_mpc)+' q0 = '+repr(q0))


    def q(r):
        return q0*(r/Re)**(-p)*np.exp(-b*(r/Re)**(1/n))


    def L1(r):
        return q0*Re**2*n*b**(n*(p-2))*gammaincc(n*(2-p),b*(r/Re)**(1/n))*gamma(n*(2-p))


    f0 = -4*pi*G*L1(0)
    print('f0 = '+repr(f0))


    def L2(r):
        return q0*Re**3*n*b**(n*(p-3))*(gamma(n*(3-p))*gammaincc(n*(3-p),0)-gamma(n*(3-p))*gammaincc(n*(3-p),b*(r/Re)**(1/n)))

    def f(r):
        return -4*pi*G*(1/r*L2(r)+L1(r))

    def fun(t,r):
        return np.sqrt(abs(V0**2+2*f0-2*f(r)))

    def V_start(r,v_r):
        return np.sqrt(v_r**2-2*f0+2*f(r))



    #plot r(t) and v(t),v(r)
    t0=0
    tf=10**16
    max_step=tf/100

    x0=np.array([1])
    sol = solve_ivp(fun, [t0, tf], x0, max_step=max_step)
    print(sol)
    vel = fun(sol.t,sol.y[0])

    fig, axes = plt.subplots(figsize=(22,10), nrows=1, ncols=2)
    ax1 = axes[0]
    ax2 = axes[1]
    ax1.plot(sol.t/(3.154*10**13), (sol.y[0]/(3.08*10**19))*(abs(np.cos(alf*0.0175))))

    ax1.set_xlabel('time, millions yaers')
    ax1.set_ylabel('distance, kpc')
    ax1.set_title('R(t)')



    #ax2.plot(sol.y[0]/(3.08*10**19), (vel/1000)*(abs(np.sin(alf*0.0175))))
    ax2.plot(sol.t/(3.154*10**13), (vel/1000)*(abs(np.sin(alf*0.0175))))

    ax2.set_xlabel('time, millions yaers')
    ax2.set_ylabel('Velocity, km/sec')
    ax2.set_title('V(t)')

    plt.show()

    return vel












