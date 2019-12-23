#!/usr/bin/python3
# -*- coding: UTF-8 -*-
import numpy as np
from astropy.io import fits
from astropy.table import Table
from astropy.wcs import WCS
from astropy import units as u
from astropy import constants as const
from astropy.coordinates import SkyCoord
from astropy.nddata import Cutout2D
import astropy.io.fits.hdu.hdulist
from scipy.special import gamma
from scipy.special import gammaincc
from scipy.integrate import solve_ivp, RK23
from scipy.optimize import curve_fit
import os

pi=np.pi
G=const.G.value
M0=const.M_sun.value

n = 3.23984320 #sersic n
Re = 5.44587833 * 3.08e19 #effective radius
I0_sersic = 4967.29722 #I0

def main(V0=None,alf=None,flag=None):
    """
    Input:
        V0 - stars velocity
        alf - angular of trajectory
    Output:
        date of trajectory and velocity
        plot image 
    """
    if V0 is None:
        V0 = float(input("input V0(m/sec): "))
    if alf is None:
        alf = float(input("input alf(degrees): "))
    if flag is None:
        flag = float(input("1 - plot image, 0 - not: "))
    print(flag)

    path=os.path.dirname(os.path.abspath(__file__) )
    image_hst = fits.open(os.path.join(path,'hst_14219_34_wfc3_ir_f110w_drz.fits'),memmap=True)
    hst_header = image_hst['PRIMARY'].header
    arc_pix = hst_header['D001SCAL'] 
    zeropoint = hst_header['PHOTZPT']
    phot = hst_header['PHOTFLAM']
    phot2 = hst_header['PHOTPLAM']

    zeropoint = -2.5*np.log10(phot)-21.10-5*np.log10(phot2)+18.692 #hubble zero point
    m = -2.5*np.log10(I0_sersic)+zeropoint
    m_c = m+2.5*np.log10(arc_pix**2)
    m_abs = m - 5*np.log10((53.6e6)/10)-0.1
    m_abs_m = m_abs+2.5*np.log10(arc_pix**2*0.260**2*9.5e38)
    I0_c = 10**((m_abs_m-3.82)/-2.5)*(3.08e16)**2
    I0 = 3.8*10**26*10**((m_abs_m-3.82)/-2.5)
    Re_pc = Re/3.08e16 
    

    b = 2*n-1./3+4./(405*n)+46./(25525*n**2)+131./(1148175*n**3)-2194697./(30690717750*n**4) #polinome
    p = 1.0 - 0.6097/n+0.05563/n**2 #polinome
    M_L = 0.36
    q0_mpc = M_L*I0_c*b**(n*(1-p))*gamma(2*n)/(2*Re_pc*gamma(n*(3-p)))
    q0 = q0_mpc*M0/(3.08e16)**3

    def q(r):
        """
        rho function
        """
        return q0*(r/Re)**(-p)*np.exp(-b*(r/Re)**(1/n))

    def L1(r):
        return q0*Re**2*n*b**(n*(p-2))*gammaincc(n*(2-p),b*(r/Re)**(1/n))*gamma(n*(2-p))

    f0 = -4*pi*G*L1(0)

    def L2(r):
        return q0*Re**3*n*b**(n*(p-3))*(gamma(n*(3-p))*gammaincc(n*(3-p),0)-gamma(n*(3-p))*gammaincc(n*(3-p),b*(r/Re)**(1/n)))

    def f(r):
        """
        potencial fuction
        """
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
    vel = fun(sol.t,sol.y[0])
    if (flag==1):  
        import matplotlib.pyplot as plt     
        print(sol)
        fig, axes = plt.subplots(figsize=(22,10), nrows=1, ncols=2)
        ax1 = axes[0]
        ax2 = axes[1]
        ax1.plot(sol.t/(3.154e13), (sol.y[0]/(3.08e19))*(abs(np.cos(alf*0.0175))))
        ax1.set_xlabel('time, millions yaers')
        ax1.set_ylabel('distance, kpc')
        ax1.set_title('R(t)')
        ax2.plot(sol.t/(3.154e13), (vel/1000)*(abs(np.sin(alf*0.0175))))
        ax2.set_xlabel('time, millions yaers')
        ax2.set_ylabel('Velocity, km/sec')
        ax2.set_title('V(t)')
        plt.show()
    if (sol.success == False):
        return
    else:        
        return vel














