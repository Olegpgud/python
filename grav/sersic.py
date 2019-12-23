import matplotlib.pyplot as plt
import numpy as np
from astropy.io import fits
from astropy.table import Table
from astropy.wcs import WCS
from astropy import units as u
from astropy.coordinates import SkyCoord
from astropy.nddata import Cutout2D
import astropy.io.fits.hdu.hdulist
from scipy import optimize
from scipy.optimize import curve_fit
from lmfit import Model
from lmfit import Minimizer, Parameters, report_fit

tbl = Table.read('out_NGC7619.fits')
r=tbl['SMA']
I_DN=tbl['INTENS']
I_err = tbl['INT_ERR']


image_hst = fits.open('hst_14219_34_wfc3_ir_f110w_drz.fits',memmap=True)
hst_header = image_hst['PRIMARY'].header
arc_pix = hst_header['D001ISCL']
r_kpc = r*arc_pix*0.260
I = I_DN


def I_R(x,I0,Re,n,bkg):
    return bkg+I0*np.exp(-(2*n-1./3+4./(405*n)+46./(25525*n**2)+131./(1148175*n**3)-2194697./(30690717750*n**4))*(x/Re)**(1./n))

def fcn2min(params):
    Re = params['Re'].value
    n = params['n'].value
    I0 = params['I0'].value
    bkg = params['bkg'].value
    model = bkg + I0*np.exp(-(2*n-1./3+4./(405*n)+46./(25525*n**2)+131./(1148175*n**3)-2194697./(30690717750*n**4))*(r_kpc/Re)**(1./n))
    return (I - model) / I_err

params = Parameters()
params.add('Re', value=5)
params.add('n', value=4)
params.add('I0', value=1000)
params.add('bkg', value=10)

minner = Minimizer(fcn2min, params)
result = minner.minimize(method = 'nelder')

report_fit(result)

plt.errorbar(r_kpc, I, yerr=I_err, label='data')
plt.plot(r_kpc, I_R(r_kpc,result.params['I0'].value,result.params['Re'].value,result.params['n'].value,result.params['bkg'].value), 'r-')
plt.xlim(xmin=0,xmax=15)
#plt.xscale('log')
plt.yscale('log')
plt.show()







