"""
Prepare table of gray Planck-mean absorption opacities as a function of the 
stellar color temperature. This is obtained by performing a Planck-weighted 
average of the frequency-dependent opacities of the Draine 2003 dust model, 
with the Planck-weight at the color temperature of the star. This code could
also work with any frequency-dependent opacity table. 
Paper Reference: Weingartner & Draine 2001
Reference for file: https://www.astro.princeton.edu/~draine/dust/dustmix.html - 

Instructions for other dust models :
By default, we choose the Milky Way, R_V = 5.5 grain model.                     
To choose another model from the above link provide input file as the file 
name in the ftp location of the table. E.g. for Milky Way, R_V = 3.1 model, 
open correponding link and provide its filename 
(i.e. kext_albedo_WD_MW_3.1_60_D03.all) as -infile. This will be automatically
downloaded from the link and used to create the table.

Python Requirements: numpy, scipy, astropy, argparse, and urllib

Author: Shyam Harimohan Menon (2021)
"""

import numpy as np
from scipy import integrate
from astropy.modeling import models
from astropy import units as u
import argparse
import timeit
import os

def downloadFile(url,infile):
    """
    Option to download frequency-dependent opacity table if file does not exist
    in provided directory.
    Arguments
        url : string
            URL to download file from
        infile: string
            File to store table in.

    """
    import urllib.request
    print("Retrieving table from {}".format(url))
    stat = urllib.request.urlretrieve(url,infile)
    print("Table retrieved and saved in {}".format(infile))
    return 



def planck_draine(Tstar,nu,kappa):
    """
    Returns planck-weighted gray-opacity for a given value of temperature
    See equation 15 of Menon+ 22 (VETTAM Method paper)
    Arguments:
        Tstar : float
            The temperature to calculate the Planck function at. 
        nu: ndarray
            The frequency values sampled in the input opacity table
        kappa: ndarray
            The opacity values corresponding to the frequencies sampled
    Returns:
        result: float
            Gray Planck opacity at given Tstar
    """
    bb = models.BlackBody(temperature=Tstar*u.K)
    num = bb.evaluate(nu,Tstar,scale=1)*kappa
    result = integrate.trapz(num,nu.value)/bb.bolometric_flux.value
    return result

def rosseland_draine(Tstar,nu,kappa):
    """
    Returns Rosseland-weighted gray-opacity for a given value of temperature
    See equation 21 of Menon+ 22 (VETTAM Method paper)
    Arguments:
        Tstar : float
            The temperature to calculate the Planck function at. 
        nu: ndarray
            The frequency values sampled in the input opacity table
        kappa: ndarray
            The opacity values corresponding to the frequencies sampled
    Returns:
        result: float
            Gray Rosseland opacity at given Tstar
    """
    
    bb = models.BlackBody(temperature=Tstar*u.K)
    h_planck, kb = 6.62606896e-27, 1.3807e-16
    exp_factor = np.exp(h_planck * nu.value/(kb*Tstar))
    exp_factor1 = np.exp((h_planck * nu.value)/(kb*Tstar) - 1)
    dB_dT = bb.evaluate(nu,Tstar,scale=1) * (exp_factor/exp_factor1) * h_planck * nu.value * (1./(kb * Tstar**2))
    num = integrate.trapz(dB_dT*(1./kappa),nu.value)
    denom = integrate.trapz(dB_dT,nu.value)
    result = num/denom
    return result


def prepareTable(infile='kext_albedo_WD_MW_5.5A_30_D03.all',
    noBins=1000,Tmin=3500.,Tmax=1.e5,log=False,
    outfile='DraineStellarOpacities.dat'):

    """
    Prepare table as discussed above. 
    Parameters:
        infile: string
            File where frequency-dependent opacities stored. 
        noBins: integer
            Number of color temperature bins to calculate gray opacity for
        Tmin: float
            Minimum color temperature in table
        Tmax: float 
            Maximum color temperature in table
        log : Boolean 
            Set true for log-spaced bins instead of linear bins
        outfile: string
            Output filename desired for table
    Return:
        None
    """
    file = os.path.abspath(infile)

    #Read wavelength and opacities from table
    lamda = np.loadtxt(file,usecols=0,skiprows=80) #Wavelength
    kappa = np.loadtxt(file,usecols=4,skiprows=80) #Opacity
    
    #Convert to frequencies
    nu = (lamda * u.micron).to(u.Hz, equivalencies=u.spectral())
    #Choose bins
    if(log):
        T = np.logspace(np.log10(Tmin),np.log10(Tmax),noBins)
    else:
        T = np.linspace(Tmin,Tmax,noBins)

    # Vectorise function calls
    vfunc = np.vectorize(planck_draine,excluded=['nu','kappa'])
    kappa_star = vfunc(Tstar=T,nu=nu,kappa=kappa)

    write_data = np.column_stack((T,kappa_star))
    outfile = os.path.abspath(outfile)
    np.savetxt(outfile,write_data,fmt='%5.3f')

    #Rosseland opacities
    vfunc = np.vectorize(rosseland_draine,excluded=['nu','kappa'])
    kappa_star = vfunc(Tstar=T,nu=nu,kappa=kappa)

    outfile = outfile + '_Rosseland'
    write_data = np.column_stack((T,kappa_star))
    outfile = os.path.abspath(outfile)
    np.savetxt(outfile,write_data,fmt='%5.3f')

    return

if __name__ == "__main__":

    # time the script
    start_time = timeit.default_timer()

    #Parsing Arguments
    ############################################################################

    ap = argparse.ArgumentParser(description=
        'Command Line Inputs for preparing gray-opacity table.')
    ap.add_argument('-infile',action='store',type=str,
        default='kext_albedo_WD_MW_5.5A_30_D03.all',
        help = 'Input file containing frequency-dependent opacities. Default \
./kext_albedo_WD_MW_5.5A_30_D03.all')
    ap.add_argument('-outfile',action='store',type=str,
        default='DraineStellarOpacities.dat',
        help = 'Output file to store prepared table. \
Default -  ./DraineStellarOpacities.dat')
    ap.add_argument('-noBins',action='store',type=int,default=1000,
        help='Number of color temperature bins required. Default 1000.')
    ap.add_argument('-Tmin',action='store',type=float,default=3500.,
        help='Minimum color temperature. Default is 3500K (Hayashi Limit)')
    ap.add_argument('-Tmax',action='store',type=float,default=1.e5,
        help='Maximum color temperature. Default is 10^5 K.') 
    ap.add_argument('-log',action='store_true',
        help='Set log-spaced color temperature bins. Default False')
    args = vars(ap.parse_args())
    ############################################################################

    infile = os.path.abspath(args['infile'])
    #Download input file if not present
    if(not os.path.exists(infile)):
        basename = os.path.basename(infile)
        #Assuming file would be present in Draine's ftp directory of dust models
        url = "ftp://ftp.astro.princeton.edu/draine/dust/dustmix/"+basename
        print("File not present. Attempting to download from Draine dust models\
 url: {}".format(url))
        downloadFile(url,infile)

    binType = 'Linearly'
    if(args['log']):
        binType = 'Log'

    print("Calculating gray-opacities from frequency-dependent opacities stored\
 in {}, for {} {}-spaced color temperatures in the range {} to {}.\
 Output table to be stored in {}.".format(args['infile'],args['noBins'],
    binType,args['Tmin'],args['Tmax'],args['outfile']))

    print("====================== START: Preparing Table ======================")
    prepareTable(args['infile'],args['noBins'],args['Tmin'],args['Tmax'],
        args['log'],args['outfile'])
    print("====================== END: Preparing Table ======================")

    # time the script
    stop_time = timeit.default_timer()
    total_time = stop_time - start_time
    print("***************** time to finish = "+
        str(total_time)+"s *****************")





    