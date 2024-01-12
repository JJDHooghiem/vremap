#!/usr/bin/env python
import sys
from hsvremapcon.load import load_TM5 
from hsvremapcon.write import write2netcdf
from hsvremapcon.remapsigmalib import *
from hsvremapcon.levels import *

# --------------------------------------------------------------------
# USE OF NOAA GML DATA
# 
# These data are made freely available to the public and the scientific
# community in the belief that their wide dissemination will lead to
# greater understanding and new scientific insights. To ensure that GML
# receives fair credit for their work please include relevant citation
# text in publications. We encourage users to contact the data providers,
# who can provide detailed information about the measurements and
# scientific insight.  In cases where the data are central to a
# publication, coauthorship for data providers may be appropriate.
# 
# 
# 
# Contact:  Xin Lan (xin.lan@noaa.gov)
# 
# File Creation: Fri Jan  5 03:55:24 2024
# 
# 
# --------------------------------------------------------------------
# 
# 
# See gml.noaa.gov/ccgg/trends/ for additional details.
# 
# The uncertainty in the global annual mean is estimated using a monte carlo
# technique that computes 200 global annual averages, each time using a
# slightly different set of measurement records from the NOAA GML cooperative
# air sampling network.  The reported uncertainty is the mean of the standard
# deviations for each annual average using this technique. Please see
# Conway et al., 1994, JGR, vol. 99, no. D11. for a complete discussion.
#
# CO2 expressed as a mole fraction in dry air, micromol/mol, abbreviated as ppm
#
# NOTE: In general, the data presented for the last year are subject to change, 
# depending on recalibration of the reference gas mixtures used, and other quality
# control procedures. Occasionally, earlier years may also be changed for the same
# reasons.  Usually these changes are minor.
# 
# year     mean      unc
import numpy as np 
from scipy.interpolate import interp1d
years=np.array([  1979,  1980,  1981,  1982,  1983,  1984,  1985,  1986,  1987,  1988,  1989,  1990,  1991,  1992,  1993,  1994,  1995,  1996,  1997,  1998,  1999,  2000,  2001,  2002,  2003,  2004,  2005,  2006,  2007,  2008,  2009,  2010,  2011,  2012,  2013,  2014,  2015,  2016,  2017,  2018,  2019,  2020,  2021,  2022, ]) 
co2= np.array([336.85,    338.91,    340.11,    340.86,    342.53,    344.07,    345.54,    346.97,    348.68,    351.16,    352.79,    354.06,    355.39,    356.09,    356.83,    358.33,    360.17,    361.93,    363.04,    365.70,    367.79,    368.96,    370.57,    372.58,    375.14,    376.95,    378.98,    381.15,    382.90,    385.02,    386.50,    388.76,    390.63,    392.65,    395.40,    397.34,    399.65,    403.07,    405.22,    407.61,    410.07,    412.44,    414.70,    417.07,    ])

f=interp1d(years,co2)


def main():
    if sys.argv[1] in ['--help','-h']:
        print('Usage')
        print('     tm52tm5 index infile target_levels year const')
        print('         index is an integer giving the position in the rm array with the desired tracer field to be used')
        print('         infile is the input netcdf file')
        print('         supported target levels are: %s' % list(hybrid_sigma_a.keys()))
        print('         outfile is the output netcdf file')
        exit()
    index = sys.argv[1]
    infile = sys.argv[2]
    targetlevel = sys.argv[3]
    thisyear =   sys.argv[4]
    try:
        taryear =   sys.argv[5]
    except:
        taryear=thisyear
    datain = load_TM5(infile,index)

    dataout=remapsigma(datain,targetlevel)
    # we can make this optional in the future
    checkmass=True
    if checkmass:
        print("----------------------  Summary ----------------------")
        dims=np.shape(datain['tracermassfield'])
        area_g=globarea(dims[2],dims[1])/grav
        for key in datain.keys():
            if key !='sp':
                print("input %s was %s kg" % (key,np.sum(datain[key]))) 
                print("output %s was %s kg" % (key,np.sum(dataout[key]))) 
        massout=pa2mass(  area_g[np.newaxis,:,:]*sigma2p( hybrid_sigma_a['L%s' % dims[0]], hybrid_sigma_b['L%s' % dims[0]] , datain['sp'] )) 
        print("input from surface pressure %s kg" % np.sum(massout) )
        massout=pa2mass(  area_g[np.newaxis,:,:]*sigma2p( hybrid_sigma_a[targetlevel], hybrid_sigma_b[targetlevel] , datain['sp'] )) 
        print("output from surface pressure %s kg" % np.sum(massout) )
    mass2mole=28.96/44.0
    fac=f(taryear)/f(thisyear)*mass2mole
    massratio=dataout['tracermassfield']/dataout['massfield']
    ofile = 'co2_init_%s.nc' % taryear
    write2netcdf(ofile,fac*massratio,targetlevel,date='%s0101' % taryear,reverse=True) 
    return
