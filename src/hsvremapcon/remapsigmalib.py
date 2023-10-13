#!/usr/bin/env python3
""" Tool for vertical mass conserved remapping between different hybrid pressure sigma coefficients. Copyright (C) JJD Hooghiem.

This program is free software: you can redistribute it and/or modify it under the
terms of the GNU General Public License as published by the Free Software Foundation,
version 3. This program is distributed in the hope that it will be useful, but
WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with this
program. If not, see <http://www.gnu.org/licenses/>.
"""

import netCDF4 as nc 
import numpy as np
import sys
from scipy.interpolate import interp1d
import numpy as np
from hsvremapcon.levels import * 

# global mean gravitational acceleration
grav = 9.80665
# the earth radius in meters
radius = 6.371e6  
deg2rad = np.pi / 180.

def globarea(im=360, jm=180):
    """ Function calculates the surface area according to TM5 definitions"""
    dxx = 360.0 / im * deg2rad 
    dyy = 180.0 / jm * deg2rad 
    lat = np.arange(-90 * deg2rad, 90 * deg2rad, dyy)
    dxy = dxx * (np.sin(lat + dyy) - np.sin(lat)) * radius ** 2
    area = np.resize(np.repeat(dxy, im, axis=0) , [jm, im])
    return area

def mass2pa(massfield):
    ''' Compute from full level mass arrays to half level partial pressure.
    Based on the simple pressure equation:
    p= F/A 
    with F the force perpendicular to the area A. For air: we have F=mg. so that
    p = g/A * m 
    We ignore the small decrease of g with altitude and its variation over the earth; g/A is a constant. 
    Therefor, it is not actually used, as we will devide back in the end at the same column. 
    '''
    a  =np.shape(massfield)
    out=np.zeros((a[0]+1,a[1],a[2]),order='F') 
    out[1:,:,:]=massfield
    out=np.cumsum(out,axis=0)
    return out 

def interpolate(surface_pressure,field_in,field_out,gridin,gridout):
    inshape=np.shape(field_in)
    oushape=np.shape(field_out)
    if not (inshape[1]==oushape[1])&(inshape[2]==oushape[2]):
        print(inshape)
        print(oushape)
        print('input fields lat lon dont match bye bye')
        return exit()
    for i in range(inshape[1]): # latitude
        for j in range(inshape[2]): # longitude 
            p_in=surface_pressure[i,j]*hybrid_sigma_b[gridin] + hybrid_sigma_a[gridin]
            p_out=surface_pressure[i,j]*hybrid_sigma_b[gridout] + hybrid_sigma_a[gridout]
            field_out[:,i,j]=interp1d(p_in,field_in[:,i,j],kind='quadratic')(p_out)
            massdif=np.sum(np.diff(field_out[:,i,j]))-np.sum(np.diff(field_in[:,i,j]))
            # print("interpolated mass difference %s for latlon (%s,%s)" %(massdif,i,j))
    return 

def pa2mass(psuedo_p):
    '''
    inverse operation of mass2p
    psuedo_partal pressure, we are actually ommiting the total column constant area gravity, so we are not truly working with partial pressure.
    '''
    massfield=np.diff(psuedo_p,axis=0)
    return massfield

def sigma2p(at,bt,sp):
    """
    computes atmospheric pressure based on sp and the hybrid sigma coefficients at and ab
    """
    return at[:,np.newaxis,np.newaxis]+bt[:,np.newaxis,np.newaxis]*sp[np.newaxis,:,:]

def remapsigma(datain,targetlevel):
    """
    core function that performs the remapping of the input data to output data based on the desired targetlevels
    """
    dataout=out_init(targetlevel,datain)
    for key,data in datain.items():
        if key!='sp':
            inlevels=detect_levels(data)
            print("interpolating %s from %s to %s" % (key,inlevels,targetlevel) )
            interpolate(datain['sp'],mass2pa(data),dataout[key],inlevels,targetlevel)
            dataout[key]=pa2mass(dataout[key]) 
    return dataout
