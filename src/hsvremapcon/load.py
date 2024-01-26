""" Tool for vertical mass conserved remapping between different hybrid pressure sigma coefficients. Copyright (C) JJD Hooghiem.

This program is free software: you can redistribute it and/or modify it under the
terms of the GNU General Public License as published by the Free Software Foundation,
version 3. This program is distributed in the hope that it will be useful, but
WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with this
program. If not, see <http://www.gnu.org/licenses/>.

Convention is to have the dataset to be interpolated with dimensions:
    levels, latitudes, longitudes.
    We are interpolating on levels, so these are put Fortran style memory layout for contigious access to memory. the top of the atmosphere is the zeroth index of the levels array, -1 is the index of the surface layer.  
"""
import netCDF4 as nc 
import numpy as np

def load_TM5(infile,i):
    '''
    Special function to read tm5 restart files. 
    infile, string, filename
    i, index in the tracer mass  array
    '''
    with nc.Dataset(infile) as ncf:
        sp = np.array(ncf['sp'][:,:],order='F')
        mass_in = np.array(ncf['m'][::-1,:,:],order='F')
        tracermassfield = np.array(ncf['rm'][:,::-1,:,:],order='F')
        names=nc.chartostring(ncf['names'][:])
    if 'co2_bg' in names:
        tracermassfield =np.sum(tracermassfield,axis=0)
    else:
        tracermassfield=tracermassfield[0]
    return {'massfield': mass_in, 'tracermassfield' : tracermassfield  , 'sp' : sp  } 
