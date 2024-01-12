""" Tool for vertical mass conserved remapping between different hybrid pressure sigma coefficients. Copyright (C) JJD Hooghiem.

This program is free software: you can redistribute it and/or modify it under the
terms of the GNU General Public License as published by the Free Software Foundation,
version 3. This program is distributed in the hope that it will be useful, but
WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with this
program. If not, see <http://www.gnu.org/licenses/>.
"""

import numpy as np
import netCDF4 as nc
from hsvremapcon.levels import hybrid_sigma_a, hybrid_sigma_b
import datetime as dt 
def calcmidbin(gridres,offset=-90,tot=180):
    nlatbins=tot/gridres
    binsize=tot/nlatbins
    return np.arange(nlatbins) * binsize +offset + binsize /2

def write2netcdf(ofile,dataset,targetgrid='L91',okey='co2',date='20000101',reverse=False):
    '''CF-compliant netcdf
      ofile is the target netcdf filename
      dataset is the data to be written
      targetgrid is the output vertical level system
      key name of the dataset
      date is the date of the dataset
    '''
    with nc.Dataset(ofile,'w') as ncf:
        # get the shape of input dataset
        ll=np.shape(dataset)

        # dimensions
        tdim=ncf.createDimension('time',None)
        latdim=ncf.createDimension('lat',ll[1])
        londim=ncf.createDimension('lon',ll[2])
        leveldim=ncf.createDimension('lev',ll[0])
        hleveldim=ncf.createDimension('nhyi',ll[0]+1)
        mleveldim=ncf.createDimension('nhym',ll[0])

        levelvar=ncf.createVariable('lev','f8',('lev',))
        levelvar[:] = 1+np.arange(ll[0])
        levelvar.setncattr('standard_name',"hybrid_sigma_pressure")
        levelvar.setncattr('long_name',"hybrid level at layer midpoints")
        levelvar.setncattr('formula',"hyam hybm (mlev=hyam+hybm*aps)")
        levelvar.setncattr('formula_terms',"ap: hyam b: hybm ps: aps")
        levelvar.setncattr('units',"level")
        if reverse:
            levelvar.setncattr('positive',"up")
        else:
            levelvar.setncattr('positive',"down")

        avar = ncf.createVariable('hyai','f8', ('nhyi',))
        avar.setncattr('long_name','hybrid A coefficient at layer interface')
        avar.setncattr('units','Pa')

        bvar = ncf.createVariable('hybi','f8', ('nhyi',))
        bvar.setncattr('long_name','hybrid B coefficient at layer interface')
        bvar.setncattr('units','1')
        if reverse:
            avar[:] = hybrid_sigma_a[targetgrid][::-1]
            bvar[:] =  hybrid_sigma_b[targetgrid][::-1]
        else: 
            avar[:] = hybrid_sigma_a[targetgrid] 
            bvar[:] =  hybrid_sigma_b[targetgrid] 

        date=dt.datetime.strptime(date,'%Y%m%d')


        time = ncf.createVariable('time','f8', ('time',))
        time[:] = np.array([(date-dt.datetime(2000,1,1)).days ])
        time.setncattr('long_name','time')
        time.setncattr('units','days since 2000-1-1 00:00:00'  )
        time.setncattr('axis','T')
        time.setncattr('calendar','standard')

        latvar=ncf.createVariable('lat','f8',('lat',))
        latvar[:] = calcmidbin(180/ll[1])
        latvar.setncattr('units','degrees_north')
        latvar.setncattr('axis','Y')
        latvar.setncattr('long_name','latitude')

        lonvar=ncf.createVariable('lon','f8',('lon',))
        lonvar[:] = calcmidbin(360/ll[2],offset=-180,tot=360)
        lonvar.setncattr('units','degrees_east')
        lonvar.setncattr('axis','X')
        lonvar.setncattr('long_name','longitude')

        tracervar = ncf.createVariable('co2','f8', ('time','lev','lat','lon',))
        if reverse:
            tracervar[:] = dataset[np.newaxis,::-1,:,:]
        else:
            tracervar[:] = dataset[np.newaxis,:,:,:]
        tracervar.setncattr('units','kg kg-1')
        tracervar.setncattr('long_name','co2 mass fraction in air')
        tracervar.setncattr('Parameter ID','210061')
    return
