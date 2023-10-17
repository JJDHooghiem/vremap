#!/bin/sh

# first remap the restart file ('TM5_restart_20000108_0000_glb300x200.nc') vertically
# Note that the nc file is not provided by this repository
tm52cf 1 TM5_restart_20000108_0000_glb300x200.nc L91 temp.nc


# use an existing grib file containing the grid definition of the target
# this has the same layer definition as the above produce "temp.nc"
template=ICMGGECE4INIUA-co2
# note that the grib file 'ICMGGECE4INIUA-co2' is not included  
# here in this repo

# Use cdo to perform the horizontal remapping 
cdo -b 64 -P 8 --eccodes -f grb -remapbil,${template} temp.nc ICMGGECE4INIUA-co2_TM5

# clean up of files no longer needed
rm temp.nc
