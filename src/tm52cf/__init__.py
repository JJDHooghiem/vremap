#!/usr/bin/env python
import sys
from hsvremapcon.load import load_TM5 
from hsvremapcon.write import write2netcdf
from hsvremapcon.remapsigmalib import *
from hsvremapcon.levels import *

def main():
    if sys.argv[1] in ['--help','-h']:
        print('Usage')
        print('     tm52cf index infile target_levels outfile')
        print('         index is an integer giving the position in the rm array with the desired tracer field to be used')
        print('         infile is the input netcdf file')
        print('         supported target levels are: %s' % list(hybrid_sigma_a.keys()))
        print('         outfile is the output netcdf file')
        exit()
    index = sys.argv[1]
    infile = sys.argv[2]
    targetlevel = sys.argv[3]
    ofile = sys.argv[4]

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
    date=infile.split('_')[2] 
    massratio=dataout['tracermassfield']/dataout['massfield']
    write2netcdf(ofile,massratio,targetlevel,date=date) 
    return
