# Vertical remapping tool for atmospheric trace gases

## Introduction 

This tool was developed to do vertical mass conserved remapping of CO~2~ fields from the so-called [TM5](https://tm5.sourceforge.net/) restart files to be used in [OpenIFS](https://confluence.ecmwf.int/display/OIFS/). Currently, the tool is limited to CO~2~ and TM5 restart files. The modular design aims to generalize this to other trace gases and or other input and output data types. Horizontal remapping is not provided since such operations on [netCDF](https://www.unidata.ucar.edu/software/netcdf/) can be achieved readily achieved by using tools like [CDO](https://code.mpimet.mpg.de/projects/cdo) or [NCO](https://nco.sourceforge.net/).

This tool focusses on exact mass conserved tracer remapping between different hybrid-sigma systems, while maintaining the concentration gradient. This is achieved by using quadratic spline interpolation using the partial pressure field of the tracer, defined at the half-level interfaces, as a vertical coordinate and the partial pressure of the trace gas of interest. Mass is exactly conserved in both upsampling (i.e. from a lower resolution vertical coordinate system to a higher) or when downsampling. The concentration gradient is approximated in the upsampling case.   

## Installation

Clone the repository and cd into the main directory. Build with
```sh
pip install . 
```
This will install a small library and a command line tool called tm52cf. 

This repo is as of today not present on the python package index. 

## Using the tool

Type 
```sh
tm52cf <--help|-h>
```
For a short help message.

```sh
tm52cf index infile target_levels outfile
```
Here `index` is the desired location of the tracer in tracer dimension of the TM5 restart file you wish to remap. The `infile` is the restart file itself. `target_levels` is the desired levels. For currently supported target levels see the help message from the tool above or the `src/hsvremapcon/levels.py`.

## Contributing

## Licensing
This code was published under the GNU General Public License V 3.0 of which you should have recieved a copy with the source code. 
