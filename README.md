# Time_Series
This site contains the FORTRAN time-series analysis routines used in the publication: "Leaf phenology amplitude derived from MODIS NDVI and EVI: Maps of leaf phenology synchrony for Meso- and South America", which has been accepted for publication in Geoscience Data Journal. 

The FORTRAN time-series analysis routines were compiled for use with python (using f2py) and called from a general time series python control script. Here, we provide:
* the FORTRAN time-series analysis subroutines used by the python script (as single files)
   (a) Spectral_Analysis_routines_v1p1_noWRITES.f90, for processing large datasets
   (b) Spectral_Analysis_routines_v1p1.f90, for testing and debugging 
* FORTRAN control program to undertake a time series analysis (Spectral_Analysis_main_v1p1.f90)
* FORTRAN netCDF routines to input the time series data (Spectral_Analysis_READ_netCDF_v1p1.f90)

We also include:
* Zipped tar file with make files (Spectral_Analysis_make.tar.gz)
* Zipped tar file with example point and gridded time series, runtime input data and the corresponding output (Spectral_Analysis_Test_Data_Output.tar.gz)

The compilation of the FORTRAN source code requires netCDF libraries. The paths to the netCDF module and libraries in the Makefiles will need to be repointed. The Makefiles also use makedepf90 to generate the dependencies (see https://pkgs.org/download/makedepf90 or https://anaconda.org/conda-forge/makedepf90).

The point example corresponds to site 10 in the paper (Table 2 and Figure 8).

Garry Hayman, Centre for Ecology & Hydrology, garr@ceh.ac.uk
