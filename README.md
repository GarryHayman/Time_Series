# Time_Series
This contains the FORTRAN time-series analysis routines used in the publication: "Leaf phenology amplitude derived from MODIS NDVI and EVI: Maps of leaf phenology synchrony for Meso- and South America", which has been accepted for publication in Geoscience Data Journal. 

The FORTRAN time-series analysis routines were compiled for use with python and called from a general time series python control script. Here, we provide:
* FORTRAN time-series analysis routines (as a single file)
   (a) Spectral_Analysis_routines_v1p1.f90, for processing large datasets
   (b) Spectral_Analysis_routines_v1p1_noWRITES.f90, for testing and debugging 

To follow shortly
-----------------
* FORTRAN control program to undertake a time series analysis
* FORTRAN netCDF routines to input the time series data
* Example time series and runtime input data are provided, together with the corresponding output.

Garry Hayman, Centre for Ecology & Hydrology, garr@ceh.ac.uk
