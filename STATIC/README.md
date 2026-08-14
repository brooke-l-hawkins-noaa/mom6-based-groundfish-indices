## Description

Static files necessary for running scripts to process groundfish indices

Description

Static files necessary for running scripts to process the longitudinal extent of groundfish indices based on 1) isobaths defined by bathymetry data and 2) distance from shore.

MOM6 distance from shore on native grid: shomom6-nep_dist2coast_CCS.nc

MOM6 distance from shore on the regrid: mom6-nep_regrid_dist2coast_CCS.nc

MOM6 bathymetry for isobaths: ocean_static.deptho.nc

To access/download bathymetry file for MOM6 directly from the thredds server:

``` r
library(reticulate)

# Import Python modules into R variables
xr   <- import("xarray")
np   <- import("numpy")
pd   <- import("pandas")
time <- import("time")

url <- "https://psl.noaa.gov/thredds/fileServer/Projects/CEFI/regional_mom6/cefi_derivative/northeast_pacific/full_domain/hindcast/monthly/regrid/r20260701/static/ocean_static.deptho.nc"
 dest_file <- "ocean_static.deptho.nc"

download.file(url, destfile = dest_file, method = "libcurl", mode = "wb")
```
