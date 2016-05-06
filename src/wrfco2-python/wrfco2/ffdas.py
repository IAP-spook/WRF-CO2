#!/usr/bin/env python
"""
ffdas.py
FFDAS preprocessing Python module
C. Martin - 3/2016
"""
import numpy as np
import sys, time
from netCDF4 import Dataset

def prep(infile,outfile):
  """
  Function to preprocess FFDAS netCDF files
  Converts from 24 hourly variables to one
  flux variable varying in t dimension
  """
  
  # read in input netCDF file
  datain = Dataset(infile,'r')

  # open output netCDF file
  dataout = Dataset(outfile, "w", format="NETCDF4")

  # get dimensions from input file
  lat1 = datain.variables['latitude'][:]
  lon1 = datain.variables['longitude'][:]

  # get data from input file
  flux = np.zeros((24,len(lat1),len(lon1)))

  # read in each strange hourly variable and concatenate to one array
  for i in range(1,25):
    exec("tmpflux = datain.variables['flux_h%02d'][:]" % i)
    flux[i-1,:,:] = tmpflux

  # set up output file
  timeout = dataout.createDimension("time", None)
  lat2 = dataout.createDimension("latitude", len(lat1))
  lon2 = dataout.createDimension("longitude", len(lon1))
  times = dataout.createVariable("hour", "i2", ("time",))
  latitudes = dataout.createVariable("latitude","f4",("latitude",))
  longitudes = dataout.createVariable("longitude","f4",("longitude",))
  outflux = dataout.createVariable("flux","f8",("time","latitude","longitude" ,))
  timearr = np.arange(1,25,1)

  # add some metadata
  dataout.description = "Converted Hourly FFDAS flux netCDF file"
  dataout.history = "Created " + time.ctime(time.time())
  dataout.source = "convert_ffdas_hrly.py - C. Martin - Univ. of MD - 2/2016"
  latitudes.units = "degrees north"
  longitudes.units = "degrees east"
  outflux.units = "kgC/cell/h"
  times.units = "hour of day" 

  # write to file
  latitudes[:] = lat1
  longitudes[:] = lon1
  times[:] = timearr
  outflux[:] = flux

  # close files
  datain.close()
  dataout.close()



if __name__ == "__main__":
  if len(sys.argv) !=3:
    sys.exit('Usage: %s inputfile outputfile' % sys.argv[0])
  prep(sys.argv[1],sys.argv[2])
