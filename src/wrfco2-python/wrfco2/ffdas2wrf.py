#!/usr/bin/env python
"""
ffdas2wrf.py
Regrid FFDAS input .nc files
to a different resolution and domain
and output in a format for WRF flux input
C. Martin - 3/2016
"""
from netCDF4 import Dataset
import numpy as np
import sys
from scipy import interpolate
import time

def regrid2(infile,geofilename,datestr):
	"""
	Regrid FFDAS input .nc files to the specified
	WRF domain and outputs the proper file
	C. Martin - 7/2016
	- This version will get domain grid information from WRF geogrid files
	"""
	# read in input netCDF file
	datain = Dataset(infile,'r')
	# read in input geogrid file
	geofile = Dataset(geofilename,'r')
	# get domain info from Geogrid file
	latwrf = geofile.variables['XLAT_M'][0]
	lonwrf = geofile.variables['XLONG_M'][0]
	# get the four corners of the domain
	lat_s = latwrf[0,0];lat_n = latwrf[-1,-1]
	lon_w = lonwrf[0,0];lon_e = lonwrf[-1,-1]
	# get the resolution in degrees
	res = (lat_n - lat_s)/len(latwrf[0,:])
	# get input lat and lon
	lat1 = datain.variables['latitude'][:]
	lon1 = datain.variables['longitude'][:]
	# indices of data within WRF domain
	lat_inds = np.where(( lat1 >= lat_s-0.1) & (lat1 <= lat_n+0.1))
	lon_inds = np.where(( lon1 >= lon_w-0.1) & (lon1 <= lon_e+0.1))
	# subset the data
	flux = datain.variables['flux'][:]
	fluxsubset = flux[:,lat_inds[0],:]
	fluxsubset = fluxsubset[:,:,lon_inds[0]]
	lons, lats = np.meshgrid(lon1[lon_inds[0]], lat1[lat_inds[0]])
	# the FFDAS data is in units kg/cell/hr
	# must convert it to mol/km^2/hr for WRF
	mass_C = 12.
	area_div = 11**2 # assuming 11km for each 0.1 degree and squaring it for km^2
	# open the output file for writing
	# get the path to the output file
	# wrf looks for wrfchemi_dXX_YYYY-MM-DD_HH:MM:SS
	yearstr = datestr[0:4]; monstr = datestr[4:6]; daystr = datestr[6:8]
	geodir = geofilename.split('/')
	del geodir[-1]
	wrfdir = ''
	if geodir[0] == '':
		for a in range(len(geodir)-1):
	  		wrfdir = wrfdir+'/'+geodir[a+1]
	else:
		for a in range(len(geodir)):
			wrfdir = wrfdir+'/'+geodir[a]
	domain = geofilename[-5:]
	domain = domain[0:2]
	emissfile = 'wrfchemi_d'+domain+'_'+yearstr+'-'+monstr+'-'+daystr+'_00:00:00' # NOTE change this later to be more flexible with start time?
	# open output netCDF file
	dataout = Dataset(wrfdir+'/'+emissfile, "w", format="NETCDF4")
	# set up the regridded data
	newflux = np.zeros((24,1,len(latwrf[:,0]),len(latwrf[0])))
	# regrid the data
	for i in range(0,24):
		regridded = interpolate.griddata((lons.flatten(),lats.flatten()),fluxsubset[i].flatten(),(lonwrf,latwrf),method='linear')
		newflux[i,0,:,:] = (regridded*1000./mass_C)/area_div # convert to g and then moles then divide to get km^2
	# set up output file
	timeout = dataout.createDimension("Time", None)
	StrLength = dataout.createDimension("StrLength", 19)
	lat2 = dataout.createDimension("latitude", len(latwrf[:,0]))
	lon2 = dataout.createDimension("longitude", len(latwrf[0]))
	ez = dataout.createDimension("emissions_zdim",1)
	times = dataout.createVariable("Times", "S1", ("Time","StrLength"))
	co2 = dataout.createVariable("E_CO2","f4",("Time","emissions_zdim","latitude","longitude"))
	co2.setncattr("Sector","PMCH")
	co2.setncattr("FieldType", 104)
	timearr = []
	for i in range(24):
		exec("timestr1 = yearstr+'-'+monstr+'-'+daystr+'_%02d:00:00' % i")
		timearr.append(list(timestr1))
	# add some metadata
	dataout.description = "Regridded FFDAS flux netCDF file - "+str(res)+" degree"
	dataout.history = "Created " + time.ctime(time.time())
	dataout.source = "ffdas2wrf.py - C. Martin - Univ. of MD - 7/2016"
	co2.units = "mol/km^2/hr"
	dataout.setncattr("MMINLU", "USGS")
	dataout.setncattr("NUM_LAND_CAT", 24)
	# write to file
	times[:] = timearr
	co2[:] = newflux
	# close files
	datain.close()
	geofile.close()
	dataout.close()
		


def regrid(infile,outfile,lat_s,lat_n,lon_w,lon_e,res,datestr):
  # read in input netCDF file
  datain = Dataset(infile,'r')

  # open output netCDF file
  dataout = Dataset(outfile, "w", format="NETCDF4")

  # get dimensions from input file
  lat1 = datain.variables['latitude'][:]
  lon1 = datain.variables['longitude'][:]

  # indices
  lat_inds = np.where(( lat1 >= lat_s-0.1) & (lat1 <= lat_n+0.1))
  lon_inds = np.where(( lon1 >= lon_w-0.1) & (lon1 <= lon_e+0.1))
  
  # subset the data
  flux = datain.variables['flux'][:]
  fluxsubset = flux[:,lat_inds[0],:]
  fluxsubset = fluxsubset[:,:,lon_inds[0]]

  # regrid the data
  newlat = np.arange(lat_s+res/2.,lat_n-res/2.,res)
  #newlon = np.arange(lon_w+res/2.,lon_e-res/2.,res)
  newlon = np.arange(lon_w+res/2.,lon_e+res,res)

  # the FFDAS data is in units kg/cell/hr
  # must convert it to mol/km^2/hr for WRF
  mass_C = 12.
  area_div = (0.1/res)**2 # change in resolution squared to divide by (I think this is right?)

  lons, lats = np.meshgrid(lon1[lon_inds[0]], lat1[lat_inds[0]])
  newlons, newlats = np.meshgrid(newlon, newlat)
  newflux = np.zeros((24,1,len(newlat),len(newlon)))
  ppm400 = np.zeros((24,1,len(newlat),len(newlon)))

  for i in range(0,24):
    regridded = interpolate.griddata((lons.flatten(),lats.flatten()),fluxsubset[i].flatten(),(newlons,newlats),method='linear')
    newflux[i,0,:,:] = (regridded*1000./mass_C)/area_div # convert to g and then moles then divide by new area
    ppm400[i,0,:,:] = 400.


  # set up output file
  timeout = dataout.createDimension("Time", None)
  StrLength = dataout.createDimension("StrLength", 19)
  lat2 = dataout.createDimension("latitude", len(newlat))
  #lat2 = dataout.createDimension("south_north", len(newlat))
  #lon2 = dataout.createDimension("west_east", len(newlon))
  lon2 = dataout.createDimension("longitude", len(newlon))
  ez = dataout.createDimension("emissions_zdim",1)
  times = dataout.createVariable("Times", "S1", ("Time","StrLength"))
  #latitudes = dataout.createVariable("latitude","f4",("latitude",))
  #longitudes = dataout.createVariable("longitude","f4",("longitude",))
  co2 = dataout.createVariable("E_CO2","f4",("Time","emissions_zdim","latitude","longitude"))
  #co2_2 = dataout.createVariable("CO2_ANT","f4",("Time","emissions_zdim","latitude","longitude"))
  #co2 = dataout.createVariable("E_CO2","f4",("Time","emissions_zdim","south_north","west_east"))

  co2.setncattr("Sector","PMCH")
  co2.setncattr("FieldType", 104)
  #co2_2.setncattr("Sector","PMCH")
  #co2_2.setncattr("FieldType", 104)

  timearr = []

  yearstr = datestr[0:4]; monstr = datestr[4:6]; daystr = datestr[6:8]

  for i in range(24):
    exec("timestr1 = yearstr+'-'+monstr+'-'+daystr+'_%02d:00:00' % i")
    timearr.append(list(timestr1))

  # add some metadata
  dataout.description = "Regridded FFDAS flux netCDF file - "+str(res)+" degree"
  dataout.history = "Created " + time.ctime(time.time())
  dataout.source = "convert_ffdas_hrly.py - C. Martin - Univ. of MD - 2/2016"
  co2.units = "kgC/cell/h"
  #co2_2.units = "kgC/cell/h"

  dataout.setncattr("MMINLU", "USGS")
  dataout.setncattr("NUM_LAND_CAT", 24)

  # write to file
  times[:] = timearr
  co2[:] = newflux
  #co2_2[:] = ppm400

  # close files
  datain.close()
  dataout.close()
  
""" main script below """

if __name__ == "__main__":
  if len(sys.argv) != 9:
    sys.exit("Usage: %s input output lat_s lat_n lon_w lon_e res YYYYMMDD" % sys.argv[0])
  regrid(sys.argv[1],sys.argv[2],float(sys.argv[3]),float(sys.argv[4]),float(sys.argv[6]),float(sys.argv[7]),sys.argv[8])

