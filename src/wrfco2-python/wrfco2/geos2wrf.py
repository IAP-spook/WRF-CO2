#!/usr/bin/env python
"""
geos2wrf.py
C. Martin - 3/2016
"""
import sys,time
'''
def readWRFnamelist(workdir):
  """ Reads WRF namelist using f90nml module
      to simplify things for the user
  """
  import f90nml
  nml = f90nml.read(workdir+'/namelist.input')
  start_year = nml['time_control']['start_year'][0]
'''
def geos2wrf(inputdir,workdir,domain):
  """ Reads in appropriate GEOSChem output,
  regrids and interpolates it.
  Modifies WRF boundary and initial conditions files
  to match GEOSChem input.
  """
  from netCDF4 import Dataset,date2index
  import numpy as np
  from scipy import interpolate
  from datetime import datetime,timedelta
  from cdo import Cdo 
  print("--------------------------------------------------")
  print("GEOS2WRF - C.Martin - 3/2016")
  print("NOTE: cdo must be on your path or errors will ensue!")
  print("Reading in WRF files for domain "+domain)
  print("from :"+workdir)
  wrfbdy = Dataset(workdir+'/wrfbdy_'+domain,'r+') # read wrf boundary file
  wrfin = Dataset(workdir+'/wrfinput_'+domain,'r+') # read wrf ic file

  # get basic information from these files
  wrft00 = ''.join(wrfin.variables['Times'][0])
  wrflats = wrfin.variables['XLAT'][0][:]
  wrflons = wrfin.variables['XLONG'][0][:]
  wrfeta = wrfin.variables['ZNU'][0][:]
  # four corners of domain
  wrflat_s = wrflats[0,len(wrflats)-1]
  wrflat_n = wrflats[len(wrflats)-1,len(wrflons[0])-1]
  wrflon_w = wrflons[0,0]
  wrflon_e = wrflons[len(wrflats)-1,len(wrflons[0])-1]
  print("Latitudinal extent: "+str(wrflat_s)+" to "+str(wrflat_n)+" degrees N")
  print("Longitudinal extent: "+str(wrflon_w)+" to "+str(wrflon_e)+" degrees E")
  # gridpoint size and number of gridpoints
  dx = (wrflon_e-wrflon_w)/len(wrflons[0])
  dy = (wrflat_n-wrflat_s)/len(wrflats)
  print("DX: "+str(dx)+" DY: "+str(dy))
  print("NX: "+str(len(wrflons[0]))+" NY: "+str(len(wrflats)))
  print(str(len(wrfeta))+" vertical levels")
  wrftimes = wrfbdy.variables['Times'][:]
  wrftimes = [''.join(wrftimes[t]) for t in range(len(wrftimes))]
  print("WRF simulation is for between: "+wrftimes[0]+" and "+wrftimes[len(wrftimes)-1])
  wrftimes = [datetime.strptime(wrftimes[t],"%Y-%m-%d_%H:%M:%S") for t in range(len(wrftimes))]

  # initial conditions first
  print(" * Initial Conditions *")
  geosfile = inputdir+'/'+'GEOSChem_'+datetime.strftime(wrftimes[0],"%Y%m%d")+'.nc'
  print("Interpolating GEOSChem output on to WRF grid...")
  # write the config file for CDO interpolation
  f = open(inputdir+'/WRFinterp.cdo','w')
  f.write('gridtype = lonlat\nxsize = '+str(len(wrflons[0]))+'\nysize = '+ str(len(wrflats))+'\nxfirst = '+str(wrflon_w)+'\nxinc = '+str((wrflon_e-wrflon_w)/len(wrflons[0]))+'\nyfirst ='+str(wrflat_s)+'\nyinc = '+str((wrflat_n-wrflat_s)/len(wrflats)))
  f.close()
  cdo = Cdo()
  cdo.remapbil(inputdir+'/WRFinterp.cdo',input=geosfile,output=geosfile+'.interp')
  geos_in = Dataset(geosfile+'.interp')
  print("Using initial conditions from: "+geosfile)
#  geos_t = date2index(wrftimes[0],geos_in.variables['time']) # get the index of the matching timestep
  geos_t = 0
  #geoslats = geos_in.variables['lat'][:]
  #geoslons = geos_in.variables['lon'][:]
  #geoslevs = geos_in.variables['lev'][:]
  geosco2 = geos_in.variables['co2'][:]
  wrfco2 = np.zeros((len(wrfeta),len(wrflats),len(wrflons[0])))
  z1,j1 = np.meshgrid(np.arange(0,len(wrflats)),np.linspace(1,0,47))
  b4 = (z1.flatten(),j1.flatten())
  z2,j2 = np.meshgrid(np.arange(0,len(wrflats)),wrfeta)
  after = (z2,j2)
  # interpolate through meridional slices
  for i in range(len(wrflons[0])):
    regridded = interpolate.griddata(b4,geosco2[geos_t,:,:,i].flatten(),after,method='linear')
    wrfco2[:,:,i] = regridded
  wrfin.variables['CO2_ANT'][0] = wrfco2 
  wrfin.close()
  geos_in.close()
  print("Initial conditions successfully modified!")

  # now boundary conditions
  if domain == 'd01':
    print(" * Domain is top level d01: * ") 
    print(" * Calculating Lateral Boundary Conditions * ") 
    # get boundary width 
    CO2XE = wrfbdy.variables['CO2_ANT_BXE'][:]
    CO2YE = wrfbdy.variables['CO2_ANT_BYE'][:]
    xdim = len(CO2XE[0,0,0,:]); ydim = len(CO2YE[0,0,0,:]); zdim = len(CO2XE[0,0,:,0]); bdywidth = len(CO2XE[0,:,0,0])
    starttime = wrftimes[0] ; endtime = wrftimes[len(wrftimes)-1]
    # empty arrays for placing output in
    CO2_XE = np.zeros((len(wrftimes),bdywidth,zdim,xdim))
    CO2_XS = np.zeros((len(wrftimes),bdywidth,zdim,xdim))
    CO2_YE = np.zeros((len(wrftimes),bdywidth,zdim,ydim))
    CO2_YS = np.zeros((len(wrftimes),bdywidth,zdim,ydim))
    # coordinates for interpolation
    zx1,jx1 = np.meshgrid(np.arange(0,len(wrflats)),np.linspace(1,0,47))
    xb4 = (zx1.flatten(),jx1.flatten())
    zx2,jx2 = np.meshgrid(np.arange(0,len(wrflats)),wrfeta)
    xafter = (zx2,jx2)
    zy1,jy1 = np.meshgrid(np.arange(0,len(wrflons[0])),np.linspace(1,0,47))
    yb4 = (zy1.flatten(),jy1.flatten())
    zy2,jy2 = np.meshgrid(np.arange(0,len(wrflons[0])),wrfeta)
    yafter = (zy2,jy2)
    # for each timestep in the boundary condition file
    for t in range(len(wrftimes)):
      print("Processing Lateral Boundary Conditions for :"+str(wrftimes[t]))
      if (t == 0): # open a new input file if its the first timestep
        geosfile = inputdir+'/'+'GEOSChem_'+datetime.strftime(wrftimes[0],"%Y%m%d")+'.nc'
        print('Opening file: '+geosfile)
        f = open(inputdir+'/WRFinterpbdy.cdo','w')
        f.write('gridtype = lonlat\nxsize = '+str(len(wrflons[0])+(2*bdywidth))+'\nysize = '+ str(len(wrflats)+(2*bdywidth))+'\nxfirst = '+str(wrflon_w-(bdywidth*dx))+'\nxinc = '+str(dx)+'\nyfirst ='+str(wrflat_s-(bdywidth*dy))+'\nyinc = '+str(dy))
        f.close()
        cdo.remapbil(inputdir+'/WRFinterpbdy.cdo',input=geosfile,output=geosfile+'.interpbdy')
        geos_in = Dataset(geosfile+'.interpbdy')
        geosco2 = geos_in.variables['co2'][:]
      else:
        if (wrftimes[t].date() > wrftimes[t-1].date()): # also open a new input file if its a new day
          geosfile = inputdir+'/'+'GEOSChem_'+datetime.strftime(wrftimes[t],"%Y%m%d")+'.nc'
          print('Opening file: '+geosfile)
          f = open(inputdir+'/WRFinterpbdy.cdo','w')
          f.write('gridtype = lonlat\nxsize = '+str(len(wrflons[0])+(2*bdywidth))+'\nysize = '+ str(len(wrflats)+(2*bdywidth))+'\nxfirst = '+str(wrflon_w-(bdywidth*dx))+'\nxinc = '+str(dx)+'\nyfirst ='+str(wrflat_s-(bdywidth*dy))+'\nyinc = '+str(dy))
          f.close()
          cdo.remapbil(inputdir+'/WRFinterpbdy.cdo',input=geosfile,output=geosfile+'.interpbdy')
          geos_in = Dataset(geosfile+'.interpbdy')
          geosco2 = geos_in.variables['co2'][:]
      #geos_t = date2index(wrftimes[t],geos_in.variables['time']) # get the index of the matching timestep
      geos_t = t-(int(t/24)*24)
      # do west side XS
      for i in range(bdywidth):
        regridded = interpolate.griddata(xb4,geosco2[geos_t,:,bdywidth:-bdywidth,i].flatten(),xafter,method='linear')
        CO2_XS[t,i,:,:] = regridded
      # do east side XE
      for i in range(len(wrflons)-bdywidth+1,len(wrflons)+1):
        regridded = interpolate.griddata(xb4,geosco2[geos_t,:,bdywidth:-bdywidth,i].flatten(),xafter,method='linear')
        CO2_XE[t,i-(len(wrflons)-bdywidth+1),:,:] = regridded
      # do south side YS
      for i in range(bdywidth):
        regridded = interpolate.griddata(yb4,geosco2[geos_t,:,i,bdywidth:-bdywidth].flatten(),yafter,method='linear')
        CO2_YS[t,i,:,:] = regridded
      # do north side YE
      for i in range(len(wrflats)-bdywidth+1,len(wrflats)+1):
        regridded = interpolate.griddata(yb4,geosco2[geos_t,:,i,bdywidth:-bdywidth].flatten(),yafter,method='linear')
        CO2_YE[t,i-(len(wrflats)-bdywidth+1),:,:] = regridded
    # write boundary conditions to file
    print(" * Writing Boundary Conditions to File * ")
    wrfbdy.variables['CO2_ANT_BXS'][:] = CO2_XS 
    wrfbdy.variables['CO2_ANT_BXE'][:] = CO2_XE 
    wrfbdy.variables['CO2_ANT_BYS'][:] = CO2_YS 
    wrfbdy.variables['CO2_ANT_BYE'][:] = CO2_YE 
    # calculate tendencies
    nsecs = (wrftimes[1]-wrftimes[0]).total_seconds()
    for t in range(len(wrftimes)-1):
      wrfbdy.variables['CO2_ANT_BTXS'][t] = (CO2_XS[t+1,:,:,:] - CO2_XS[t,:,:,:])/nsecs
      wrfbdy.variables['CO2_ANT_BTXE'][t] = (CO2_XE[t+1,:,:,:] - CO2_XE[t,:,:,:])/nsecs
      wrfbdy.variables['CO2_ANT_BTYS'][t] = (CO2_YS[t+1,:,:,:] - CO2_YS[t,:,:,:])/nsecs
      wrfbdy.variables['CO2_ANT_BTYE'][t] = (CO2_YE[t+1,:,:,:] - CO2_YE[t,:,:,:])/nsecs
    # for now, I'm assuming the last tendency is same as next to last, not sure how to do this correctly...
    t = len(wrftimes)-1
    wrfbdy.variables['CO2_ANT_BTXS'][t] = (CO2_XS[t,:,:,:] - CO2_XS[t-1,:,:,:])/nsecs
    wrfbdy.variables['CO2_ANT_BTXE'][t] = (CO2_XE[t,:,:,:] - CO2_XE[t-1,:,:,:])/nsecs
    wrfbdy.variables['CO2_ANT_BTYS'][t] = (CO2_YS[t,:,:,:] - CO2_YS[t-1,:,:,:])/nsecs
    wrfbdy.variables['CO2_ANT_BTYE'][t] = (CO2_YE[t,:,:,:] - CO2_YE[t-1,:,:,:])/nsecs
    print(" * GEOS2WRF COMPLETE * ") 
  else:
    print(" * Domain is nest: * ") 
    print(" * Lateral Boundary Conditions Not Generated* ") 
    print(" * GEOS2WRF COMPLETE * ") 
    sys.exit()

if __name__ == "__main__":
  if len(sys.argv) !=4:
    sys.exit('Usage: %s inputdir workdir dXX' % sys.argv[0])
  geos2wrf(sys.argv[1],sys.argv[2],sys.argv[3])
