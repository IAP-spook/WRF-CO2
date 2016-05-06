"""
wrfplot.py
Python module to plot WRF output easily
C. Martin - UMD - 3/2016
"""

def basicsurfmap(file,varstr,outputdir):
  """
  Creates plots for each timestep of a file
  for the specified variable and outputs them
  to a specified directory
  """
  import matplotlib
  matplotlib.use('Agg')
  import numpy as np
  import sys
  from netCDF4 import Dataset
  import matplotlib.pyplot as plt
  from mpl_toolkits.basemap import Basemap
  # read in file and vars
  f = Dataset(file,'r')
  wrfeta = f.variables['ZNU'][0][:]
  times = f.variables['Times'][:]
  wrflats = f.variables['XLAT'][0][:]
  wrflons = f.variables['XLONG'][0][:]
  var = f.variables[varstr][:]
  # four corners of domain
  print wrflats.shape
  print wrflons.shape
  print wrflons[0].shape
  wrflat_s = wrflats[0,len(wrflats)-1]
  wrflat_n = wrflats[len(wrflats)-1,len(wrflons[0])-1]
  wrflon_w = wrflons[0,0]
  wrflon_e = wrflons[len(wrflats)-1,len(wrflons[0])-1]
 
  z = 0 # assuming lowest level of model

  # set up map
  map = Basemap(projection='merc',llcrnrlon=wrflon_w,urcrnrlon=wrflon_e,llcrnrlat=wrflat_s,urcrnrlat=wrflat_n,resolution='i')
  map.drawstates()
  map.drawcounties()
  map.drawcoastlines()
  x,y = map(wrflons,wrflats)
 
  # loop through times
  for t in range(len(times)):
    timestr = ''.join(times[t,:])
    map.drawstates()
    map.drawcounties()
    map.drawcoastlines()
    plt1 = map.pcolormesh(x,y,var[t,z,:,:],vmin=np.amin(var),vmax=np.amax(var))
    colorbar = map.colorbar(plt1,"right", size="5%",pad="2%")
    colorbar.set_label(f.variables[varstr].description+' '+f.variables[varstr].units)   
    plt.title('WRF output valid: '+timestr)
    plt.savefig(outputdir+'/%03d_' % (t) +timestr+'_'+varstr+'.png')
    plt.clf()

def CO2nWind(file,outputdir):
  """
  Outputs plots of CO2 in the lowest level and 10m winds
  """
  import matplotlib
  matplotlib.use('Agg')
  import numpy as np
  import sys
  from netCDF4 import Dataset
  import matplotlib.pyplot as plt
  from mpl_toolkits.basemap import Basemap
  # read in file and vars
  f = Dataset(file,'r')
  wrfeta = f.variables['ZNU'][0][:]
  times = f.variables['Times'][:]
  wrflats = f.variables['XLAT'][0][:]
  wrflons = f.variables['XLONG'][0][:]
  var = f.variables['CO2_ANT'][:,0,:,:]
  u = f.variables['U'][:,0,:,:]
  v = f.variables['V'][:,0,:,:]
  # destagger u/v
  u = (u[:,:,:-1] + u[:,:,1:])/2.
  v = (v[:,:-1,:] + v[:,1:,:])/2.

  # four corners of domain
  wrflat_s = wrflats[0,len(wrflats)-1]
  wrflat_n = wrflats[len(wrflats)-1,len(wrflons[0])-1]
  wrflon_w = wrflons[0,0]
  wrflon_e = wrflons[len(wrflats)-1,len(wrflons[0])-1]

  z = 0 # assuming lowest level of model

  # set up map
  map = Basemap(projection='merc',llcrnrlon=wrflon_w,urcrnrlon=wrflon_e,llcrnrlat=wrflat_s,urcrnrlat=wrflat_n,resolution='i')
  map.drawstates()
  map.drawcounties()
  map.drawcoastlines()
  x,y = map(wrflons,wrflats)

  # loop through times
  for t in range(len(times)):
    timestr = ''.join(times[t,:])
    map.drawstates(color='gray',linewidth=1)
    map.drawcounties(color='white')
    map.drawcoastlines(color='gray',linewidth=1)
    plt1 = map.pcolormesh(x,y,var[t,:,:],vmin=380,vmax=450)
    #plt1 = map.pcolormesh(x,y,var[t,:,:],vmin=np.amin(var),vmax=np.amax(var))
    winds = map.barbs(x[::20,::20],y[::20,::20],u[t,::20,::20]*1.94,v[t,::20,::20]*1.94,length=6,color='white') # *1.94 to convert m/s to knots (barb convention)
    colorbar = map.colorbar(plt1,"right", size="5%",pad="2%")
    colorbar.set_label(f.variables['CO2_ANT'].description+' '+f.variables['CO2_ANT'].units)
    plt.title('WRF output valid: '+timestr)
    plt.savefig(outputdir+'/%03d_' % (t) +timestr+'_CO2_wind.png')
    plt.clf()

def findpoint(latvar,lonvar,lat0,lon0):
    '''
    Find closest point in a set of (lat,lon) points to specified point
    latvar - 2D latitude variable from an open netCDF dataset
    lonvar - 2D longitude variable from an open netCDF dataset
    lat0,lon0 - query point
    Returns iy,ix such that the square of the tunnel distance
    between (latval[it,ix],lonval[iy,ix]) and (lat0,lon0)
    is minimum.
    From Unidata python workshop
    '''
    from math import pi
    import numpy as np
    from numpy import cos,sin
    rad_factor = pi/180.0 # for trignometry, need angles in radians
    # Read latitude and longitude from file into numpy arrays
    latvals = latvar[:] * rad_factor
    lonvals = lonvar[:] * rad_factor
    ny,nx = latvals.shape
    lat0_rad = lat0 * rad_factor
    lon0_rad = lon0 * rad_factor
    # Compute numpy arrays for all values, no loops
    clat,clon = cos(latvals),cos(lonvals)
    slat,slon = sin(latvals),sin(lonvals)
    delX = cos(lat0_rad)*cos(lon0_rad) - clat*clon
    delY = cos(lat0_rad)*sin(lon0_rad) - clat*slon
    delZ = sin(lat0_rad) - slat;
    dist_sq = delX**2 + delY**2 + delZ**2
    minindex_1d = dist_sq.argmin()  # 1D index of minimum element
    iy_min,ix_min = np.unravel_index(minindex_1d, latvals.shape)
    return iy_min,ix_min

def UMDmodelVcss(modelfile,cssfile,plotdir):
    ''' Plot comparison of model run to observed data '''
    import matplotlib
    matplotlib.use('Agg')
    import numpy as np
    import sys
    from netCDF4 import Dataset
    import matplotlib.pyplot as plt
    from datetime import datetime
    # read in file and vars
    f = Dataset(modelfile,'r')
    wrfeta = f.variables['ZNU'][0][:]
    times = f.variables['Times'][:]
    wrflats = f.variables['XLAT'][0][:]
    wrflons = f.variables['XLONG'][0][:]
    var = f.variables['CO2_ANT'][:,0,:,:]
    z = 0 # assuming lowest level of model
    # College Park lat/lon
    UMD = [38.99,-76.94]
    UMDx,UMDy = findpoint(wrflats,wrflons,UMD[0],UMD[1])
   
    # read in CO2 data from SENSE format
    date,time,co2css,ch4css,h2ocss = np.genfromtxt(cssfile,missing_values=-9999.00,filling_values=np.nan,usecols=(0,1,2,3,4),unpack=True)
    co2css[co2css==-9999]=np.nan                 # above doesn't work; explicitly define here
    tim = ["%8d%06d" % (date[t],time[t]) for t in range(len(date))]
    tim = [datetime.strptime(tim[t], "%Y%m%d%H%M%S") for t in range(len(tim))]

    # convert wrf times to datetime objects
    times = [''.join(times[t,:]) for t in range(len(times))]
    times = [datetime.strptime(times[t],"%Y-%m-%d_%H:%M:%S") for t in range(len(times))]
    
    # plot data
    plt.plot_date(times,var[:,UMDy,UMDx],color='red',label='WRF')
    plt.plot_date(tim,co2css,color='black',label='LGR')
    plt.ylim([380,430])
    plt.legend()
    plt.savefig(plotdir+'/UMDLGRvWRF.png') 

def theta_co2_profiles(file,plotdir):
    ''' Plot two panels at multiple (4) time series
        one of CO2 vertical profile and one of
        theta to get PBLH / inversions
    '''
    import matplotlib
    matplotlib.use('Agg')
    import numpy as np
    import sys
    from netCDF4 import Dataset
    import matplotlib.pyplot as plt
    from datetime import datetime, timedelta
    
    # read in file and vars
    f = Dataset(file,'r')
    wrfeta = f.variables['ZNU'][:]
    times = f.variables['Times'][:]
    wrflats = f.variables['XLAT'][0][:]
    wrflons = f.variables['XLONG'][0][:]
    
    # convert wrf times to datetime objects
    times = [''.join(times[t,:]) for t in range(len(times))]
    times = [datetime.strptime(times[t],"%Y-%m-%d_%H:%M:%S") for t in range(len(times))]
    times = [times[t] - timedelta(hours=5) for t in range(len(times))] # convert to EST

    # profile for college park
    UMD = [38.99,-76.94]
    UMDx,UMDy = findpoint(wrflats,wrflons,UMD[0],UMD[1])

    co2 = f.variables['CO2_ANT'][:,:,UMDy,UMDx]
    P = f.variables['P'][:,:,UMDy,UMDx]
    T = f.variables['T'][:,:,UMDy,UMDx] # potential temperature
    T = T + 300 # WRF outputs temp like this, weird I know
    PB = f.variables['PB'][:,:,UMDy,UMDx]
    P = P + PB # pressure at each eta level is P + perturbation pressure PB
    
    #### set up figure
    # times to plot
    ts = [0,12,24,36,48,60,72,84]
    ts = [0,24,48,72]
    colors = ['red','orange','yellow','green','blue','purple','black','gray']
    colors = ['red','blue','green','black']
    fig, (ax1, ax2) = plt.subplots(1,2, sharey=True)
    for t in range(len(ts)):
		ax1.plot(co2[ts[t],:],P[ts[t],:]/100,color=colors[t],label=times[ts[t]].strftime('%m/%d %H%M LST'))
		ax2.plot(T[ts[t],:],P[ts[t],:]/100,color=colors[t],label=times[ts[t]].strftime('%m/%d %H%M LST'))
    ax1.grid(True)
    ylabels = [1000,975,950,925,900,850]
    ax1.set_ylim([150,1025])
    ax1.set_xlim([380,450])
    ax1.invert_yaxis()
    ax1.set_ylabel('Pressure (hPa)')
    ax1.set_xlabel(f.variables['CO2_ANT'].description+' '+f.variables['CO2_ANT'].units)
    ax2.set_xlim([270,320])
    #plt.yticks(ylabels,ylabels)
    ax2.set_xlabel('Potential Temperature (Kelvin)')
    ax2.grid(True)
    plt.legend()
    plt.suptitle('WRF-CO2 output for College Park, MD')
    
    plt.savefig(plotdir+'/theta_co2_profiles.png')
   
def anim_profiles(file,plotdir):
    ''' Plot animated profiles of CO2/theta with a complete timeseries of
    	surface CO2 below it
    '''
    import matplotlib
    matplotlib.use('Agg')
    import numpy as np
    import sys
    from netCDF4 import Dataset
    import matplotlib.pyplot as plt
    from datetime import datetime, timedelta
    import matplotlib.dates as mdates
    
    # read in file and vars
    f = Dataset(file,'r')
    wrfeta = f.variables['ZNU'][:]
    times = f.variables['Times'][:]
    wrflats = f.variables['XLAT'][0][:]
    wrflons = f.variables['XLONG'][0][:]
    
    # convert wrf times to datetime objects
    times = [''.join(times[t,:]) for t in range(len(times))]
    times = [datetime.strptime(times[t],"%Y-%m-%d_%H:%M:%S") for t in range(len(times))]
    times = [times[t] - timedelta(hours=5) for t in range(len(times))] # convert to EST

    # profile for college park
    UMD = [38.99,-76.94]
    UMDx,UMDy = findpoint(wrflats,wrflons,UMD[0],UMD[1])

    co2 = f.variables['CO2_ANT'][:,:,UMDy,UMDx]
    P = f.variables['P'][:,:,UMDy,UMDx] # base state pressure
    PBLH = f.variables['PBLH'][:,UMDy,UMDx] # PBL Height
    T = f.variables['T'][:,:,UMDy,UMDx] # potential temperature
    T = T + 300 # WRF outputs temp like this, weird I know
    PB = f.variables['PB'][:,:,UMDy,UMDx] # perturbation pressure
    P = P + PB # pressure at each eta level is P + perturbation pressure PB
    for t in range(len(times)):
    	timestr = times[t].strftime("%Y-%m-%d_%H:%M:%S")
    	fig = plt.figure(figsize=(11,8.5))
    	ax1 = plt.subplot2grid((4,2),(0,0),rowspan=2) # CO2 profile
    	ax2 = plt.subplot2grid((4,2),(0,1),rowspan=2,sharey=ax1) # theta profile
    	ax3 = plt.subplot2grid((4,2),(2,0),colspan=2) # surface co2 time series
    	ax4 = plt.subplot2grid((4,2),(3,0),colspan=2,sharex=ax3) # PBLH time series
    	# CO2 profile
    	ax1.plot(co2[t,:],P[t,:]/100.,color='black')
    	ax1.set_ylim([150,1025])
    	ax1.set_xlim([380,450])
    	ax1.invert_yaxis()
    	ax1.grid(True)
    	ax1.set_title('CO2 (ppm)')
    	ax1.set_ylabel('Pressure (hPa)')
    	# theta profile
    	ax2.plot(T[t,:],P[t,:]/100.,color='black')
    	ax2.set_ylim([150,1025])
    	ax2.set_xlim([270,320])
    	ax2.invert_yaxis()
    	ax2.grid(True)
    	ax2.set_title('Potential Temperature (K)')
    	# surface CO2 time series
    	ax3.plot_date(times,co2[:,0],linestyle='-',marker='',color='blue')
    	ax3.plot_date(times[t],co2[t,0],marker='o',color='red')
    	ax3.set_ylim([380,450])
    	plt.setp( ax3.get_xticklabels(), visible=False)
    	ax3.grid(True)
    	ax3.set_ylabel('Surface CO2 (ppm)')
    	# PBLH time series
    	ax4.plot_date(times,PBLH,linestyle='-',marker='',color='blue')
    	ax4.plot_date(times[t],PBLH[t],marker='o',color='red')
    	ax4.grid(True)
    	ax4.xaxis.set_major_formatter(mdates.DateFormatter('%H:%M'))
    	ax4.set_ylabel('PBL height (m)')
    	
    	# save figure and restart loop
    	#fig.autofmt_xdate()
    	plt.suptitle('WRF-CO2 output valid: '+timestr+' LST for College Park, MD')
    	plt.savefig(plotdir+'/%03d_' % (t) +timestr+'_anim_profiles.png')
    	plt.clf()
