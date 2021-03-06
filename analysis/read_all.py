#!/usr/bin/env python
# read_all.py
# use netCDF4-python to read in multiple WRF-CO2 files for analysis
# C. Martin
from netCDF4 import Dataset,MFDataset
from datetime import datetime, timedelta
import os.path
import numpy as np
rootoutput = '/homes/cmart90/WRF-CO2/output'
startday = '20160328'; startday=datetime.strptime(startday,'%Y%m%d')
endday = '20160503'; endday=datetime.strptime(endday,'%Y%m%d')
files = []
daynow = startday
while daynow <= endday:
  fname = rootoutput+'/'+daynow.strftime('%Y%m%d')+'/wrfout_d01_'+daynow.strftime('%Y-%m-%d_%H:%M:%S')
  if os.path.isfile(fname):
     files.append(fname)
  daynow = daynow + timedelta(days=1)

#data = MFDataset(files)

# plots?
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap

# lats and lons
#wrflats = data.variables['XLAT'][0][:]
#wrflons = data.variables['XLONG'][0][:]
#times = data.variables['Times'][:]

from wrfco2.plots import findpoint
UMD = [38.99,-76.94]

#times = [''.join(times[t,:]) for t in range(len(times))]
#times = [datetime.strptime(times[t],"%Y-%m-%d_%H:%M:%S") for t in range(len(times))]
#timesofday = [datetime.strptime(times[t].strftime("%H%M%S"),"%H%M%S") for t in range(len(times))]
#doy = [times[t].strftime("%j/365.") for t in range(len(times))] # day of year for dot color

#co2 = data.variables['CO2_ANT'][:,0,UMDy,UMDx]
fig1 = plt.figure()
ax1 = fig1.add_subplot(111)
a = 0

co2avearr = []

for file in files:
  print(file)
  data = Dataset(file)
  wrflats = data.variables['XLAT'][0][:]
  wrflons = data.variables['XLONG'][0][:]
  UMDx,UMDy = findpoint(wrflats,wrflons,UMD[0],UMD[1])
  times = data.variables['Times'][:]
  times = [''.join(times[t,:]) for t in range(len(times))]
  times = [datetime.strptime(times[t],"%Y-%m-%d_%H:%M:%S") for t in range(len(times))]
  timesofday = [datetime.strptime(times[t].strftime("%H%M%S"),"%H%M%S") for t in range(len(times))]
  doy = [times[t].strftime("%j") for t in range(len(times))] # day of year for dot color
  co2 = data.variables['CO2_ANT'][:,0,UMDy,UMDx]
  plt.plot(timesofday[:-1],co2[:-1],color=str(a),linestyle='-',markersize=5,label=doy[0])
  a = a + 0.02
  co2avearr.append(co2)

co2avearr = np.array(co2avearr)
co2ave = []
for i in range(len(co2avearr[0,:])):
  co2ave.append(np.mean(co2avearr[:,i]))

#plt.legend()
plt.plot(timesofday[:-1],co2ave[:-1],color='red')
ax1.set_xlabel('UTC Time of Day')
ax1.set_ylabel('CO2 concentration (ppm)')
ax1.set_title('College Park, MD WRF-CO2 - '+startday.strftime('%Y-%m-%d')+' - '+endday.strftime('%Y-%m-%d'))
import matplotlib.dates as mdates
ax1.xaxis.set_major_formatter(mdates.DateFormatter('%H:%M'))

plt.show()
