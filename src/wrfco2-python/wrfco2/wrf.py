import logging
import numpy as np
import os

log = logging.getLogger(__name__)

def getMet(model,date,inittime,outdir,length='All'):
#  """ Downloads meteorology to <outdir> from
#  a specified model <model> and for a specified
#  initialization time <date><inittime> for WRF
#  """
    log = getFnLog()
    import urllib2
    import ssl # skip ssl vertification, fixes error
    gcontext = ssl.SSLContext(ssl.PROTOCOL_TLSv1) # fixes error
#  try:
    if model == 'HRRRv2_Utah':
      if length == 'All':
        length = 24
      inittime = 0
      # HRRR (at least in the Utah archive) only has the initialization times 
      filename_template = 'hrrr.t{0:02}z.wrfprsf{1:02}.grib2'
      for Hr in range(0,length,1):
        filename=filename_template.format(Hr,inittime)
        log.info('Downloading {0} ...'.format(filename))
        log.info('Downloading file: '+date.strftime(\
                'https://api.mesowest.utah.edu/archive/%Y%m%d/models/hrrr/')\
                +filename)
        f = urllib2.urlopen(date.strftime(\
                'https://api.mesowest.utah.edu/archive/%Y%m%d/models/hrrr/')\
                +filename,context=gcontext)
        print(date.strftime(\
                'https://api.mesowest.utah.edu/archive/%Y%m%d/models/hrrr/')\
                +filename)
        data = f.read()
        with open(outdir+'/'+filename,'wb') as output:
            output.write(data)
    elif model == 'HRRRv2':

      if length == 'All':
        length = 24
      inittime = 0 # lets just use analysis
      filename_template = 'hrrr.t{0:02}z.wrfprsf{1:02}.grib2'
      for Hr in range(0,length,1):
        filename=filename_template.format(Hr,inittime)
        log.info('Downloading {0} ...'.format(filename))
        log.info('Downloading file: '+date.strftime(\
                'http://www.ftp.ncep.noaa.gov/data/nccf/nonoperational/com/hrrr/prod/hrrr.%Y%m%d/')\
                +filename)
        f = urllib2.urlopen(date.strftime(\
                'http://www.ftp.ncep.noaa.gov/data/nccf/nonoperational/com/hrrr/prod/hrrr.%Y%m%d/')\
                +filename)
        data = f.read()
        with open(outdir+'/'+filename,'wb') as output:
           output.write(data)
      # also download 0z of today to make WRF happy
      filename = filename_template.format(0,0)
      from datetime import timedelta
      date = date+timedelta(days=1)
      log.info('Downloading {0} ...'.format(filename))
      log.info('Downloading file: '+date.strftime(\
              'http://www.ftp.ncep.noaa.gov/data/nccf/nonoperational/com/hrrr/prod/hrrr.%Y%m%d/')\
              +filename)
      f = urllib2.urlopen(date.strftime(\
              'http://www.ftp.ncep.noaa.gov/data/nccf/nonoperational/com/hrrr/prod/hrrr.%Y%m%d/')\
              +filename,context=gcontext)
      data = f.read()
      with open(outdir+'/hrrr.t24z.wrfprsf00.grib2','wb') as output:
          output.write(data)

    else:
      log.warn('Model not currently supported, yell at C. Martin to code more!')

def linkWPS(dest,WPSdir):
   """ Link WPS programs to the working directory """
   log = getFnLog()
   files = ['geogrid.exe',
            'link_grib.csh',
            'metgrid.exe',
            'ungrib.exe',
            'Vtable',]
   for f in files:
        filename=WPSdir+'/'+f
        try:
          os.symlink(filename,dest+'/'+f)
        except:
          pass
   log.info('Linked WPS files to working directory')

def linkWRF(dest,WRFdir):
   """ Link WRF files to the working directory """
   log = getFnLog()
   from glob import glob
   files = glob(WRFdir+'/*')

   for f in files:
      try:
        os.symlink(f,dest+'/'+os.path.basename(f))
      except:
        pass

   log.info('Linked WRF files to working directory')

def writeNml(configdir,dest,startdate,enddate):
   """ Reads in a template for both namelist.wps and namelist.input,
   replaces placeholder variables with values needed for this run """
   log = getFnLog()

   # namelist.wps 
   infilename = configdir+'/namelist.wps'
   log.info('Creating namelist.wps from template in {0}'.format(infilename))

   infile = open(infilename)
   outfile = open(dest+'/namelist.wps', 'w')

   replacements = {
        '{startdate}':startdate.strftime('%Y-%m-%d_%H:%M:%S'),
        '{enddate}':enddate.strftime('%Y-%m-%d_%H:%M:%S')}
   for line in infile:
       for src, target in replacements.iteritems():
           line = line.replace(src,target)
       outfile.write(line)
   infile.close()
   outfile.close()

   # namelist.input 
   infilename = configdir+'/namelist.input'
   log.info('Creating namelist.input from template in {0}'.format(infilename))

   infile = open(infilename)
   outfile = open(dest+'/namelist.input', 'w')

   replacements = {
        '{start_yr}':str(startdate.year),
        '{start_mn}':str(startdate.month),
        '{start_dy}':str(startdate.day),
        '{start_hh}':str(startdate.hour),
        '{start_mm}':str(startdate.minute),
        '{start_ss}':str(startdate.second),
        '{end_yr}':str(enddate.year),
        '{end_mn}':str(enddate.month),
        '{end_dy}':str(enddate.day),
        '{end_hh}':str(enddate.hour),
        '{end_mm}':str(enddate.minute),
        '{end_ss}':str(enddate.second),
        }
   for line in infile:
       for src, target in replacements.iteritems():
           line = line.replace(src,target)
       outfile.write(line)
   infile.close()
   outfile.close()

# below from T. Sluka
def getFnLog():
    '''Get a logger named specifically for the function it was called from.'''
    import inspect
    return logging.getLogger(__name__+'.'+inspect.stack()[1][3]+'()')
