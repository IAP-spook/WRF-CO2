from setuptools import setup

setup(name='wrfco2',
      version='0.1',
      description='WRF-CO2 Python Framework Package',
      url='http://aosc.umd.edu/~cmartin/wrfco2',
      author='Cory Martin',
      author_email='cmart90@umd.edu',
      license='MIT',
      packages=['wrfco2'],
      install_requires=[
          'logging','numpy','netCDF4','matplotlib','scipy','cdo','Basemap'
          ],
      zip_safe=False)
