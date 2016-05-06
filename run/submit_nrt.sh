#!/bin/bash
#SBATCH --nodes=4
#SBATCH -t 3:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=cmart90@umd.edu

module load netcdf-fortran/intel/2015.0.3.032/intelmpi-mt/netcdf/4.3.3.1/4.4.2
module load netcdf-fortran/intel/2015.0.3.032/intelmpi-mt/netcdf/4.3.3.1/4.4.2
module load cdo

export NETCDF=$NETCDF_FORTRAN_ROOT
export LD_LIBRARY_PATH=$NETCDF_LIBDIR:$NETCDF_FORTRAN_LIBDIR:$LD_LIBRARY_PATH
export MPI_LIB=
ulimit -s unlimited
export WRF_CHEM=1
export NRTdir=/homes/cmart90/WRF-CO2

source ~cmart90/.bashrc.mine
cd $NRTdir

python $NRTdir/run/run_nrt
