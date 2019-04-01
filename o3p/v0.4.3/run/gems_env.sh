#!/bin/bash


# Find out if we have HOSTTYPE (Harvard SGIs don"t!)
if [ -"${HOSTTYPE}" = -"" ]; then
    export HOSTTYPE=`uname`
fi

# Determine the name of the machine name we are running on
 export MACHINE=`uname -a | awk "{print $2}"`

if [ ! ${FLAVOR} ]; then
   if [ $HOSTTYPE = "i386-linux" ]; then
      export FLAVOR=linux
   fi
   if [ $HOSTTYPE = "i386" ]; then
      export FLAVOR=linux
   fi
   if [ $HOSTTYPE = "SunOS" ]; then
      export FLAVOR=solaris
   fi
   if [ $HOSTTYPE = "alpha" ]; then
      export FLAVOR=alpha
   fi

   # On JOVE set the flavor of Linux to that of Intel Fortran compiler
   if [ $HOSTTYPE = "Linux" ]; then
      export FLAVOR=linux
   fi

   # On the Harvard SGIs we can have to HOSTTYPEs, depending on
   # whether it was set by the OS or by `uname` above
   if [ $HOSTTYPE = "IRIX64" -o $HOSTTYPE = "iris4d" ]; then
      export FLAVOR=sgi32
   fi

   if [ $HOSTTYPE = "x86_64" ]; then
      export FLAVOR=linux64
   fi
fi

export O3P_VER=v0.4.3
export PGEVERSION=$O3P_VER
export PGEHOME=$O3P_HOME/$O3P_VER
export LIDORTHOME=$O3P_HOME/share/lidort
export PATH=$PATH0  # initialize $PATH
export PATH=$PATH:$O3P_HOME/$O3P_VER/run
export PATH=$PATH:/opt/hdf/h4h5tools-2.2.2-linux-static/bin:$O3P_HOME/$O3P_VER/bin:.


export SHARE_HOME=$GEMS_HOME/share
export OMIUTIL=$SHARE_HOME/extlib
export DATDIR=$O3P_HOME/dat
export GEMSBIN=$PGEHOME/bin 

# Source the central toolkit set-up script; this gives access to
# binaries required during compilation
#source $OMIUTIL/toolkit/bin/linux64/pgs-dev-env.csh
source $OMIUTIL/toolkit/bin/linux64/pgs-dev-env.ksh

# Re-determine the name of the machine name we are running on
export MACHINE=`uname -a | awk "{print $2}"`

# Environment variables for this program/project
export SAOPGE=GEMS_O3P
export PGEIODIR=$PGEHOME/out
export PGSMSG=$PGEHOME/msg
export MKFDIR=$PGEHOME/make


# OMI environment variables (included for OMI compatibility)
export UTDATA=/data/mambo/tkurosu/omi/UTdata/
#export PGS_PC_INFO_FILE=`pwd`/../src/$SAOPGE.pcf
export PGS_PC_INFO_FILE=$PGEHOME/run/conf/$SAOPGE.pcf
export RADLUN=1100
export IRRLUN=1104

# Harvard SGI r=uire LD_LIBRARY_PATH, because of prebuilt 
# HDF5 shared object libraries
if [ -"$FLAVOR" = -"sgi32" ]; then
   export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$OMIUTIL/lib/$FLAVOR
fi

echo "PGEHOME $PGEHOME"  
echo "Building 'SOMIPROF.exe for" $FLAVOR 
