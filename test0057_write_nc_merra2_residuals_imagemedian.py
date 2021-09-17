import pandas as pd
import numpy as np
import sys

import os
import matplotlib.pyplot as plt

from netCDF4 import Dataset

import datetime

daterange = pd.date_range(datetime.datetime(2021, 3, 22, 0, 45), 
        datetime.datetime(2021, 3, 22, 0, 45), freq='H')

# ncfp = '/data2/L2_GEMS/val_1008/softcal/'
ncfp = '/data/private/soodal/20210322_residuals/residuals/'

import glob



for dt in daterange:
    # select file from the date
    yyyy = str(dt.year).zfill(4)
    mm = str(dt.month).zfill(2)
    dd = str(dt.day).zfill(2)
    hh = str(dt.hour).zfill(2)
    mi = str(dt.minute).zfill(2)
    dtfl = glob.glob(ncfp + 'GK2B_GEMS_L2_O3P_' + yyyy + mm + dd + '_' + hh + mi +
        '_winliminit300_prec000000_climML_b4x4_o2021*.nc')
    print(dtfl)

    if len(dtfl) > 0:
        
        dtfl.sort()
        ncfn = dtfl[-1]

        # read data
        nc = Dataset(ncfn)
        print(ncfn)

        All_wavelength_ResidualsOfFit = nc['/Data Fields/AllWavelResidualsOfFit'][:]

        sza = nc['/Geolocation Fields/SolarZenithAngle']
        vza = nc['/Geolocation Fields/ViewingZenithAngle']
        raa = nc['/Geolocation Fields/RelativeAzimuthAngle']
        lat = nc['/Geolocation Fields/Latitude'][:]
        lon = nc['/Geolocation Fields/Longitude'][:]

        fitspec = nc['/Data Fields/Fitspec'][:]
        simrad = nc['/Data Fields/SimulatedRadiances'][:]

        # get EffectiveCloudFractionUV
        ecf = nc['/Data Fields/EffectiveCloudFractionUV'][:]

        # cld filtering
        ecf = np.where(ecf < 0, np.nan, 
                (np.where(ecf > 0.20, np.nan, ecf)))

        awrf = np.where(All_wavelength_ResidualsOfFit < -1E29, np.nan, All_wavelength_ResidualsOfFit)
        fitspec = np.where(fitspec < -1E29, np.nan, fitspec)

        dim = awrf.shape
        print('awrf.shape: ', dim)

        fitspec_dim = fitspec.shape
        print('fitspec.shape: ', fitspec_dim)

        # for iy in range(dim[1]):
            # for jx in range(dim[2]):
                # if not np.isfinite(ecf[iy, jx]):
                    # awrf[:, iy, jx] = np.nan

        wavel = nc['/Data Fields/Wavelengths'][:]
        wavel = wavel[wavel > 0]

        residuals_imagemedian = np.nanmedian(awrf, axis=2)
        print(residuals_imagemedian.shape)

        # write residuals
        residual_ncfn = ('/data/private/soodal/GEMS_MERRA2_prof_rtm_residuals/GK2B_GEMS_L1C_residuals_' + yyyy + mm + dd 
                        + 'T' + hh + mi+ 'Z_imagemedian.nc')

        try: ncfile.close(residual_ncfn)  # just to be safe, make sure dataset is not already open.
        except: pass



        # write residual data file
        ncfile = Dataset(residual_ncfn, mode='w', format='NETCDF4_CLASSIC')

        image_dim = ncfile.createDimension('image', dim[2])     # image axis(longitudinal)
        spatial_dim = ncfile.createDimension('spatial', dim[1])    # spatial axis(latitudinal)
        wavel_dim = ncfile.createDimension('wavelengths', 204) # wavelength axis
        for dim_item in ncfile.dimensions.items():
                print(dim_item)
        
        residuals = ncfile.createVariable('residuals', np.float64, 
                ('wavelengths', 'spatial', 'image'), fill_value=-1.E30)
        residuals.long_name = 'residuals from VLIDORT run with MERRA2 ozone profile'

        residuals[:, :, :] = -1.E30
        for i in range(dim[2]):
            print('image: ', i)
            residuals[:, :, i] = residuals_imagemedian
        # residuals.units = 
        ncfile.close()

        # os.system('scp -P18742 ' + outfp+outfn + ' soodal@164.125.38.179:/mnt/d/fig/softcal &')

