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
ncfp = '/data/soodal_data/20210322_residuals/softcal/'

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

        for iy in range(dim[1]):
            for jx in range(dim[2]):
                if not np.isfinite(ecf[iy, jx]):
                    awrf[:, iy, jx] = np.nan

        wavel = nc['/Data Fields/Wavelengths'][:]
        wavel = wavel[wavel > 0]

        relativeResiduals = awrf/fitspec
        print(relativeResiduals.shape)

        y = range(dim[1]) # spatial
        x = range(dim[2]) # image

        print('dim:', dim)

        z = wavel[0:198]

        xv, yv = np.meshgrid(x, y, indexing='xy')
        zv, yv = np.meshgrid(z, y, indexing='xy')

        relativeResiduals_imagemedian = np.nanmedian(relativeResiduals, axis=2)
        relativeResiduals_imagemedian = relativeResiduals_imagemedian.T
        relativeResiduals_imagemedian = relativeResiduals_imagemedian[:, 0:198]
        print(relativeResiduals_imagemedian.shape)
        print(zv.shape)
        print(yv.shape)
        fig, ax = plt.subplots()
        c1 = ax.contourf(zv, yv, relativeResiduals_imagemedian*100, 
                np.arange(100),
                cmap=plt.cm.jet) 
        ax.text(0.02, 0.95, yyyy+mm+dd+'T'+hh+mi+'Z', fontsize=10, transform = ax.transAxes)
        ax.set_title(yyyy+mm+dd+ '_'+hh+mi)

        cbar = plt.colorbar(c1)

        outfp = './fig/relativeResiduals_imagemedian/'
        if not os.path.exists(outfp):
            os.mkdir(outfp)
            
        outfn = yyyy + mm + dd + hh + mi + '_awrf_imagemedian.png'
        fig.savefig(outfp + outfn)
        plt.close()

        # for iw in range(dim[2]):
            # for jy in range(3):
                # mean = np.nanmean(awrf[(jy+1)*100, :, iw])
                # median = np.nanmedian(awrf[(jy+1)*100, :, iw])
                # print((jy+1)*100, 'mean', mean, 'median', median)
            # i03 = str(iw).zfill(3)
            # print(i03)

            # plt.contourf(xv, yv, awrf[:, :, iw], 
                    # np.array([-1, -0.9, -0.8, -0.7, -0.6, -0.5, -0.4,- 0.3, -0.2, 
                        # -0.1, 0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0])/30., 
                    # cmap=plt.cm.bwr) 
                    # # vmin=0.0, vmax=1.0)
            # plt.title('wavelength index='+i03)
            # plt.colorbar()

            # os.system('scp -P18742 ' + outfp+outfn + ' soodal@164.125.38.179:/mnt/d/fig/softcal &')

