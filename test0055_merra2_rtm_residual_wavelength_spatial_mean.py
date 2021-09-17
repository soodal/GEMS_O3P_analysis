import pandas as pd
import numpy as np
import sys

import os
import matplotlib.pyplot as plt

from netCDF4 import Dataset

import datetime

daterange = pd.date_range(datetime.datetime(2020, 8, 1, 0, 45), 
        datetime.datetime(2020,8, 11, 00, 45), freq='H')

ncfp = '/data2/L2_GEMS/val_1008/softcal/'

import glob


awrf_total = np.zeros([512,204])

for dt in daterange:
    # select file from the date
    yyyy = str(dt.year).zfill(4)
    mm = str(dt.month).zfill(2)
    dd = str(dt.day).zfill(2)
    hh = str(dt.hour).zfill(2)
    mi = str(dt.minute).zfill(2)
    dtfl = glob.glob(ncfp + 'GK2B_GEMS_L2_O3P_' + yyyy + mm + dd + '_' + hh + mi +
        '_winliminit310_prec000000_climML_b4x4_o2021*.nc4')

    if len(dtfl) > 0:
        
        dtfl.sort()
        ncfn = dtfl[-1]
        print(ncfn)

        # read data
        nc = Dataset(ncfn)

        All_wavelength_ResidualsOfFit = nc['/Data Fields/All_wavelength_ResidualsOfFit'][:]

        sza = nc['/Geolocation Fields/SolarZenithAngle']
        vza = nc['/Geolocation Fields/ViewingZenithAngle']
        raa = nc['/Geolocation Fields/RelativeAzimuthAngle']

        df = nc['/Data Fields']

        # get EffectiveCloudFractionUV
        ecf = nc['/Data Fields/EffectiveCloudFractionUV'][:]
        # print(type(ecf))
        # print(ecf.shape)

        # cld filtering
        ecf = np.where(ecf < 0, np.nan, 
                (np.where(ecf > 0.20, np.nan, ecf)))


        # print(np.nansum(ecf))

        # All_wavelength_ResidualsOfFit[nanidx] = np.nan
        # print(np.nansum(All_wavelength_ResidualsOfFit))
        awrf = np.where(All_wavelength_ResidualsOfFit < -1E29, np.nan, All_wavelength_ResidualsOfFit)
        print(awrf.shape)

        dim = awrf.shape
        print('awrf.shape: ', dim)

        for iy in range(dim[0]):
            for jx in range(dim[1]):
                if not np.isfinite(ecf[iy, jx]):
                    awrf[iy, jx, :] = np.nan

        y = range(dim[0]) # spatial
        x = range(dim[1]) # image
        z = range(dim[2]) # wavelength
        xv, yv = np.meshgrid(x, y, indexing='xy')
        xv, yv = np.meshgrid(x, y, indexing='xy')

        zv, yv = np.meshgrid(z, y, indexing='xy')

        print(xv[:, 0])
        print(yv[:, 0])

        print(xv[0, :])
        print(yv[0, :])

        print(zv[:, 0])
        print(yv[:, 0])

        print(zv[0, :])
        print(yv[0, :])

        awrf_imagemean = np.nanmean(awrf, axis=1)
        plt.contourf(zv, yv, awrf_imagemean, 
                np.array([-1, -0.9, -0.8, -0.7, -0.6, -0.5, -0.4,- 0.3, -0.2, 
                    -0.1, 0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0])/30., 
                cmap=plt.cm.bwr) 
        plt.title(yyyy+mm+dd+ '_'+hh+mi)
        plt.colorbar()

        outfp = './fig/awrf_imagemean/'
        if not os.path.exists(outfp):
            os.mkdir(outfp)
            
        outfn = yyyy + mm + dd + hh + mi + '_awrf_imagemean.png'
        plt.savefig(outfp + outfn)
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

