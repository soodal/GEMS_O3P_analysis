import pandas as pd
import numpy as np
import sys
import os
import matplotlib.pyplot as plt
from netCDF4 import Dataset
import datetime
import glob
import seaborn as sns

daterange = pd.date_range(datetime.datetime(2021, 3, 22, 0, 45), 
        datetime.datetime(2021, 3, 22, 0, 45), freq='H')

ncfp = '/data2/L2_GEMS/val_1008/softcal/'
ncfp = '/data/private/soodal/20210322_residuals/residuals/'


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
        print(ncfn)

        # read data
        nc = Dataset(ncfn)

        All_wavelength_ResidualsOfFit = nc['/Data Fields/AllWavelResidualsOfFit'][:]

        sza = nc['/Geolocation Fields/SolarZenithAngle'][:]
        vza = nc['/Geolocation Fields/ViewingZenithAngle'][:]
        raa = nc['/Geolocation Fields/RelativeAzimuthAngle'][:]
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

        wavel = nc['/Data Fields/Wavelengths'][:]
        wavel = wavel[wavel > 0]

        relativeResiduals = awrf/fitspec

        for iy in range(dim[1]):
            for jx in range(dim[2]):
                if not np.isfinite(ecf[iy, jx]):
                    awrf[:, iy, jx] = np.nan

        for iy in range(50, dim[1], 50):
            for jx in range(30, dim[2], 40):
                image_xidx = jx
                spatial_yidx = iy
                pix_lat = lat[iy, jx]
                pix_lon = lon[iy, jx]

                plt.style.use('default')
                plt.rcParams['figure.figsize'] = (7, 5)
                plt.rcParams['font.size'] = 12

                black_color = '#000000'
                orange_color = '#e4a000'
                skyblue_color = '#56b4e8'
                bluishgreen_color = '#009f73'
                yellow_color = '#f0e442'
                blue_color = '#0072b1'
                vermilion_color = '#d65e00'
                redish_purploe = '#cc79a7'

                fig, ax1 = plt.subplots()
                line1 = ax1.plot(wavel[0:198], relativeResiduals[0:198, spatial_yidx, image_xidx]*100., 
                        marker='o', 
                        ms=5,
                        color=skyblue_color,
                        linestyle='solid', 
                        label='Relative Differences')
                ax1.grid(True)
                ax1.set_title('GEMS Spectral Radiance with RTM, Differences')
                ax1.text(0.05, 0.85, 'lon='+str(pix_lon), fontsize=10, transform = ax1.transAxes)
                ax1.text(0.05, 0.80, 'lat='+str(pix_lat), fontsize=10, transform = ax1.transAxes)
                ax1.text(0.02, 0.95, yyyy+mm+dd+'T'+hh+mi+'Z', fontsize=10, transform = ax1.transAxes)

                ax1.set_xlabel('Wavelengths[nm]')
                ax1.set_ylabel('Differences[%]')
                ax1.legend(loc=(0.05, 0.6))
                
                ax2 = ax1.twinx()
                line2 = ax2.plot(wavel[0:198], simrad[0:198, spatial_yidx, image_xidx], 
                        marker='+', 
                        ms=5,
                        color=bluishgreen_color,
                        linestyle='dotted',
                        label='RTM from MERRA2 profile')
                line3 = ax2.plot(wavel[0:198], fitspec[0:198, spatial_yidx, image_xidx], 
                        marker='x', 
                        ms=5,
                        color=yellow_color,
                        linestyle='dashed',
                        label='GEMS Observed Radiances')
                ax2.set_ylabel('Radiances[W cm^-2 cm^-1 sr-^1]')
                ax2.legend(loc=(0.5, 0.4))

                # lines = line1 + line2 + line3
                # labels = [l.get_label() for l in lines]
                # plt.legend(lines, labels, loc=9)

                outfp = './fig/awrf_pixel/'
                if not os.path.exists(outfp):
                    os.mkdir(outfp)
                    
                outfn = (yyyy + mm + dd + hh + mi + 
                        '_awrf_x'+str(image_xidx).zfill(3) + 
                        '_y'+str(spatial_yidx).zfill(3) + '.png')
                print(outfp + outfn)

                fig.savefig(outfp + outfn)
                plt.close()
