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


# o3p file path setting
ncfp = '/data2/L2_GEMS/val_1008/softcal/'

projects = ['model_300_310', 'model_300_340', 'model_300_290']

subs = ['raw']

for sub in subs:
    fp1 = '/data/private/soodal/merra2_residual/softcal/' + sub + '/' + projects[0] # 300_310
    fp2 = '/data/private/soodal/merra2_residual/softcal/' + sub + '/' + projects[1] # 300_340
    fp3 = '/data/private/soodal/merra2_residual/softcal/' + sub + '/' + projects[2] # not applied
    fp4 = '/data2/L2_GEMS/nier_L1C' # not applied



    for dt in daterange:
        # select file from the date
        yyyy = str(dt.year).zfill(4)
        mm = str(dt.month).zfill(2)
        dd = str(dt.day).zfill(2)
        hh = str(dt.hour).zfill(2)
        mi = str(dt.minute).zfill(2)

        fl1 = glob.glob(fp1 + '/GK2_GEMS_L2_O3P_' + yyyy + mm + dd + '_' + hh + mi +
                '_winliminit300' + '_prec000000_climML_b4x4_o2021*.nc')

        fl2 = glob.glob(fp2 + '/GK2_GEMS_L2_O3P_' + yyyy + mm + dd + '_' + hh + mi +
                '_winliminit300' + '_prec000000_climML_b4x4_o2021*.nc')

        fl3 = glob.glob(fp3 + '/GK2_GEMS_L2_O3P_' + yyyy + mm + dd + '_' + hh + mi +
                '_winliminit300' + '_prec000000_climML_b4x4_o2021*.nc')

        fl4 = glob.glob(fp4 + '/GK2B_GEMS_L2_O3P_' + yyyy + mm + dd + '_' + hh + mi +
                '_4x4.nc4')

        print(fp1)
        print(fp2)
        print(fp3)
        print(fp4)

        print(fl1)
        print(fl2)
        print(fl3)
        print(fl4)

        if len(fl1) > 0 and len(fl2) > 0 and len(fl3) > 0:
            
            fl1.sort()
            fn1 = fl1[-1]
            fl2.sort()
            fn2 = fl2[-1]
            fl3.sort()
            fn3 = fl3[-1]
            fl4.sort()
            fn4 = fl4[-1]

            print(fn1)
            print(fn2)
            print(fn3)
            print(fn4)

            # read data
            nc1 = Dataset(fn1)
            nc2 = Dataset(fn2)
            nc3 = Dataset(fn3)
            nc4 = Dataset(fn4)

            # All_wavelength_ResidualsOfFit1 = nc1['/Data Fields/AllWavelResidualsOfFit'][:]
            # All_wavelength_ResidualsOfFit2 = nc1['/Data Fields/AllWavelResidualsOfFit'][:]
            # All_wavelength_ResidualsOfFit3 = nc1['/Data Fields/AllWavelResidualsOfFit'][:]

            o3_1 = nc1['/Data Fields/O3'][:]
            o3_2 = nc2['/Data Fields/O3'][:]
            o3_3 = nc3['/Data Fields/O3'][:]
            o3_4 = nc4['/Data Fields/O3'][:]

            o3ap_1 = nc1['/Data Fields/O3Apriori'][:]
            o3ap_2 = nc2['/Data Fields/O3Apriori'][:]
            o3ap_3 = nc3['/Data Fields/O3Apriori'][:]
            o3ap_4 = nc4['/Data Fields/O3Apriori'][:]

            p_1 = nc1['/Geolocation Fields/Pressure'][:]
            p_2 = nc2['/Geolocation Fields/Pressure'][:]
            p_3 = nc3['/Geolocation Fields/Pressure'][:]
            p_4 = nc4['/Geolocation Fields/Pressure'][:]

            sza1 = nc1['/Geolocation Fields/SolarZenithAngle'][:]
            vza1 = nc1['/Geolocation Fields/ViewingZenithAngle'][:]
            raa1 = nc1['/Geolocation Fields/RelativeAzimuthAngle'][:]
            lat1 = nc1['/Geolocation Fields/Latitude'][:]
            lon1 = nc1['/Geolocation Fields/Longitude'][:]

            sza2 = nc2['/Geolocation Fields/SolarZenithAngle'][:]
            vza2 = nc2['/Geolocation Fields/ViewingZenithAngle'][:]
            raa2 = nc2['/Geolocation Fields/RelativeAzimuthAngle'][:]
            lat2 = nc2['/Geolocation Fields/Latitude'][:]
            lon2 = nc2['/Geolocation Fields/Longitude'][:]

            sza3 = nc3['/Geolocation Fields/SolarZenithAngle'][:]
            vza3 = nc3['/Geolocation Fields/ViewingZenithAngle'][:]
            raa3 = nc3['/Geolocation Fields/RelativeAzimuthAngle'][:]
            lat3 = nc3['/Geolocation Fields/Latitude'][:]
            lon3 = nc3['/Geolocation Fields/Longitude'][:]

            sza4 = nc4['/Geolocation Fields/SolarZenithAngle'][:]
            vza4 = nc4['/Geolocation Fields/ViewingZenithAngle'][:]
            raa4 = nc4['/Geolocation Fields/RelativeAzimuthAngle'][:]
            lat4 = nc4['/Geolocation Fields/Latitude'][:]
            lon4 = nc4['/Geolocation Fields/Longitude'][:]

            fitspec1 = nc1['/Data Fields/Fitspec'][:]
            simrad1 = nc1['/Data Fields/SimulatedRadiances'][:]

            fitspec2 = nc2['/Data Fields/Fitspec'][:]
            simrad2 = nc2['/Data Fields/SimulatedRadiances'][:]

            fitspec3 = nc3['/Data Fields/Fitspec'][:]
            simrad3 = nc3['/Data Fields/SimulatedRadiances'][:]

            # fitspec4 = nc4['/Data Fields/Fitspec'][:]
            # simrad4 = nc4['/Data Fields/SimulatedRadiances'][:]

            # get EffectiveCloudFractionUV
            ecf1 = nc1['/Data Fields/EffectiveCloudFractionUV'][:]
            ecf2 = nc2['/Data Fields/EffectiveCloudFractionUV'][:]
            ecf3 = nc3['/Data Fields/EffectiveCloudFractionUV'][:]
            ecf4 = nc4['/Data Fields/EffectiveCloudFractionUV'][:]

            # cld filtering
            ecf1 = np.where(ecf1 < 0, np.nan, 
                    (np.where(ecf1 > 0.20, np.nan, ecf1)))
            ecf2 = np.where(ecf2 < 0, np.nan, 
                    (np.where(ecf2 > 0.20, np.nan, ecf2)))
            ecf3 = np.where(ecf3 < 0, np.nan, 
                    (np.where(ecf3 > 0.20, np.nan, ecf3)))
            ecf4 = np.where(ecf4 < 0, np.nan, 
                    (np.where(ecf4 > 0.20, np.nan, ecf4)))


            # awrf1 = np.where(All_wavelength_ResidualsOfFit1 < -1E29, np.nan, All_wavelength_ResidualsOfFit1)
            # fitspec1 = np.where(fitspec1 < -1E29, np.nan, fitspec1)

            # awrf2 = np.where(All_wavelength_ResidualsOfFit2 < -1E29, np.nan, All_wavelength_ResidualsOfFit2)
            # fitspec2 = np.where(fitspec2 < -1E29, np.nan, fitspec2)

            # awrf3 = np.where(All_wavelength_ResidualsOfFit3 < -1E29, np.nan, All_wavelength_ResidualsOfFit3)
            # fitspec3 = np.where(fitspec3 < -1E29, np.nan, fitspec3)

            dim1 = simrad1.shape
            print('simrad1.shape: ', dim1)
            dim2 = simrad2.shape
            print('simrad2.shape: ', dim2)
            dim3 = simrad3.shape
            print('simrad3.shape: ', dim3)
            # dim4 = simrad4.shape
            # print('simrad4.shape: ', dim4)

            fitspec_dim1 = fitspec1.shape
            print('fitspec1.shape: ', fitspec_dim1)
            fitspec_dim2 = fitspec2.shape
            print('fitspec2.shape: ', fitspec_dim2)
            fitspec_dim3 = fitspec3.shape
            print('fitspec3.shape: ', fitspec_dim3)
            # fitspec_dim4 = fitspec4.shape
            # print('fitspec4.shape: ', fitspec_dim4)

            wavel1 = nc1['/Data Fields/Wavelengths'][:]
            wavel1 = wavel1[wavel1 > 0]
            wavel2 = nc2['/Data Fields/Wavelengths'][:]
            wavel2 = wavel2[wavel2 > 0]
            wavel3 = nc3['/Data Fields/Wavelengths'][:]
            wavel3 = wavel3[wavel3 > 0]
            wavel4 = nc4['/Data Fields/Wavelengths'][:]
            wavel4 = wavel4[wavel4 > 0]

            # relativeResiduals1 = awrf1/fitspec1
            # relativeResiduals2 = awrf2/fitspec2
            # relativeResiduals3 = awrf3/fitspec3
            # relativeResiduals4 = awrf4/fitspec4

            # for iy in range(dim1[1]):
                # for jx in range(dim1[2]):
                    # if not np.isfinite(ecf1[iy, jx]):
                        # awrf1[:, iy, jx] = np.nan
                        # awrf2[:, iy, jx] = np.nan
                        # awrf3[:, iy, jx] = np.nan

            for iy in range(50, dim1[1], 50):
                for jx in range(30, dim1[2], 40):
                    image_xidx = jx
                    spatial_yidx = iy

                    pix_lat = lat1[iy, jx]
                    pix_lon = lon1[iy, jx]

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

                    line1 = ax1.plot(o3_1[:, spatial_yidx, image_xidx], p_1[:-1, spatial_yidx, image_xidx], 
                            marker='o', 
                            ms=5,
                            color=skyblue_color,
                            linestyle='solid', 
                            label='Model 300-310 GEMS 310-340')
                    ax1.grid(True)
                    plt.yscale('log')
                    plt.ylim((1000, 0.1))
                    ax1.set_title('GEMS L2 O3 profile')
                    ax1.text(0.05, 0.85, 'lon='+str(pix_lon), fontsize=10, transform = ax1.transAxes)
                    ax1.text(0.05, 0.80, 'lat='+str(pix_lat), fontsize=10, transform = ax1.transAxes)
                    ax1.text(0.02, 0.95, yyyy+mm+dd+'T'+hh+mi+'Z', fontsize=10, transform = ax1.transAxes)

                    ax1.set_xlabel('O3[DU]')
                    ax1.set_ylabel('Pressure[hPa]')
                    # ax1.legend(loc=(0.05, 0.6))
                    
                    line2 = ax1.plot(o3_2[:, spatial_yidx, image_xidx], p_2[:-1, spatial_yidx, image_xidx], 
                            marker='+', 
                            ms=5,
                            color=bluishgreen_color,
                            linestyle='dotted',
                            label='Model 300-340')

                    line3 = ax1.plot(o3_3[:, spatial_yidx, image_xidx], p_3[:-1, spatial_yidx, image_xidx], 
                            marker='x', 
                            ms=5,
                            color=yellow_color,
                            linestyle='dashed',
                            label='GEMS 300-340')

                    line4 = ax1.plot(o3_4[:, spatial_yidx, image_xidx], p_4[:-1, spatial_yidx, image_xidx], 
                            marker='.', 
                            ms=2,
                            color=vermilion_color,
                            linestyle='dashed',
                            label='GEMS 310-340')

                    line5 = ax1.plot(o3ap_3[:, spatial_yidx, image_xidx], p_3[:-1, spatial_yidx, image_xidx], 
                            marker='.', 
                            ms=2,
                            color=black_color,
                            linestyle='dashdot',
                            label='Apriori')

                    ax1.legend(loc=(0.5, 0.8))

                    outfp = './fig/GEMS_O3P_profile_model_300_310/' + sub + '/'
                    if not os.path.exists(outfp):
                        os.mkdir(outfp)
                        
                    outfn = (yyyy + mm + dd + hh + mi + 
                            '_o3profile_x'+str(image_xidx).zfill(3) + 
                            '_y'+str(spatial_yidx).zfill(3) + '.png')
                    print(outfp + outfn)

                    fig.savefig(outfp + outfn)
                    plt.close()
