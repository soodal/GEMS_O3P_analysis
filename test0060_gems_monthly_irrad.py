#!/home/soodal/anaconda3/bin/python3
import pandas as pd
import numpy as np
import sys

import os
import matplotlib.pyplot as plt

from netCDF4 import Dataset

import datetime


# date 
daterange = pd.date_range(datetime.datetime(2021, 2, 15, 0, 45), 
        datetime.datetime(2021, 3, 21, 0, 45), freq='D')

# input irradiance file list
input_irr_fp = '/data2/L1C_GEMS/L1C_nier/4x4/'

import glob

def local_func():
    var = "Test"
    if 'var' in locals():
        print ('var variable exists')
    else:
        print ('var variable does not exist in the local namespace')

for dt in daterange:
    # select file from the date
    yyyy = str(dt.year).zfill(4)
    mm = str(dt.month).zfill(2)
    dd = str(dt.day).zfill(2)
    hh = str(dt.hour).zfill(2)
    mi = str(dt.minute).zfill(2)
    input_fl = glob.glob(input_irr_fp + 'GK2_GEMS_IRR_' + yyyy + mm + dd + '_4x4.nc')

    if len(input_fl) > 0:
        input_fl.sort()
        ncfn = input_fl[-1]
        print(ncfn)

        # read data
        nc = Dataset(ncfn)

        bad_pixel_mask = nc['/bad_pixel_mask'][:]
        image_pixel_values = nc['/image_pixel_values'][:]
        wavelength = nc['/wavelength'][:]
        wavelength_bfcal = nc['/wavelength_bfcal'][:]
g       # print(image_pixel_values.shape)
        if 'image_pixel_values_total' in locals():
            image_pixel_values_total = np.concatenate(
                    (image_pixel_values_total, 
                        image_pixel_values.reshape(1, 1033, 512)), axis=0)
            bad_pixel_mask_total = np.concatenate(
                    (bad_pixel_mask_total, 
                        bad_pixel_mask.reshape(1, 1033, 512)), axis=0)
            wavelength_total = np.concatenate(
                    (wavelength_total, 
                        wavelength.reshape(1, 1033, 512)), axis=0)
            wavelength_bfcal_total = np.concatenate(
                    (wavelength_bfcal_total, 
                        wavelength_bfcal.reshape(1, 1033, 512)), axis=0)
        else:
            image_pixel_values_total = image_pixel_values.reshape(1, 1033, 512)
            bad_pixel_mask_total = bad_pixel_mask.reshape(1, 1033, 512)
            wavelength_total = wavelength.reshape(1, 1033, 512)
            wavelength_bfcal_total = wavelength_bfcal.reshape(1, 1033, 512)
        print(image_pixel_values_total.shape)



bad_pixel_mask_mean = np.nanmean(bad_pixel_mask_total, axis=0)
image_pixel_values_mean = np.nanmean(image_pixel_values_total, axis=0)
print(image_pixel_values_mean[100, :])
wavelength_mean = np.nanmean(wavelength_total, axis=0)
print(wavelength_mean[100, :])
wavelength_bfcal_mean = np.nanmean(wavelength_bfcal_total, axis=0)
print(wavelength_mean[100, :])

print(bad_pixel_mask.shape)
print(image_pixel_values_mean.shape)
print(wavelength_mean.shape)
print(wavelength_bfcal_mean.shape)

mean_irr_fp = '/data/private/soodal/IRR_mean/'
if not os.path.exists(mean_irr_fp):
    os.makedirs(mean_irr_fp)

# using last day 
mean_irr_fn = 'GK2_GEMS_IRR_' + yyyy + mm + dd + '_mean31_4x4.nc'

# write residual data file
ncfile = Dataset(mean_irr_fp + '/' + mean_irr_fn, mode='w', format='NETCDF4_CLASSIC')

dim_image_y = ncfile.createDimension('dim_image_y', 512)     # image axis(longitudinal)
dim_image_band = ncfile.createDimension('dim_image_band', 1033)     # image axis(longitudinal)

image_pixel_values = ncfile.createVariable('image_pixel_values', np.float32, 
        ('dim_image_band', 'dim_image_y'))
# image_pixel_values.long_name = 'image pixel values'
image_pixel_values.units = '[W/cm2/cm/sr]'
image_pixel_values[:, :] = -999
image_pixel_values.FillValue="-999f"
image_pixel_values[:, :] = image_pixel_values_mean[:, :]

wavelength = ncfile.createVariable('wavelength', np.float32, 
        ('dim_image_band', 'dim_image_y'))
# wavelength.long_name = 'wavelength'
wavelength.units = '[nm]'
# wavelength[:, :] = -1.0E30
wavelength[:, :]= wavelength_mean[:, :]

wavelength_bfcal = ncfile.createVariable('wavelength_bfcal', np.float32, 
        ('dim_image_band', 'dim_image_y'))
# wavelength_bfcal.long_name = 'wavelength_bfcal'
wavelength_bfcal.units = '[nm]'
# wavelength_bfcal[:, :] = -1.0E30
wavelength_bfcal[:, :] = wavelength_bfcal_mean[:, :]

bad_pixel_mask = ncfile.createVariable('bad_pixel_mask', np.int16, 
        ('dim_image_band', 'dim_image_y'))
bad_pixel_mask.FillValue="-999s"
# bad_pixel_mask.long_name = 'bad pixel mask'
# bad_pixel_mask[:, :] = -999
bad_pixel_mask[:, :] = bad_pixel_mask_mean[:, :]

print(image_pixel_values[100, :])
ncfile.close()

# os.system('scp -P18742 ' + outfp+outfn + ' soodal@164.125.38.179:/mnt/d/fig/softcal &')

