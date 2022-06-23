import os
import numpy as np
import datetime
import pandas as pd
from netCDF4 import Dataset

daterange = pd.date_range(datetime.datetime(2021, 1, 1), 
        datetime.datetime(2021, 12, 31), freq='D')
for d in daterange:
    yyyy = d.strftime('%Y')
    mm = d.strftime('%m')
    dd = d.strftime('%d')
# read 
    fn = f'/data/WORKED/TSIS_IRR/GK2_GEMS_TSIS-1_{yyyy}{mm}{dd}.nc'
    rootgrp = Dataset(fn)

    ipv = rootgrp['image_pixel_values'][:]
    ipv_shape = ipv.shape
    print(ipv_shape)

    wvl = rootgrp['wavelength'][:]
    wvl_shape = wvl.shape
    print(wvl_shape)
    bpm = rootgrp['bad_pixel_mask'][:]
    bpm_shape = bpm.shape
    print(bpm_shape)

    rootgrp.close()

# calculate
    new_shape = (1033, 512)
    ipv_4x4 = np.zeros(new_shape, dtype=np.float32)
    wvl_4x4 = np.zeros(new_shape, dtype=np.float32)
    bpm_4x4 = np.zeros(new_shape, dtype=np.int16)


    for iy in range(new_shape[1]):
        ipv_4x4[:, iy] = np.mean(ipv[:, iy*4:(iy+1)*4], axis=1)
        wvl_4x4[:, iy] = np.mean(wvl[:, iy*4:(iy+1)*4], axis=1)


# write

    outfp = '/data/WORKED/TSIS_IRR/4x4'
    if not os.path.isdir(outfp):
        os.makedirs(outfp)

    ncfile = Dataset(outfp + '/' + f'GK2_GEMS_TSIS-1_{yyyy}{mm}{dd}_BIN4x4.nc', 
            mode='w', 
            format='NETCDF4_CLASSIC')

    image_y = ncfile.createDimension('dim_image_y', 512) # spatial axis
    image_band = ncfile.createDimension('dim_image_band', 1033) # wavelength axis

    for dim in ncfile.dimensions.items():
        print(dim)

    bad_pixel_mask = ncfile.createVariable('bad_pixel_mask', np.int16, ('dim_image_band', 'dim_image_y'))
    bad_pixel_mask.FillValue = "-999s"
    bad_pixel_mask[:] = bpm_4x4[:]


    wavelength_bfcal = ncfile.createVariable('wavelength_bfcal', np.float32, ('dim_image_band', 'dim_image_y'))
    wavelength_bfcal.units = '[nm]'
    wavelength_bfcal[:] = wvl_4x4[:]

    wavelength = ncfile.createVariable('wavelength', np.float32, ('dim_image_band', 'dim_image_y'))
    wavelength.units = '[nm]'
    wavelength[:] = wvl_4x4[:]

    image_pixel_values = ncfile.createVariable('image_pixel_values', np.float32, ('dim_image_band', 'dim_image_y'))
    image_pixel_values.FillValue = '-999f'
    image_pixel_values.units = '[W/cm2/cm/sr]'
    image_pixel_values[:] = ipv_4x4[:]

    ncfile.close()

