import os
import numpy as np
import datetime
import pandas as pd
from netCDF4 import Dataset

daterange = pd.date_range(datetime.datetime(2022, 3, 1), 
        datetime.datetime(2022, 3, 2), freq='D')
for d in daterange:
    yyyy = d.strftime('%Y')
    mm = d.strftime('%m')
    dd = d.strftime('%d')
# read 
    fn = f'/data/WORKED/TSIS_IRR/GK2_GEMS_TSIS-1_{yyyy}{mm}{dd}.nc'
    rootgrp = Dataset(fn)

    try: 
        created_by = rootgrp.Created_by
        do_write_created_by = True 
    except AttributeError:
        do_write_created_by = False

    try: 
        created_at = rootgrp.Created_at
        do_write_created_at = True
    except AttributeError:
        do_write_created_at = False

    ipv = rootgrp['image_pixel_values'][:]
    ipv_shape = ipv.shape

    wvl = rootgrp['wavelength'][:]
    wvl_shape = wvl.shape

    bpm = rootgrp['bad_pixel_mask'][:]
    bpm_shape = bpm.shape
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

    outfn = outfp + '/' + f'GK2_GEMS_TSIS-1_{yyyy}{mm}{dd}_BIN4x4.nc'
    ncfile = Dataset(outfn, 
            mode='w', 
            format='NETCDF4_CLASSIC')

    image_y = ncfile.createDimension('dim_image_y', 512) # spatial axis
    image_band = ncfile.createDimension('dim_image_band', 1033) # wavelength axis

    # for dim in ncfile.dimensions.items():
        # print(dim)

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

    attr_binned_by = ncfile.setncattr('Binned_by', 'Dae Sung Choi')
    attr_binned_at = ncfile.setncattr('Binned_at', datetime.datetime.utcnow().replace(tzinfo=datetime.timezone.utc).isoformat())

    if do_write_created_by:
        attr_created_by = ncfile.setncattr('Created_by', created_by)

    if do_write_created_at:
        attr_created_at = ncfile.setncattr('Created_at', created_at)

    ncfile.close()

    print(fn + ' -> ' + outfn + ' DONE.')

