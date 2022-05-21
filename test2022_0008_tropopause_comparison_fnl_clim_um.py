import numpy as np
from netCDF4 import Dataset
import pandas as pd

import datetime

import matplotlib.pyplot as plt
import os

dr = pd.date_range(datetime.datetime(2021, 3, 1, 00), 
        datetime.datetime(2021, 3, 31, 18), freq='1d')


for idr in dr:
    yyyy = str(idr.year).zfill(4)
    mm = idr.strftime('%m')
    dd = idr.strftime('%d')
    hh = idr.strftime('%H')
    mi = idr.strftime('%M')

    #--------------------------------------------------------------------------
    # fn = f'/data/MODEL/FNL/2021/fnl_{yyyy}{mm}{dd}_{hh}_{mi}.grib2.nc'
    # dataset = Dataset(fn)
    # lat_0 = dataset['lat_0'][:] # latitude for 1d
    # lon_0 = dataset['lon_0'][:] # latitude for 1d

    # # make lon lat to 2d 
    # fnl_lon2d, fnl_lat2d = np.meshgrid(lat_0, lon_0)

    # # get tropopause from nc
    # # PRES_P0_L1_GLL0 = dataset['PRES_P0_L1_GLL0'][:] # Ground or water surface [Pa]
    # PRES_P0_L7_GLL0 = dataset['PRES_P0_L7_GLL0'][:] # Tropopause Pressure
    # # TMP_P0_L1_GLL0 = dataset['TMP_P0_L1_GLL0'][:] # Surface Temperature [K]
    # # TMP_P0_L100_GLL0 = dataset['TMP_P0_L100_GLL0'][:] # Isobaric surface Pressure 34 level

    #--------------------------------------------------------------------------
    # # read Surface Temperature to check North-South order
    # fnl_st_clim_fn = f'/OPER/SYSTEMS/ALGRTH/GEMS/L2/data/o3p/ATMOS/fnl13.75LST/fnlst/fnlstavg{mm}.dat'
    # f = open(fnl_st_clim_fn, 'r')
    # fnl_sfct = np.zeros([180, 360])
    # for jy in range(180):
        # line = f.readline()
        # for ix in range(360):
            # fnl_sfct[jy, ix] = float(line[ix*3:(ix+1)*3])

    # fnl_lon_dat = np.arange(360, dtype=np.float32) - 180.
    # fnl_lat_dat = np.arange(180, dtype=np.float32) - 90. # converted to Ascending order in PREPARE_FNL_GEUN

    # fnl_lon_dat_2d, fnl_lat_dat_2d = np.meshgrid(fnl_lon_dat, fnl_lat_dat)

    # plt.contourf(fnl_lon_dat_2d, fnl_lat_dat_2d, fnl_sfct, levels=np.arange(250, 350, 1))
    # plt.show()

    #--------------------------------------------------------------------------
    fnl_tp_clim_fn = f'/OPER/SYSTEMS/ALGRTH/GEMS/L2/data/o3p/ATMOS/fnl13.75LST/fnltp/fnltpavg{mm}.dat'
    f = open(fnl_tp_clim_fn, 'r')
    fnl_tpres = np.zeros([180, 360])
    for jy in range(180):
        line = f.readline()
        for ix in range(360):
            fnl_tpres[jy, ix] = float(line[ix*3:(ix+1)*3])

    fnl_lon_dat = np.arange(360, dtype=np.float32) - 180. + 0.5
    fnl_lat_dat = np.arange(180, dtype=np.float32) - 90. + 0.5 # converted to Ascending order in PREPARE_FNL_GEUN

    fnl_lon_dat_2d, fnl_lat_dat_2d = np.meshgrid(fnl_lon_dat, fnl_lat_dat)

    # plt.contourf(fnl_lon_dat_2d, fnl_lat_dat_2d, fnl_tpres, levels=np.arange(50, 350, 3))
    # plt.show()
    # quit()

    #--------------------------------------------------------------------------
    # read um dat file
    um_tpres = np.zeros([769, 1024])
    fnl_tpres_in_umgrid = np.zeros([769, 1024])

    # # read Surface Temperature to check North-South order
    # umdatfn = f'/OPER/SYSTEMS/ALGRTH/GEMS/L2/data/o3p/ATMOS/umatmos/umst/umst_{yyyy}{mm}{dd}.dat'
    # um_sfct = np.zeros([769, 1024])
    # f = open(umdatfn, 'r')
    # for jy in range(769):
        # line = f.readline()
        # for ix in range(1024):
            # um_sfct[jy, ix] = float(line[ix*3:(ix+1)*3])

    # um_lon_dat = np.arange(1024, dtype=np.float32)/1024.*360. - 180
    # um_lat_dat = np.arange(769, dtype=np.float32)/768.*180. - 90

    # um_lon_dat_2d, um_lat_dat_2d = np.meshgrid(um_lon_dat, um_lat_dat)

    # plt.contourf(um_lon_dat_2d, um_lat_dat_2d, um_sfct, levels=np.arange(200, 310, 1))
    # plt.show()
    # quit()


    # calculated from a day before 18 utc 6 hour prediction temperature profile
    umdatfn = f'/OPER/SYSTEMS/ALGRTH/GEMS/L2/data/o3p/ATMOS/umatmos/umtp/umtp_{yyyy}{mm}{dd}.dat'
    f = open(umdatfn, 'r')
    for jy in range(769):
        line = f.readline()
        for ix in range(1024):
            um_tpres[jy, ix] = float(line[ix*3:(ix+1)*3])

    temp = np.zeros([769, 1024])

    # fnl_tpres from dat file that saved in -180~180 range to 0~360
    # -90~90

    # make lon lat for um
    um_lon_dat = np.arange(1024, dtype=np.float32)/1024.*360. - 180
    um_lat_dat = np.arange(769, dtype=np.float32)/768.*180. - 90

    um_lon_dat_2d, um_lat_dat_2d = np.meshgrid(um_lon_dat, um_lat_dat)

    fnl_col_idx_for_um = np.zeros([769, 1024], dtype=np.int32)

    print(fnl_lon_dat_2d[0, 0], fnl_lon_dat_2d[0, -1], fnl_lon_dat_2d[-1, 0], fnl_lon_dat_2d[-1, -1])
    print(fnl_lat_dat_2d[0, 0], fnl_lat_dat_2d[0, -1], fnl_lat_dat_2d[-1, 0], fnl_lat_dat_2d[-1, -1])
    print(um_lon_dat_2d[0, 0], um_lon_dat_2d[0, -1], um_lon_dat_2d[-1, 0], um_lon_dat_2d[-1, -1])
    print(um_lat_dat_2d[0, 0], um_lat_dat_2d[0, -1], um_lat_dat_2d[-1, 0], um_lat_dat_2d[-1, -1])

    # fnl is in the um boundary

    # collocate FNL to UM
    cllc_idx_fn = '/data/private/soodal/um_tp_comparison/fnl_dat_collocation_idx_for_um_dat.npy'
    if not os.path.exists(cllc_idx_fn):
        for iy in np.arange(0, 769):
            for ix in np.arange(0, 1024):
                distance = (fnl_lon_dat_2d[:] - um_lon_dat_2d[iy, ix])**2 + (fnl_lat_dat_2d[:] - um_lat_dat_2d[iy, ix])**2
                idx = np.argmin(distance)
                fnl_col_idx_for_um[iy, ix] = idx

        np.save(cllc_idx_fn, fnl_col_idx_for_um)
    else:
        fnl_col_idx_for_um = np.load(cllc_idx_fn)

    print('-'*80)
    print(fnl_col_idx_for_um.shape)
    print('-'*80)
    fnl_col_yidx_for_um, fnl_col_xidx_for_um = np.unravel_index(fnl_col_idx_for_um, (180, 360))


    print('-'*80)
    print('fnl_col_yidx_for_um')
    print(fnl_col_yidx_for_um)
    print('-'*80)
    print('fnl_col_xidx_for_um')
    print(fnl_col_xidx_for_um)
    print('-'*80)

    for iy in np.arange(0, 769):
        for ix in np.arange(0, 1024):
            fnl_tpres_in_umgrid[iy, ix] = fnl_tpres[fnl_col_yidx_for_um[iy, ix],
                                                    fnl_col_xidx_for_um[iy, ix]]

    # scatter plot for every point
    plt.scatter(um_lat_dat_2d[:], um_tpres[:], s=1, alpha=0.1, label='UM')
    plt.scatter(um_lat_dat_2d[:], fnl_tpres_in_umgrid[:], s=1, 
            alpha=0.1, label='FNL monthly')
    plt.ylim(500, 50)
    plt.grid()
    plt.legend()
    plt.title('Tropopause Comparison UM(WMO Thermal TP) with FNL')
    plt.xlabel('Latitude')
    plt.ylabel('Pressure')
    plt.savefig(f'/home/soodal/works/GEMS_O3P_analysis/plot/um_tp_comparison/{yyyy}{mm}{dd}_01_all_scatterplot.png')
    plt.close()

    plt.scatter(um_lat_dat_2d[:], um_tpres[:] - fnl_tpres_in_umgrid[:], s=1, 
            alpha=0.1, label='UM - FNL')
    plt.ylim(-150, 250)
    plt.grid()
    plt.legend()
    plt.title('Tropopause Pressure Difference : UM(WMO Thermal TP) - FNL')
    plt.xlabel('Latitude')
    plt.ylabel(r'$\nabla$Pressure')
    plt.savefig(f'/home/soodal/works/GEMS_O3P_analysis/plot/um_tp_comparison/{yyyy}{mm}{dd}_02_diff_scatterplot.png')
    plt.close()


    # freq = np.zeros([30, 18])

    # delta_tp = PRES




