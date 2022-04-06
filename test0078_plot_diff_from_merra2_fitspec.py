import os
import numpy as np
from netCDF4 import Dataset

from matplotlib import pyplot as plt

# import numpy as np
from scipy.optimize import curve_fit
# import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import datetime

import pandas as pd
def func(X, asza, bsza, csza, avza, bvza, cvza):
    xsza, xvza = X
    xsza = np.sin(np.deg2rad(xsza))
    xvza = np.sin(np.deg2rad(xvza))

    return asza * xsza**2 + bsza*xsza + csza + avza * xvza**2 + bvza*xvza + cvza

def file_search(dirname):
    filenames = os.listdir(dirname)
    full_filenames = []
    for filename in filenames:
        full_filename = os.path.join(dirname, filename)
        full_filenames.append(full_filename)
    return full_filenames


datetimelist = pd.date_range(start='2020-08-01T0345', end='2021-06-30T0745', freq='h')
for a_datetime in datetimelist:
    yyyymmdd_hhmi = a_datetime.strftime('%Y%m%d_%H%M')
    yesterday_datetime = a_datetime - datetime.timedelta(days=1)
    yyyymmdd_hhmi_yesterday = yesterday_datetime.strftime('%Y%m%d_%H%M')
    irr_file = '/data/private/soodal/ln/GEMS/IRR/GK2_GEMS_IRR_' + yyyymmdd_hhmi_yesterday + '_4x4.nc'
    l1c_file = '/data/private/soodal/ln/GEMS/L1C/GK2_GEMS_L1C_' + yyyymmdd_hhmi + '_4x4.nc'
    diff_file = '/data/private/soodal/softcal_test/residual/gems_merra2_' + yyyymmdd_hhmi + '_raw.nc'
    irr_file_exist = os.path.exists(irr_file)
    l1c_file_exist = os.path.exists(l1c_file)
    diff_file_exist = os.path.exists(diff_file)

    if l1c_file_exist and irr_file_exist and diff_file_exist:
        # read diff_file
        rootgrp = Dataset(diff_file, "r", format="NETCDF4")
        # sza = rootgrp['solarzenithangle'][:]
        # vza = rootgrp['viewingzenithangle'][:]
        # wvl = rootgrp['wavelengths'][:]
        diff = rootgrp['fitspec_difference_from_merra2'][:]
        # ecf = rootgrp['effectivecloudfraction'][:]
        # print(ecf.shape)
        rootgrp.close()

        # read l1c_file
        rootgrp = Dataset(l1c_file, "r", format="NETCDF4")
        # sza = rootgrp['solarzenithangle'][:]
        # vza = rootgrp['viewingzenithangle'][:]
        # wvl = rootgrp['wavelengths'][:]
        diff = rootgrp['image_pixel_values'][:]
        rootgrp.close()

        # read irr_file
        rootgrp = Dataset(irr_file, "r", format="NETCDF4")
        # sza = rootgrp['solarzenithangle'][:]
        # vza = rootgrp['viewingzenithangle'][:]
        # wvl = rootgrp['wavelengths'][:]
        diff = rootgrp['image_pixel_values'][:]
        rootgrp.close()
        


filelist = file_search('/data/private/soodal/softcal_test/residual/')
print(filelist)

sza_list = []
vza_list = []
diff_list = []
ecf_list = []
wvl_list = []

for file in filelist:
    rootgrp = Dataset(file, "r", format="NETCDF4")
    sza = rootgrp['solarzenithangle'][:]
    vza = rootgrp['viewingzenithangle'][:]
    wvl = rootgrp['wavelengths'][:]
    diff = rootgrp['fitspec_difference_from_merra2'][:]
    ecf = rootgrp['effectivecloudfraction'][:]
    print(ecf.shape)

    diff_list.append(diff)
    sza_list.append(sza)
    vza_list.append(vza)
    ecf_list.append(ecf)
    wvl_list.append(wvl)
    # print(wvl.shape)
    # quit()


iy = 59
# for iy in range(512):
sy = str(iy).zfill(3)
iw = 10
    # for iw in range(198):

sw = str(iw).zfill(3)



p0 = [ 1.96326535e+00, -1.97739394e+00, -3.55075300e+03,  1.25250876e+01, 
         -2.03555294e+01,  3.55964733e+03]

# for ispatial in range(sza_list[ifile].shape[0]):
for ispatial in range(157, 512):
    # for iwav in range(diff_list[ifile].shape[0]):
    sza_wavs = []
    vza_wavs = []
    diff_wavs = []
    wvl_wavs = []
    popt_wavs = []
    for iwav in range(198):
        sza_files = np.array([])
        vza_files = np.array([])
        diff_files = np.array([])
        wvl_files = np.array([])
        popt_files = np.array([])
        for ifile in range(len(sza_list)):
            sf = str(ifile).zfill(2)
            # print(sza_list[ifile].shape)
            # print(diff_list[ifile].shape)

            sza_vza_idx = np.logical_and(sza_list[ifile][ispatial, :] > -1e29, vza_list[ifile][10, :] > -1e29)
            ecf_idx = np.logical_and(sza_vza_idx, ecf_list[ifile][ispatial, :] < 0.2)
            x_finiteidx = np.logical_and(ecf_idx, diff_list[ifile][iwav, ispatial, :] > -1e20)
            # print(sza_list[ifile][10, :])
            # print(sza_list[ifile][10, :] > -1e29)
            # print(ecf_list[ifile][10, :])
            # print(ecf_list[ifile][10, :] < 0.2)
            # print(x_finiteidx)
            # print(x_finiteidx.shape)

            # type(sza_list[ifile])
            # print(diff_list[ifile].shape)
            # print(diff_list[ifile][:, 10, x_finiteidx])
            sza_files = np.concatenate((sza_files, sza_list[ifile][ispatial, x_finiteidx].ravel()), axis=0)
            vza_files = np.concatenate((vza_files, vza_list[ifile][ispatial, x_finiteidx].ravel()), axis=0)
            diff_files = np.concatenate((diff_files, diff_list[ifile][iwav, ispatial, x_finiteidx].ravel()), axis=0)
            wvl_files = np.concatenate((wvl_files, wvl_list[ifile][iwav, ispatial].ravel()), axis=0)

# =============================================================================

        if len(sza_files) > 1:
            # print(p0)
            idx = np.logical_and(np.isfinite(sza_files), np.isfinite(vza_files))
            idx = np.logical_and(idx, np.isfinite(diff_files))
            popt, pcov = curve_fit(func, (sza_files[idx], vza_files[idx]), diff_files[idx], p0)
# print(popt)
# print(pcov)
            # print(popt.shape)
            popt_files = np.concatenate((popt_files, popt), axis=0)
            # print('popt_files:', popt_files)
            # print('popt_files.shape:', popt_files.shape)
            # if iwav == 10:
                # quit()


            model_diff = (np.sin(np.deg2rad(sza_files))**2 * popt[0] + np.sin(np.deg2rad(sza_files)) * popt[1] + popt[2] +
                np.sin(np.deg2rad(vza_files)) **2 *  popt[3] + np.sin(np.deg2rad(vza_files)) * popt[4] + popt[5])


# Plot the 3D figure of the fitted function and the residuals.
            fig = plt.figure()
# # ax = fig.gca(projection='3d') !< for trisurf
            ax = fig.add_subplot(projection='3d')
# # ax.plot_trisurf(sza_files, vza_files, diff_files, cmap='plasma')
            ax.scatter(sza_files, vza_files, diff_files, cmap='plasma')
            ax.set_xlabel('SZA')
            ax.set_ylabel('VZA')
            ax.set_zlabel('Difference')

            ax.plot_trisurf(sza_files, vza_files, model_diff, cmap='plasma')

# # ax.scatter(sza_files, vza_files, diff_files, cmap='plasma')
# # ax.set_zlim(-np.max(abs(diff_files)),np.max(abs(diff_files))+2)
            # plt.show()
            savefn = ('./plot/diff_modelling/wav' + str(iwav).zfill(3) 
                    + '_spatial' + str(ispatial).zfill(3) + '.png')
            print(savefn)
            plt.savefig(savefn)
            plt.close()


            if iwav == 20:
                fig2 = plt.figure()
                ax = fig.add_subplot()
                for ix in range(len(x_finiteidx)):
                    if x_finiteidx[ix]:
                        ax.plot(wvl_list[ifile][:, ispatial].ravel(), 
                                diff_list[ifile][:, ispatial, ix].ravel())
                savefn = ('./plot/diff_modelling/image' + str(ix).zfill(3) 
                        + '_spatial' + str(ispatial).zfill(3) + '.png')


        if len(sza_files) > 1 :
            # np.round(np.random.random(1)*100)
            sza_wavs.append(sza_files)
            vza_wavs.append(vza_files)
            diff_wavs.append(diff_files)
            wvl_wavs.append(wvl_files)

quit()

fig = plt.figure()
# ax = fig.gca(projection='3d') !< for trisurf
ax = fig.add_subplot(projection='3d')
# ax.plot_trisurf(sza_files, vza_files, diff_files, cmap='plasma')
ax.scatter(sza_files, vza_files, diff_files, cmap='plasma')
# ax.set_zlim(-np.max(abs(diff_files)),np.max(abs(diff_files))+2)
ax.set_xlabel('SZA')
ax.set_ylabel('VZA')
ax.set_zlabel('Difference')
plt.show()
quit()

