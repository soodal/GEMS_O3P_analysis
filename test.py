import irradpy
import os
import pandas as pd

from requests import get  # to make GET request
import datetime


# socket = irradpy.downloader.download.SocketManager()


# build url
# url = os.path.join(socket.build_remote_url(merra2_collection, date),
                # socket.build_remote_filename(merra2_collection, date, params))

class MerraDownloader():
    def __init__(self):
        self.uid = 'daeseongchoe'
        self.pw = 'tfa2Adds'
        self.socketmanager = irradpy.downloader.download.SocketManager()
        self.dm = irradpy.downloader.process.DownloadManager()
        self.dm.set_username_and_password(self.uid, self.pw)
        self.dm.download_path = './MERRA2_data/'


    # def irradpy_merra_o3_download(self, 
            # initial_year, initial_month, initial_day,
            # final_year, final_month, final_day,
            # output_directory,
            # # auth=auth,
            # # params=parameter,
            # thread_num=thread_num,
            # merge=merge_timelapse,
            # ):

        # for i, collection_name in enumerate(collection_names):
            # if not merra2_var_dicts:
                # merra2_var_dict = var_list[collection_name]
            # else:
                # merra2_var_dict = merra2_var_dicts[i]
        # if isinstance(merra2_var_dict['var_name'], list):
            # requested_params = merra2_var_dict['var_name']
        # else:
            # requested_params = [merra2_var_dict['var_name']]
        
        # if merra2_var_dict["collection"].startswith("const"):
            # requested_time = '[0:0]'
        # else:
            # requested_time = '[0:23]'
       
        # parameter = self.socketmanager.generate_url_params(
                # requested_params, requested_time, requested_lat, requested_lon)

        # self.socketmanager.subdaily_universal_download(
            # merra2_var_dict,
            # initial_year,
            # final_year,
            # initial_month=initial_month,
            # final_month=final_month,
            # initial_day=initial_day,final_day=final_day,
            # output_directory=output_directory,
            # auth={'uid':self.uid, 'password':self.pw}, 
            # # params=parameter,
            # thread_num=thread_num,
            # # merge=merge_timelapse,
            # )



    def merra_o3_download(self, year, month, day, download_path):
        if download_path is None:
            download_path = './MERRA2_data/'


        fn = 'MERRA2_400.inst3_3d_chm_Nv.' + f'{year:04}{month:02}{day:02}' + '.nc4'
        url = ('https://goldsmr5.gesdisc.eosdis.nasa.gov/data/MERRA2/'
                + 'M2I3NVCHM.5.12.4/' + f'{year:04}/{month:02}/' + fn)
        print(url)
        self.dm.download_url = url
        self.dm.download_path = download_path
        self.dm.start_download()


if __name__ == '__main__':
    md = MerraDownloader()
    # md.irradpy_merra_o3_download(2020, 8, 1, 2020, 9, 1, './MERRA2_data/', 
            # thread_num=20)
    
    initial_datetime = datetime.datetime(2020, 8, 1)
    final_datetime = datetime.datetime(2021, 6, 30)

    datetime_list = pd.date_range(
            start=initial_datetime, end=final_datetime).to_pydatetime().tolist()

    for idate in range(len(datetime_list)):
        the_datetime = datetime_list[idate]
        year = the_datetime.year
        month = the_datetime.month
        day = the_datetime.day

        print(f'{year}-{month}-{day}')

        md.merra_o3_download(year, month, day, f'/data/MERRA2/{year:04}/{month:02}/')

