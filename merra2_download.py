import irradpy
import os
import pandas as pd

from requests import get  # to make GET request
import datetime


class MerraDownloader():
    def __init__(self):
        self.uid = 'daeseongchoe'
        self.pw = 'tfa2Adds'
        self.socketmanager = irradpy.downloader.download.SocketManager()
        self.dm = irradpy.downloader.process.DownloadManager()
        self.dm.set_username_and_password(self.uid, self.pw)
        self.dm.download_path = './MERRA2_data/'

    def merra_o3_download(self, year, month, day, download_path):
        if download_path is None:
            download_path = './MERRA2_data/'

        version = 400
        if year == 2020 and month == 9:
            version=401
        fn = 'MERRA2_'+f'{version:3}' + '.inst3_3d_chm_Nv.' + \
            f'{year:04}{month:02}{day:02}' + '.nc4'
        url = ('https://goldsmr5.gesdisc.eosdis.nasa.gov/data/MERRA2/'
               + 'M2I3NVCHM.5.12.4/' + f'{year:04}/{month:02}/' + fn)
        print(url)
        self.dm.download_url = url
        self.dm.download_path = download_path
        self.dm.start_download()
        fn = 'MERRA2_'+f'{version:3}' + '.inst3_3d_chm_Nv.' + \
            f'{year:04}{month:02}{day:02}' + '.nc4.xml'
        url = ('https://goldsmr5.gesdisc.eosdis.nasa.gov/data/MERRA2/'
               + 'M2I3NVCHM.5.12.4/' + f'{year:04}/{month:02}/' + fn)
        print(url)
        self.dm.download_url = url
        self.dm.download_path = download_path
        self.dm.start_download()

    def merra_asm_download(self, year, month, day, download_path):
        if download_path is None:
            download_path = './MERRA2_data/'

        version = 400
        if year == 2020 and month == 9:
            version=401
        fn = 'MERRA2_'+f'{version:3}' + '.tavg3_3d_asm_Nv.' + \
            f'{year:04}{month:02}{day:02}' + '.nc4'
        url = ('https://goldsmr5.gesdisc.eosdis.nasa.gov/data/MERRA2/'
               + 'M2T3NVASM.5.12.4/' + f'{year:04}/{month:02}/' + fn)
        print(url)
        self.dm.download_url = url
        self.dm.download_path = download_path
        self.dm.start_download()
        fn = 'MERRA2_'+f'{version:3}' + '.tavg3_3d_asm_Nv.' + \
            f'{year:04}{month:02}{day:02}' + '.nc4.xml'
        url = ('https://goldsmr5.gesdisc.eosdis.nasa.gov/data/MERRA2/'
               + 'M2I3NVCHM.5.12.4/' + f'{year:04}/{month:02}/' + fn)
        print(url)
        self.dm.download_url = url
        self.dm.download_path = download_path
        self.dm.start_download()

if __name__ == '__main__':
    md = MerraDownloader()

    initial_datetime = datetime.datetime(2021, 6, 1)
    final_datetime = datetime.datetime(2021, 6, 30)

    datetime_list = pd.date_range(
        start=initial_datetime, end=final_datetime).to_pydatetime().tolist()

    for idate in range(len(datetime_list)):
        the_datetime = datetime_list[idate]
        year = the_datetime.year
        month = the_datetime.month
        day = the_datetime.day

        print(f'download start MERRA2 M2I3NVCHM {year:04}-{month:02}-{day:02}')

        # md.merra_o3_download(
            # year, month, day, f'/data/MODEL/MERRA2/{year:04}/{month:02}/')
        md.merra_asm_download(
            year, month, day, f'/data/MODEL/MERRA2/{year:04}/{month:02}/')
