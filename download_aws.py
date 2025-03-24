import os
import json
from datetime import datetime
from concurrent import futures
from concurrent.futures import ThreadPoolExecutor
from functools import partial
#from multiprocessing import Pool, RLock, freeze_support
from random import random
from threading import RLock as TRLock
import uuid
import argparse
import warnings

import pandas as pd
from multiprocess import Pool, RLock, freeze_support
import s3fs
import boto3
import botocore
import botocore.config as botoconfig
from tqdm.auto import tqdm
from tqdm.contrib.concurrent import process_map, thread_map

def aws_download_multithread_worker(save_dir: str, 
                                    bucket: str,
                                    boto3: boto3.resource, 
                                    s3_file: str,
                                    file_name: str, 
                                    progress_pos: int) -> None:

    file_size = int(boto3.Object(bucket, s3_file).content_length)
    pretty_file_name = os.path.basename(s3_file)
    file_path = os.path.join(save_dir, file_name)

    try:
        with tqdm(unit='B', unit_scale=True, unit_divisor=1024, miniters=1, desc=pretty_file_name, total=file_size, leave=None) as progress:
            boto3.Bucket(bucket).download_file(s3_file, file_path, Callback=progress.update)
            progress.close()    
            #progress.display(f"{pretty_file_name}: Finished downloading.", pos=progress_pos)
    except Exception as e:
        print(e)

def aws_download_multithread(save_dir: str, 
                             bucket: str, 
                             keys: list[str, ...], 
                             file_names: list[str, ...]) -> None:
    '''
    Thin wrapper for multithreaded downloading.
    '''
    boto3_session = boto3.Session()
    boto3_client = boto3_session.resource("s3", config=botoconfig.Config(signature_version=botocore.UNSIGNED))

    if not os.path.exists(save_dir):
        os.makedirs(save_dir)

    tqdm.set_lock(TRLock())
    try:
        with ThreadPoolExecutor(initializer=tqdm.set_lock, initargs=(tqdm.get_lock(),)) as executor:
            executor.map(partial(aws_download_multithread_worker, save_dir, bucket, boto3_client), 
                         keys, 
                         file_names, 
                         range(1, len(keys)+1, 1))
    except Exception as e:
        print(e)

def list_aws_keys(bucket: str, 
                  key_pattern: str,
                  glob_match: bool = False):
    '''
    Function to find valid s3 files based on a bucket name and a key pattern.
    Allows for unix-style glob matching.

    Parameters
    ----------
    bucket  : str
        Bucket name to use.
    key_pattern  : str
        Key pattern string to use. If glob_match is False, must be a
        path to a folder containing files (not subfolders), otherwise
        nothing will be downloaded.
        If glob_match is True, uses standard terminology.
    glob_match  : bool, optional [default: False]
        Turns on glob-style matching for key names.

    Returns
    -------
    out  : A list of valid file keys for the given bucket.
    '''
    s3fs_client = s3fs.S3FileSystem(anon=True)

    if not glob_match:
        return [key.replace(f"{bucket}/", "") for key in s3fs_client.ls(f"{bucket}/{key_pattern}")]
    else:
        return [key.replace(f"{bucket}/", "") for key in s3fs_client.glob(f"{bucket}/{key_pattern}")]

if __name__ == "__main__":

    parser = argparse.ArgumentParser(description=f'Download weather data from multiple sources on AWS')
    parser.add_argument('--satellite',
                        help='download specified satellite',
                        choices=['goes16', 
                                 'goes17', 
                                 'goes18', 
                                 'goes16-glm', 
                                 'goes17-glm', 
                                 'goes18-glm', 
                                 'himawari8', 
                                 'himawari9'],
                        type=str,
                        default=None)
    parser.add_argument('--model', 
                        help='download specified model', 
                        choices=['gfs-anl', 
                                 'hrrr-anl', 
                                 'hrrr-anl-ak', 
                                 'nam-anl', 
                                 'nam-anl-ak', 
                                 'nam-anl-hi', 
                                 'nam-anl-pr', 
                                 'nam-anl-fw', 
                                 'gfs-fcst',
                                 'hrrr-fcst', 
                                 'hrrr-fcst-ak', 
                                 'nam-fcst',
                                 'nam-fcst-ak',
                                 'nam-fcst-hi',
                                 'nam-fcst-pr',
                                 'nam-fcst-fw'], 
                        type=str,
                        default=None)
    parser.add_argument('--radar',
                        help='download data from specified radar (or "mrms" for MRMS data)',
                        type=str,
                        metavar=('XXXX'),
                        default=None)
    parser.add_argument('--start-time',  
                        help='download data starting at this time',
                        nargs=4,
                        dest='start_time',
                        type=int,
                        metavar=('yyyy', 'mm', 'dd', 'hh'),
                        default=None)
    parser.add_argument('--end-time',  
                        help='download data up until this time',
                        nargs=4,
                        dest='end_time',
                        type=int,
                        metavar=('yyyy', 'mm', 'dd', 'hh'),
                        default=None)
    parser.add_argument('--around-time',  
                        help='download data around a specified time. an hour delta for both sides can be specified.',
                        nargs=5,
                        dest='around_time',
                        type=int,
                        metavar=('yyyy', 'mm', 'dd', 'hh', 'hr_delta'),
                        default=None)
    parser.add_argument('save_directory',  
                        help='directory to save downloaded data to',
                        type=str)

    args = parser.parse_args()

    print(args)

    if args.save_directory[-1:] != "/":
        save_dir = args.save_directory + "/"
    else:
        save_dir = args.save_directory

    if args.satellite:
        sel_sat = args.satellite
        if "goes" in sel_sat:
            if "glm" in sel_sat:
                product = "glm"
                bucket = f"noaa-{sel_sat[:-4]}"
                dl_freq = '1h'
            else:
                product = "goes"
                bucket = f"noaa-{sel_sat}"
                dl_freq = '1h'
        if "himawari" in sel_sat:
            product = "hw"
            bucket = f"noaa-himawari{sel_sat[-1:]}"
            dl_freq = '10min'
    elif args.model:
        sel_model = product = args.model
        if "gfs" in sel_model:
            bucket = "noaa-gfs-bdp-pds"
            dl_freq = '1h'
        if "hrrr" in sel_model:
            bucket = "noaa-hrrr-bdp-pds"
            dl_freq = '1h'
        if "nam" in sel_model:
            bucket = "noaa-nam-pds"
            dl_freq = '1h'

    elif args.radar:
        if (args.radar == "mrms") or (args.radar == "MRMS"):
            product = "mrms"
            bucket = "noaa-mrms-pds"
            dl_freq = '1h'
        else:
            radar_id = args.radar
            product = "nexrad-l2"
            bucket = 'noaa-nexrad-level2'
            dl_freq = '1h'
    else:
        raise ValueError("Downloading either satellite or model data must be specified.")

    if args.start_time and args.end_time:
        start_ts = pd.Timestamp(args.start_time[0], args.start_time[1], args.start_time[2], args.start_time[3])
        end_ts = pd.Timestamp(args.end_time[0], args.end_time[1], args.end_time[2], args.end_time[3])
        if start_ts > end_ts:
            raise ValueError("Start time must before end time.")
        else:
            dl_range = pd.period_range(start_ts.to_period(dl_freq), end_ts.to_period(dl_freq))
    elif args.around_time and (args.start_time or args.end_time):
        raise ValueError("Please only specify start/end time bounds or a 'time around' a specified time, not both.")
    elif (args.start_time and not args.end_time) or (not args.start_time and args.end_time):
        raise ValueError("When specifying a start/end time, both must be specified.")
    elif args.around_time:
        main_ts = pd.Timestamp(args.around_time[0], args.around_time[1], args.around_time[2], args.around_time[3])
        delta = pd.Timedelta(args.around_time[4], "hr")
        start_ts = main_ts - delta
        end_ts = main_ts + delta

        dl_range = pd.period_range(start_ts.to_period(dl_freq), end_ts.to_period(dl_freq))
    else:
        raise ValueError("Please specify either start/end time bounds or a 'time around' a specified time.")

    # if args.format_for_uems:
    uems_format = False

    keys = []
    file_names = []
    init_ts = None
    stop_pop_flag = False

    if "fcst" in product:
        req_fcst_len = len(dl_range)
        current_fcst_hr = 0
        print(f"Requesting (up to) {req_fcst_len} hours of forecast data, initalized at {start_ts.strftime('%y-%m-%d %H:%M:%S')}.")

    with tqdm(miniters=0, total=len(dl_range), desc=f"Indexing '{bucket}' for requested data...", ascii=" 123456789#") as progress:
        for time in dl_range:
            if not init_ts:
                init_ts = time.to_timestamp()

            this_ts = time.to_timestamp()

            if "goes" in product:
                if "meso" in product:
                    protokey = f"ABI-L2-MCMIPM/{str(this_ts.year).zfill(4)}/{str(this_ts.day_of_year).zfill(3)}/{str(this_ts.hour).zfill(2)}/*"
                else:
                    protokey = f"ABI-L2-MCMIPC/{str(this_ts.year).zfill(4)}/{str(this_ts.day_of_year).zfill(3)}/{str(this_ts.hour).zfill(2)}/*"
                files = list_aws_keys(bucket, protokey, glob_match=True)
                for file in files:
                    keys += [file]
                    file_names += [file[file.rfind('/')+1:]]

            elif product == "glm":
                protokey = f"GLM-L2-LCFA/{str(this_ts.year).zfill(4)}/{str(this_ts.day_of_year).zfill(3)}/{str(this_ts.hour).zfill(2)}/*"
                files = list_aws_keys(bucket, protokey, glob_match=True)
                for file in files:
                    keys += [file]
                    file_names += [file[file.rfind('/')+1:]]

            elif product == "hw":
                protokey = f"AHI-L1b-FLDK/{str(this_ts.year).zfill(4)}/{str(this_ts.month).zfill(2)}/{str(this_ts.day).zfill(2)}/{str(this_ts.hour).zfill(2)}{str(this_ts.minute).zfill(2)}/*"
                files = list_aws_keys(bucket, protokey, glob_match=True)
                for file in files:
                    keys += [file]
                    file_names += [file[file.rfind('/')+1:]]

            elif "hrrr" in product:
                if "ak" in product:
                    file_id = "alaska"
                else:
                    file_id = "conus"

                if "fcst" in product:
                    protokey = f'hrrr.{str(init_ts.year).zfill(4)}{str(init_ts.month).zfill(2)}{str(init_ts.day).zfill(2)}/{file_id}/*wrfprsf*'
                else:
                    protokey = f'hrrr.{str(this_ts.year).zfill(4)}{str(this_ts.month).zfill(2)}{str(this_ts.day).zfill(2)}/{file_id}/*wrfprsf*'

                files = list_aws_keys(bucket, protokey, glob_match=True)
                for file in files:
                    if 'idx' in file:
                        continue

                    #advanced filtering to get hourly resolution
                    file_name = file[file.rfind('/')+1:]

                    file_hrstr = file_name[file_name.find('.'):file_name.find('.')+file_name[file_name.find('.')+1:].find(".")+1]
                    file_hr = int(file_hrstr[2:4])

                    file_fcsthr = int(file_name[file_name.find('wrfprsf')+7:file_name.rfind('.')][:2])

                    if "fcst" in product:
                        if (file_hr == init_ts.hour) and (file_fcsthr == current_fcst_hr):
                            keys += [file]
                            if uems_format:
                                file_names += [file[file.rfind('/')+1:].replace("hrrr.",f"{str(this_ts.year)[2:].zfill(2)}{str(this_ts.month).zfill(2)}{str(this_ts.day).zfill(2)}{str(file_hr).zfill(2)}.hrrr.")]
                            else:
                                file_names += [file[file.rfind('/')+1:].replace("hrrr.",f"hrrr.{str(this_ts.year).zfill(4)}{str(this_ts.month).zfill(2)}{str(this_ts.day).zfill(2)}.")]
                            current_fcst_hr += 1

                            if current_fcst_hr == req_fcst_len:
                                stop_pop_flag = True
                                break
                    else:
                        if (file_hr == this_ts.hour) and (file_fcsthr == 0):
                            keys += [file]
                            if uems_format:
                                file_names += [file[file.rfind('/')+1:].replace("hrrr.",f"{str(this_ts.year)[2:].zfill(2)}{str(this_ts.month).zfill(2)}{str(this_ts.day).zfill(2)}{str(file_hr).zfill(2)}.hrrr.")]
                            else:
                                file_names += [file[file.rfind('/')+1:].replace("hrrr.",f"hrrr.{str(this_ts.year).zfill(4)}{str(this_ts.month).zfill(2)}{str(this_ts.day).zfill(2)}.")]

            elif "gfs" in product:
                if "fcst" in product:
                    protokey = f'gfs.{str(init_ts.year).zfill(4)}{str(init_ts.month).zfill(2)}{str(init_ts.day).zfill(2)}/{str(init_ts.hour).zfill(2)}/atmos/*.pgrb2.0p25.f*'
                else:
                    protokey = f'gfs.{str(this_ts.year).zfill(4)}{str(this_ts.month).zfill(2)}{str(this_ts.day).zfill(2)}/{str(this_ts.hour).zfill(2)}/atmos/*.pgrb2.0p25.anl'
                files = list_aws_keys(bucket, protokey, glob_match=True)
                for file in files:
                    if 'idx' in file:
                        continue

                    #advanced filtering to get hourly resolution
                    file_name = file[file.rfind('/')+1:]

                    file_hrstr = file_name[file_name.find('.'):file_name.find('.')+file_name[file_name.find('.')+1:].find(".")+1]
                    file_hr = int(file_hrstr[2:4])

                    if "fcst" in product:
                        file_fcsthr = int(file_name[file_name.rfind('f')+1:])
                        if (file_hr == init_ts.hour) and (file_fcsthr == current_fcst_hr):
                            keys += [file]
                            file_names += [file[file.rfind('/')+1:].replace("gfs.",f"gfs.{str(this_ts.year).zfill(4)}{str(this_ts.month).zfill(2)}{str(this_ts.day).zfill(2)}.")]

                            #GFS goes from 1-hourly resolution to 3-hourly resolution at >120 forecast hours
                            if current_fcst_hr < 120:
                                current_fcst_hr += 1
                            else:
                                current_fcst_hr += 3

                            if current_fcst_hr >= req_fcst_len:
                                stop_pop_flag = True
                                break

                    else:
                        if (file_hr == this_ts.hour):
                            keys += [file]
                            file_names += [file[file.rfind('/')+1:].replace("gfs.",f"gfs.{str(this_ts.year).zfill(4)}{str(this_ts.month).zfill(2)}{str(this_ts.day).zfill(2)}.")]


            elif "nam" in product:
                if "ak" in product:
                    file_id = "alaskanest"
                elif "hi" in product:
                    file_id = "hawaiinest"
                elif "pr" in product:
                    file_id = "priconest"
                elif "fw" in product:
                    file_id = "firewxnest"
                else:
                    file_id = "conusnest"

                if "fcst" in product:
                    protokey = f'nam.{str(init_ts.year).zfill(4)}{str(init_ts.month).zfill(2)}{str(init_ts.day).zfill(2)}/*{file_id}.hiresf*'
                else:
                    protokey = f'nam.{str(this_ts.year).zfill(4)}{str(this_ts.month).zfill(2)}{str(this_ts.day).zfill(2)}/*{file_id}.hiresf*'

                files = list_aws_keys(bucket, protokey, glob_match=True)
                for file in files:
                    if 'idx' in file:
                        continue

                    #advanced filtering to get hourly resolution
                    file_name = file[file.rfind('/')+1:]

                    file_hrstr = file_name[file_name.find('.'):file_name.find('.')+file_name[file_name.find('.')+1:].find(".")+1]
                    file_hr = int(file_hrstr[2:4])

                    file_fcsthr = int(file_name[file_name.find('hiresf')+6:file_name.rfind('.')-5])

                    if "fcst" in product:
                        if (file_hr == init_ts.hour) and (file_fcsthr == current_fcst_hr):
                            keys += [file]
                            file_names += [file[file.rfind('/')+1:].replace("nam.",f"nam.{str(this_ts.year).zfill(4)}{str(this_ts.month).zfill(2)}{str(this_ts.day).zfill(2)}.")]
                            current_fcst_hr += 1

                            if current_fcst_hr == req_fcst_len:
                                stop_pop_flag = True
                                break
                    else:
                        if (file_hr == this_ts.hour) and (file_fcsthr == 0):
                            keys += [file]
                            file_names += [file[file.rfind('/')+1:].replace("nam.",f"nam.{str(this_ts.year).zfill(4)}{str(this_ts.month).zfill(2)}{str(this_ts.day).zfill(2)}.")]

            elif product == "nexrad-l2":
                protokey = f'{str(this_ts.year).zfill(4)}/{str(this_ts.month).zfill(2)}/{str(this_ts.day).zfill(2)}/{radar_id.upper()}/*'
                files = list_aws_keys(bucket, protokey, glob_match=True)
                for file in files:
                    #advanced filtering to get hourly resolution
                    file_timestr = file[file.find('_')+1:file.find('_')+7] #TODO: make this more dynamic (regex)
                    file_hr = int(file_timestr[0:2])
                    file_min = int(file_timestr[2:4])
                    file_sec = int(file_timestr[4:6])

                    if (file_hr == this_ts.hour) and ('MDM' not in file):
                        keys += [file]
                        file_names += [file[file.rfind('/')+1:]]

            elif product == "mrms":
                protokey = f'CONUS/MergedReflectivityQCComposite_00.50/{str(this_ts.year).zfill(4)}{str(this_ts.month).zfill(2)}{str(this_ts.day).zfill(2)}/*'
                files = list_aws_keys(bucket, protokey, glob_match=True)
                for file in files:
                    #advanced filtering to get hourly resolution
                    file_timestr = file[file.find('-')+1:] #TODO: make this more dynamic (regex)
                    file_hr = int(file_timestr[0:2])
                    file_min = int(file_timestr[2:4])
                    file_sec = int(file_timestr[4:6])

                    if (file_hr == this_ts.hour):
                        keys += [file]
                        file_names += [file[file.rfind('/')+1:]]



            if stop_pop_flag:
                progress.update(len(dl_range))
                break
            progress.update()

    if (not keys) or (not file_names):
        raise FileNotFoundError("No files found on AWS. Try reviewing input parameters, or change start time for forecast data.")
    else:
        print(f'{len(keys)} files found. Starting download...')

    aws_download_multithread(save_dir, bucket, keys, file_names)

