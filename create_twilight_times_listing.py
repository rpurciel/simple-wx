## The file at the following link lists units that can be used in this context:
## https://github.com/hgrecco/pint/blob/master/pint/default_en.txt
## A dimensionality check is automatically performed on any specified units

## Default units for pressure values
## (dimensionality: '[mass] / [length] / [time] ** 2')
DEFAULT_PRESSURE_UNITS = 'hectopascal'

## Default units for temperature/dewpoint values
## (dimensionality: '[temperature]')
DEFAULT_TEMP_UNITS = 'degC'

## Default units for wind speed/wind gust values
## (dimensionality: '[length] / [time]')
DEFAULT_WIND_UNITS = 'kt'


## Set to True to use MISSING_FLAG for when a wind gust speed for WRF
## cannot be calculated due to a low PBL. (default: True)
USE_MISSING_FOR_WRF_GUSTS = True

## Flag to use for any missing values. Can be a number, 'None', or
## np.nan for NaN. (default: -999)
MISSING_FLAG = -999









'''


            STATIC CODE BELOW HERE


'''

import os
import glob
import datetime
import argparse
import csv
import warnings
import pytz

import numpy as np
import skyfield.api as skapi
from skyfield import almanac

from __internal_funcs import get_nearest_feature

np.seterr(divide = 'ignore', invalid='ignore')

def create_twilight_times_listing(listing_lat: float,
                                  listng_lon: float,
                                  start_date_utc: datetime.datetime,
                                  end_date_utc: datetime.datetime,
                                  save_directory: str,
                                  **kwargs) -> None:

    wgs_lat = np.abs(listing_lat) * (skapi.N if listing_lat >=0 else skapi.S) 
    wgs_lon = np.abs(listing_lon) * (skapi.E if listing_lon >=0 else skapi.W)

    ts = skapi.load.timescale()
    eph = skapi.load('de421.bsp')
    sun = eph['Sun']
    obs_pos = skapi.wgs84.latlon(wgs_lat, wgs_lon)
    observer = eph['Earth'] + obs_pos

    alm_times = almanac.dark_twilight_day(eph, obs_pos) 

    start_time = ts.from_datetime(start_date_utc)
    end_time = ts.from_datetime(end_date_utc)
    times, events = almanac.find_discrete(start_time, end_time, alm_times)

    save_path = os.path.join(save_directory, 'TwilightTimes_Listing.txt')

    header = f"Twilight times for {listing_lat}, {listing_lon} ({get_nearest_feature(listing_lat, listing_lon, show_distance=False)})\nBetween {start_date.strftime('%Y-%m-%d %H:%M:%S %Z')} and {end_date.strftime('%Y-%m-%d %H:%M:%S %Z')}\n\n\n"

    if kwargs.get('local_timezone', False):
        local_ts = pytz.timezone(kwargs.pop('local_timezone'))
    else:
        local_ts = None

    previous_e = alm_times(start_time).item()
    with open(save_path, 'w') as listing_file:
        listing_file.write(header)
        for t, e in zip(times, events):
            pytime = t.utc_datetime()
            local_timestr = f" [{pytime.astimezone(local_ts).strftime('%Y-%m-%d %H:%M:%S %Z')}] " if local_ts is not None else ' '
            if previous_e < e:
                twilight = almanac.TWILIGHTS[e]
                line = f"{pytime.strftime('%Y-%m-%d %H:%M:%S %Z')}{local_timestr}- {'Sunrise' if twilight == 'Day' else f'{twilight} starts'}\n"
                listing_file.write(line)
            else:
                twilight = almanac.TWILIGHTS[previous_e]
                line = f"{pytime.strftime('%Y-%m-%d %H:%M:%S %Z')}{local_timestr}- {'Sunset' if twilight == 'Day' else f'{twilight} ends'}\n"
                listing_file.write(line)
                if e == 0:
                    listing_file.write("---\n")
            previous_e = e

    listing_file.close()

if __name__ == "__main__":

    parser = argparse.ArgumentParser(description=f'Create a listing data file from model data')
    parser.add_argument('--local-timezone',  
                        help='download data starting at this time',
                        dest='timezone',
                        type=str,
                        metavar=('tz'),
                        choices=pytz.all_timezones,
                        default=None)
    parser.add_argument('--location', 
                        help='specify a point to take the listing data from', 
                        nargs=2, 
                        type=float,
                        dest='loc',
                        metavar=('lat', 'lon'),
                        default=None)
    parser.add_argument('--start-date',  
                        help='download data starting at this time',
                        nargs=3,
                        type=int,
                        dest='start_date',
                        metavar=('yyyy', 'mm', 'dd'),
                        default=None)
    parser.add_argument('--end-date',  
                        help='download data starting at this time',
                        nargs=3,
                        type=int,
                        dest='end_date',
                        metavar=('yyyy', 'mm', 'dd'),
                        default=None)
    parser.add_argument('save_directory',  
                        help='directory to save listing to',
                        type=str)

    args = parser.parse_args()
    print(args)

    ## Quality Checks/Sanitizing of input arguments

    if args.loc:
        listing_lat = args.loc[0]
        listing_lon = args.loc[1]
        if (listing_lat > 90) or (listing_lat < -90):
            raise ValueError("Invalid latitude for listing point, must +90 =< x =< -90")
        if (listing_lon > 180) or (listing_lon < -180):
            raise ValueError("Invalid longitude for listing point, must +180 =< x =< -180")

    start_date = datetime.datetime(args.start_date[0],
                                  args.start_date[1],
                                  args.start_date[2],
                                  0,
                                  0,
                                  0,
                                  tzinfo=pytz.utc)

    end_date = datetime.datetime(args.end_date[0],
                          args.end_date[1],
                          args.end_date[2],
                          0,
                          0,
                          0,
                          tzinfo=pytz.utc)

    if args.save_directory[-1:] != "/":
        save_dir = args.save_directory + "/"
    else:
        save_dir = args.save_directory

    print(f"Saving to '{save_dir}'")

    user_settings = dict()
    if args.timezone:
        user_settings.update({'local_timezone': args.timezone})

    print(f"Creating twilight times listing for {listing_lat}, {listing_lon}, between {start_date.strftime('%Y-%m-%d %H:%M:%S %Z')} and {end_date.strftime('%Y-%m-%d %H:%M:%S %Z')}...")

    create_twilight_times_listing(listing_lat,
                                  listing_lon,
                                  start_date,
                                  end_date,
                                  save_dir,
                                  **user_settings) 

    print("Done!")






