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

def create_alt_az_table(listing_lat: float,
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
    moon = eph['Moon']

    obs_pos = skapi.wgs84.latlon(wgs_lat, 
                                 wgs_lon, 
                                 elevation_m=kwargs.pop('elevation_m', 0))

    location = eph['Earth'] + obs_pos

    if kwargs.get('local_timezone', False):
        local_ts = pytz.timezone(kwargs.pop('local_timezone'))
    else:
        local_ts = None

    records = []

    period = pd.period_range(pd.Timestamp(start_date_utc).to_period('1min'), 
                             pd.Timestamp(end_date_utc).to_period('1min'))

    for time in period:
        py_time = time.to_timestamp().to_pydatetime()
        utc_time = py_time.replace(tzinfo=pytz.utc)
        tz_time = local_ts.normalize(utc_time)
        obs_time = ts.from_datetime(utc_time)

        sunPos = location.at(obs_time).observe(sun)
        moonPos = location.at(obs_time).observe(moon)
        sunAlt, sunAz, sunDist = sunPos.apparent().altaz()
        moonAlt, moonAz, moonDist = moonPos.apparent().altaz()

        obs = (py_time.strftime("%Y-%m-%d %H:%M:%S"), 
               tz_time.strftime("%Y-%m-%d %H:%M:%S"), 
               sunAlt.degrees, 
               sunAz.degrees,
               sunDist.km, 
               moonAlt.degrees, 
               moonAz.degrees,
               moonDist.km)
        records += [obs]

    df = pd.DataFrame(records, columns=["Time (UTC)", f"Time ({local_tz})", "Sun Altitude (°)", "Sun Azimuth (°)", "Sun Distance (km)", "Moon Altitude (°)", "Moon Azimuth (°)", "Moon Distance (km)"])
    df.to_csv(os.path.join(save_directory, "AltAzTable.csv"), index=False)

if __name__ == "__main__":

    parser = argparse.ArgumentParser(description=f'Create a listing data file from model data')
    parser.add_argument('--local-timezone',  
                        help='download data starting at this time',
                        dest='timezone',
                        type=str,
                        metavar=('tz'),
                        default=None)
    parser.add_argument('--location', 
                        help='specify a point to take the listing data from', 
                        nargs=2, 
                        type=float,
                        dest='loc',
                        metavar=('lat', 'lon'),
                        default=None)
    parser.add_argument('--elevation', 
                        help='specify an elevation in meters at the observation point', 
                        nargs=1, 
                        type=float,
                        dest='elev_m',
                        metavar=('elev_m'),
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
    if args.elev_m:
        user_settings.update({'elevation_m': args.elev_m})

    if args.timezone:
        user_settings.update({'local_timezone': args.timezone})

    print(f"Creating alt/az table listing for {listing_lat}, {listing_lon}, between {start_date.strftime('%Y-%m-%d %H:%M:%S %Z')} and {end_date.strftime('%Y-%m-%d %H:%M:%S %Z')}...")

    create_alt_az_table(listing_lat,
                        listing_lon,
                        start_date, 
                        end_date,
                        save_dir,
                        **user_settings) 

    print("Done!")






