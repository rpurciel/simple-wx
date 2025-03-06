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
from datetime import datetime
import argparse
import csv
import warnings

import xarray as xr
import numpy as np
import pandas as pd
import metpy
from metpy.units import units
import metpy.calc as mpcalc
from tqdm.auto import tqdm
import pint
import wrf
from wrf import getvar
from netCDF4 import Dataset

from __internal_funcs import plot_towns, generate_verification_plot

np.seterr(divide = 'ignore', invalid='ignore')

NETCDF_MODELS = ['wrf']

GRIB2_MODELS = ['hrrr', 
                'gfs', 
                'era5']

def generate_header(model: str, 
                    ureg: pint.UnitRegistry, 
                    tempunit_spec: str, 
                    windunit_spec: str, 
                    pressunit_spec: str) -> tuple:

    """
    Dynamically generate the header for the output CSV, using input units.

    Inputs
        - model, str:
                Model identifier string.
        - ureg, pint.UnitRegistry:
                A UnitRegistry object to be used to create unit identifiers.
        - tempunit_spec: str,
                String specifying what unit to be used for temperature. Must 
                be found in ureg.
        - windunit_spec: str,
                String specifying what unit to be used for wind speeds. Must 
                be found in ureg.
        - pressunit_spec: str,
                String specifying what unit to be used for pressure. Must 
                be found in ureg.

    Returns
        Tuple of strings for dataframe headers.
    """

    #create unit definitions and (sanitized) strings for headers
    temp_unit = ureg.Unit(tempunit_spec)
    wind_unit = ureg.Unit(windunit_spec)
    pressure_unit = ureg.Unit(pressunit_spec)
    tempunit_str = f"{temp_unit:~C}".replace("**","").replace("/","_").replace("*",".").replace("°","deg")
    windunit_str = f"{wind_unit:~C}".replace("**","").replace("/","_").replace("*",".")
    pressunit_str = f"{pressure_unit:~C}".replace("**","").replace("/","_").replace("*",".")

    #create headers
    temp_hdr = f"temp_{tempunit_str}"
    dpt_hdr = f"dpt_{tempunit_str}"
    ws_hdr = f"ws_{windunit_str}"
    wg_hdr = f"wg_{windunit_str}"
    press_hdr = f"pressure_{pressunit_str}"

    #filter which fields to included based on model possibilities.
    if model == "hrrr":
        header = ("time_UTC", 
                  press_hdr, 
                  temp_hdr, 
                  dpt_hdr, 
                  "rh_pct", 
                  ws_hdr, 
                  "wd_deg", 
                  wg_hdr, 
                  "precip_rate_kg_m2.s", 
                  "frozen_precip_pct", 
                  "precip_type")
    elif model == "gfs":
        header = ("time_UTC", 
                  press_hdr, 
                  temp_hdr, 
                  dpt_hdr, 
                  "rh_pct", 
                  ws_hdr, 
                  "wd_deg", 
                  wg_hdr, 
                  "precip_rate_kg_m2s", 
                  "frozen_precip_pct", 
                  "precip_type")
    elif model == "era5":
        header = ("time_UTC", 
                  press_hdr, 
                  temp_hdr, 
                  dpt_hdr, 
                  "rh_pct", 
                  ws_hdr, 
                  "wd_deg", 
                  wg_hdr, 
                  "precip_rate_kg_m2s", 
                  "frozen_precip_pct", 
                  "precip_type")
    elif model == "wrf":
        header = ("time_UTC", 
                  press_hdr, 
                  temp_hdr, 
                  dpt_hdr, 
                  "rh_pct", 
                  ws_hdr, 
                  "wd_deg", 
                  wg_hdr, 
                  "precip_rate_mm_s", 
                  "frozen_precip_pct", 
                  "precip_type")

    return header

def decode_metar(metar: str) -> dict:

    WX_CODES = {'VC': 'Vicinity',
                'RE': 'Recent',
                'MI': 'Shallow',
                'PR': 'Partial',
                'BC': 'Patches',
                'DR': 'Drifting',
                'BL': 'Blowing',
                'SH': 'Showers',
                'TS': 'Thunderstorms',
                'FZ': 'Freezing',
                'DZ': 'Drizzle',
                'RA': 'Rain',
                'SN': 'Snow',
                'SG': 'Snow Grains',
                'GS': 'Graupel',
                'GR': 'Hail',
                'PL': 'Ice Pellets',
                'IC': 'Ice Crystals',
                'UP': 'Unknown Precipitation',
                'FG': 'Fog',
                'BR': 'Mist',
                'HZ': 'Haze',
                'VA': 'Volcanic Ash',
                'DU': 'Dust',
                'FU': 'Smoke',
                'SA': 'Sand',
                'PY': 'Spray',
                'SQ': 'Squall',
                'PO': 'Dust/Sand Whirls',
                'DS': 'Duststorm',
                'SS': 'Sandstorm',
                'FC': 'Tornado',}

    WX_CODE_REGEXSTR = '[-+]?(' + '|'.join('(%s)' % code for code in WX_CODES) + ')+'

    CLD_CODES = {'SKC': 'Sky Clear',
                 'NCD': 'No Cloud Detected',
                 'CLR': 'Clear blw. FL120',
                 'NSC': 'No Significant Cloud',
                 'FEW': 'Few',
                 'SCT': 'Scattered',
                 'BKN': 'Broken',
                 'OVC': 'Overcast',}

    CLD_CODE_REGEXSTR = '(' + '|'.join('(%s)' % code for code in CLD_CODES) + ')(\d{3})?((TCU)|(CB))?'

    groups = metar.split(' ')

    apt_group = [str for str in groups if re.match('[KP]\w{3}', str)]
    utc_group = [str for str in groups if re.match('\d{6}Z', str)]
    wind_group = [str for str in groups if re.match('\d{5}(G\d{2})?[(KT)|(MPS)]', str)]
    wind_var_group = [str for str in groups if re.match('\d{3}V\d{3}', str)]
    wx_groups = [str for str in groups if re.match(WX_CODE_REGEXSTR, str)]
    cld_groups = [str for str in groups if re.match(CLD_CODE_REGEXSTR, str)]
    temp_group_auto = [str for str in groups if re.match('T\d+', str)]

    metar_dict = dict()

    if apt_group:
        metar_dict['Airport'] = apt_group[0]

    if utc_group:
        utc_group = utc_group[0]
        utc_day = int(utc_group[0:2])
        utc_hour = int(utc_group[2:4])
        utc_minute = int(utc_group[4:6])

        utc_time = f'{utc_day} {str(utc_hour).zfill(2)}:{str(utc_minute).zfill(2)}'
        metar_dict['TimeUTC'] = utc_time
        metar_dict['DayUTC'] = utc_day
        metar_dict['HourUTC'] = utc_hour
        metar_dict['MinuteUTC'] = utc_minute

    # if wind_group:
    #     wind_dir = int(wind_group[0][0:2])
    #     wind_speed = int(wind_group[0][3:4])

    if temp_group_auto:
        temp_group_auto = temp_group_auto[0]
        temp_c_str = f'{temp_group_auto[1:3]}.{temp_group_auto[4]}'
        temp_c = float(temp_c_str)
        dpt_c_str = f'{temp_group_auto[5:7]}.{temp_group_auto[8]}'
        dpt_c = float(dpt_c_str)

        metar_dict['TempAuto'] = temp_c
        metar_dict['TempAuto_units'] = 'degC'
        metar_dict['DptAuto'] = dpt_c
        metar_dict['DptAuto_units'] = 'degC'

    if wx_groups:
        metar_wx = list()
        for group in wx_groups:
            wx_str = group
            wx_str = wx_str.replace('-', 'Light ')
            wx_str = wx_str.replace('+', 'Heavy ')
            for wx_code in WX_CODES:
                wx_str = wx_str.replace(wx_code, WX_CODES[wx_code] + ' ')

            if wx_str[-1] == ' ':
                wx_str = wx_str[:-1]

            metar_wx.append(wx_str)

        metar_dict['WeatherGroups'] = metar_wx
        metar_dict['WeatherString'] = metar_wx.join(', ')

    return metar_dict


if __name__ == "__main__":

    parser = argparse.ArgumentParser(description=f'Create a meteogram data file from model data')
    parser.add_argument('-l', '--loc', 
                        help='specify a point to take the meteogram data from', 
                        nargs=2, 
                        type=float,
                        metavar=('lat', 'lon'),
                        dest='loc',
                        default=None)
    parser.add_argument('-tu', '--temperature-units',
                        help=f'specify the unit to be used for temperature/dewpoint (default: {DEFAULT_TEMP_UNITS})',
                        dest='temp_units',
                        type=str,
                        default=None)
    parser.add_argument('-wu', '--wind-units',
                        help=f'specify the unit to be used for wind speed/gusts (default: {DEFAULT_WIND_UNITS})',
                        dest='wind_units',
                        type=str,
                        default=None)
    parser.add_argument('-pu', '--pressure-units',
                        help=f'specify the unit to be used for pressure (default: {DEFAULT_PRESSURE_UNITS})',
                        dest='pressure_units',
                        type=str,
                        default=None)
    parser.add_argument('-vp', '--verification-plot',
                        help='generate a plot to verify location of meteogram',
                        dest='ver_plot',
                        action='store_true')
    parser.add_argument('-pt', '--record-precip-type',
                        help='record precip type as a string (default: off)',
                        dest='rec_pt',
                        action='store_true')
    parser.add_argument('-t', '--title',
                        help='specify a title for the output data csv',
                        type=str,
                        metavar='name',
                        dest='file_name',
                        default=None)
    parser.add_argument('model',  
                        help='specify data from what model is to be used. will error or produce undefined result if data format differs from specifed model',
                        choices=['hrrr', 'gfs', 'era5', 'wrf'],
                        type=str)
    parser.add_argument('input_file_directory',  
                        help='directory to read input files from',
                        type=str)
    parser.add_argument('save_directory',  
                        help='directory to save meteogram to',
                        type=str)

    args = parser.parse_args()
    print(args)

    ## Quality Checks/Sanitizing of input arguments

    if args.input_file_directory[-1:] != "/":
        input_dir = args.input_file_directory + "/"
    else:
        input_dir = args.input_file_directory

    print(f"Reading files from '{input_dir}'")

    if args.save_directory[-1:] != "/":
        save_dir = args.save_directory + "/"
    else:
        save_dir = args.save_directory

    if args.file_name:
        save_path = os.path.join(save_dir, f'{args.file_name}.csv')
    else:
        save_path = os.path.join(save_dir, f'MeteogramData_{args.model.upper()}.csv')

    print(f"Saving to '{save_path}'")

    if args.model in NETCDF_MODELS:
        input_files = sorted(glob.glob(f'{input_dir}*'))
    elif args.model in GRIB2_MODELS:
        input_files = sorted(glob.glob(f'{input_dir}*.grib2'))

        if not input_files:
            input_files = sorted(glob.glob(f'{input_dir}*'))

    if not input_files:
        raise FileNotFoundError("No input files found in specified directory.")
    else:
        num_input_files = len(input_files)
        print(f"Successfully found {num_input_files} files.")

    unit_registry = pint.UnitRegistry()

    if args.temp_units:
        if str(unit_registry(args.temp_units).dimensionality) == '[temperature]':
            temp_units = args.temp_units
        else:
            warnings.warn(f'Cannot use unit {args.temp_units} for temperature, wrong dimensionality. Using default unit {DEFAULT_TEMP_UNITS}')
            temp_units = DEFAULT_TEMP_UNITS
    else:
        temp_units = DEFAULT_TEMP_UNITS

    if args.wind_units:
        if str(unit_registry(args.wind_units).dimensionality) == '[length] / [time]':
            wind_units = args.wind_units
        else:
            warnings.warn(f'Cannot use unit {args.wind_units} for wind, wrong dimensionality. Using default unit {DEFAULT_WIND_UNITS}')
            wind_units = DEFAULT_WIND_UNITS
    else:
        wind_units = DEFAULT_WIND_UNITS

    if args.pressure_units:
        if str(unit_registry(args.pressure_units).dimensionality) == '[mass] / [length] / [time] ** 2':
            pressure_units = args.pressure_units
        else:
            warnings.warn(f'Cannot use unit {args.pressure_units} for pressure, wrong dimensionality. Using default unit {DEFAULT_PRESSURE_UNITS}')
            pressure_units = DEFAULT_PRESSURE_UNITS
    else:
        pressure_units = DEFAULT_PRESSURE_UNITS

    if args.rec_pt:
        warnings.warn("Precipitation type recording not yet implemented.")

    if args.ver_plot:
        print("Generating requested verification plot...")
        generate_verification_plot(meteogram_lat, meteogram_lon, save_dir)

    with tqdm(miniters=0, total=num_input_files, desc=f'Making meteogram from {args.model.upper()} data...', ascii=" ✎✏✐█") as progress:
        data_rows = []
        for input_file_path in input_files:

            if args.model == "hrrr":
                #print(data_rows)
                row_data = parse_hrrr_data(input_file_path,
                                           meteogram_lat,
                                           meteogram_lon,
                                           unit_registry,
                                           temp_units,
                                           wind_units,
                                           pressure_units)
            elif args.model == "wrf":
                row_data = parse_wrf_data(input_file_path,
                                           meteogram_lat,
                                           meteogram_lon,
                                           unit_registry,
                                           temp_units,
                                           wind_units,
                                           pressure_units)

            data_rows += [row_data]
            
            progress.update()

    meteogram_df = pd.DataFrame(data_rows, columns=generate_header(args.model, unit_registry, temp_units, wind_units, pressure_units))

    meteogram_df.to_csv(save_path, index=False)

    print("Done!")






