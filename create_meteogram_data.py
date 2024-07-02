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

def parse_hrrr_data(file_path: str, 
                    lat: float, 
                    lon: float, 
                    ureg: pint.UnitRegistry, 
                    tempunit_str: str, 
                    windunit_str: str, 
                    pressunit_str: str, 
                    **kwargs) -> list:

    """
    Using an input HRRR data file, parse the file and extract the values of
    certain variables at a specified point. Units for temperature, wind, and
    pressure can be specified.

    Inputs
        - file_path, str:
                Path to HRRR data file.
        - lat, float: 
                Latitude of requested point.
        - lon, float: 
                Longitude of requested point.
        - ureg, pint.UnitRegistry:
                A UnitRegistry object to be used to create unit identifiers.
        - tempunit_str: str,
                String specifying what unit to be used for temperature. Must 
                be found in ureg.
        - windunit_str: str,
                String specifying what unit to be used for wind speeds. Must 
                be found in ureg.
        - pressunit_str: str,
                String specifying what unit to be used for pressure. Must 
                be found in ureg.

    Returns
        List of data from this point, corresponding to the headers generated
        via generate_headers function.
    """

    row_data = []

    sfc_inst_data = xr.open_dataset(file_path, engine="cfgrib", filter_by_keys={'typeOfLevel': 'surface', 'stepType': 'instant'})
    m2_data = xr.open_dataset(file_path, engine="cfgrib", filter_by_keys={'typeOfLevel': 'heightAboveGround', 'level': 2})
    m10_data = xr.open_dataset(file_path, engine="cfgrib", filter_by_keys={'typeOfLevel': 'heightAboveGround', 'level': 10})
    
    # Becuase the grib reader reads longitude from 0-360 and not -180-180
    # we have to adjust the `lon`.
    if lon < 0:
        lon += 360
    
    # Find the nearest neightbor from the maximum absolute differences...
    abslat = np.abs(sfc_inst_data.latitude-lat)
    abslon = np.abs(sfc_inst_data.longitude-lon)

    c = np.maximum(abslon, abslat)
    
    # Get the index of the nearest point
    # (Again, the values are backwards you might expect because the 
    # coordinates are in (y, x) order.)
    ([idx_y], [idx_x]) = np.where(c == np.min(c))
    
    #Now that we have the index location of the nearest point to the requested lat/lon, 
    #we can use the select function sel to get all the data in the height coordinate at a single point in the dataset.
    
    point_sfc_inst = sfc_inst_data.sel(y=idx_y, x=idx_x)
    point_2m = m2_data.sel(y=idx_y, x=idx_x)
    point_10m = m10_data.sel(y=idx_y, x=idx_x)

    # print(f"SFC pt: {point_sfc_inst.latitude.values}, {point_sfc_inst.longitude.values-360}")
    # print(f"2m pt: {point_2m.latitude.values}, {point_2m.longitude.values-360}")
    # print(f"10m pt: {point_10m.latitude.values}, {point_10m.longitude.values-360}")

    time_utc = pd.to_datetime(point_sfc_inst.time.values)
    row_data += [time_utc.strftime("%Y-%m-%d %H:%M:%S")]
    #print(time_utc)
    
    pressure_pa = units.Quantity(point_sfc_inst.sp.values, 'pascal')
    pressure_conv = pressure_pa.to(ureg(pressunit_str))
    row_data += [pressure_conv.magnitude]
    #print(pressure_hPa)
    
    temp_k = units.Quantity(point_2m.pt.values, 'K')
    temp_conv = temp_k.to(ureg(tempunit_str))
    row_data += [temp_conv.magnitude]

    dpt_k = units.Quantity(point_2m.d2m.values, 'K')
    dpt_conv = dpt_k.to(ureg(tempunit_str))
    row_data += [dpt_conv.magnitude]

    rh = units.Quantity(point_2m.r2.values, '%')
    row_data += [rh.magnitude]
    # print(rh)

    u = units.Quantity(point_10m.u10.values, 'm/s')
    v = units.Quantity(point_10m.u10.values, 'm/s')
    ws_ms = mpcalc.wind_speed(u, v)
    ws_conv = ws_ms.to(ureg(windunit_str))
    row_data += [ws_conv.magnitude]

    wd_deg = mpcalc.wind_direction(u, v)
    row_data += [wd_deg.magnitude]

    # print(ws_kts)
    # print(wd_deg)

    print(point_sfc_inst.gust)
    wg_ms = units.Quantity(point_sfc_inst.gust.values, 'm/s')
    wg_conv = wg_ms.to(ureg(windunit_str))
    row_data += [wg_conv.magnitude]

    precip_rate = units.Quantity(point_sfc_inst.prate.values, "kg m**-2 s**-1")
    row_data += [precip_rate.magnitude]
    # print(precip_rate)

    pct_frz_precip = units.Quantity(point_sfc_inst.cpofp.values, "%")
    row_data += [pct_frz_precip.magnitude]
    # print(pct_frz_precip)

    if kwargs.get('precip_type'):

        cat_ra = point_sfc_inst.crain.values
        cat_fzrn = point_sfc_inst.cfrzr.values
        cat_ip = point_sfc_inst.cicep.values
        cat_sn = point_sfc_inst.csnow.values

        model_ptypes_name = ["Rain", "Freezing Rain", "Ice Pellets", "Snow"]
        model_ptypes_present = [cat_ra, cat_fzrn, cat_ip, cat_sn]

        precip_type = "None"
        for pval, pname in zip(model_ptypes_present, model_ptypes_name):
            # print(f"{pname}: {pval}")
            if pval == 1:
                precip_type += f"{pname}, "

        if precip_type != "None":
            incl_precip_type = precip_type[:-1].replace("None", "")
            incl_pct_frz_precip = pct_frz_precip.magnitude
            incl_precip_rate = precip_rate.magnitude
        else:
            incl_precip_type = precip_type
            incl_pct_frz_precip = 0
            incl_precip_rate = 0

        # print(incl_precip_type, incl_precip_rate, incl_pct_frz_precip)

        row_data += [incl_precip_type]

    else:
        row_data += [MISSING_FLAG]

    #print(row_data)
    
    return row_data

def parse_gfs_data(file_path, **kwargs):
    raise NotImplementedError("Parsing of GFS data not yet implemented.")

def parse_era5_data(file_path, **kwargs):
    raise NotImplementedError("Parsing of ERA5 data not yet implemented.")

def parse_wrf_data(file_path: str, 
                    lat: float, 
                    lon: float, 
                    ureg: pint.UnitRegistry, 
                    tempunit_str: str, 
                    windunit_str: str, 
                    pressunit_str: str, 
                    **kwargs) -> list:
    """
    Using an input HRRR data file, parse the file and extract the values of
    certain variables at a specified point. Units for temperature, wind, and
    pressure can be specified.

    Inputs
        - file_path, str:
                Path to HRRR data file.
        - lat, float: 
                Latitude of requested point.
        - lon, float: 
                Longitude of requested point.
        - ureg, pint.UnitRegistry:
                A UnitRegistry object to be used to create unit identifiers.
        - tempunit_str: str,
                String specifying what unit to be used for temperature. Must 
                be found in ureg.
        - windunit_str: str,
                String specifying what unit to be used for wind speeds. Must 
                be found in ureg.
        - pressunit_str: str,
                String specifying what unit to be used for pressure. Must 
                be found in ureg.

    Returns
        List of data from this point, corresponding to the headers generated
        via generate_headers function.
    """

    row_data = []

    wrf_data = xr.open_dataset(file_path, engine="netcdf4")
    latitude = getvar(Dataset(file_path), "lat")
    longitude = getvar(Dataset(file_path), "lon")
    
    # Find the nearest neightbor from the maximum absolute differences...
    abslat = np.abs(latitude-lat)
    abslon = np.abs(longitude-lon)

    c = np.maximum(abslon, abslat)
    
    # Get the index of the nearest point
    # (Again, the values are backwards you might expect because the 
    # coordinates are in (y, x) order.)
    ([idx_y], [idx_x]) = np.where(c == np.min(c))
    
    #Now that we have the index location of the nearest point to the requested lat/lon, 
    #we can use the select function sel to get all the data in the height coordinate at a single point in the dataset.
    
    point_data = wrf_data.sel(south_north=idx_y, west_east=idx_x)
    #print(point_data.KPBL)
    # print(point_data.PBLH.values)
    # print(point_data.PBLH)
    # print(point_data.KPBL)

    time_utc = pd.to_datetime(point_data.XTIME.values[0])
    row_data += [time_utc.strftime("%Y-%m-%d %H:%M:%S")]
    
    pressure_pa = units.Quantity(point_data.PSFC.values, 'pascal')
    pressure_conv = pressure_pa.to(ureg(pressunit_str))
    row_data += [pressure_conv.magnitude[0]]
    
    temp_k = units.Quantity(point_data.T2.values, 'K')
    temp_conv = temp_k.to(ureg(tempunit_str))
    row_data += [temp_conv.magnitude[0]]

    rh = units.Quantity(point_data.RH02.values * 100, '%')

    dpt_c = mpcalc.dewpoint_from_relative_humidity(temp_k.to(ureg('degC')), rh)
    if np.isnan(dpt_c):
        dpt_conv = MISSING_FLAG
        row_data += [dpt_conv]
    else:
        dpt_conv = dpt_c.to(ureg(tempunit_str))
        row_data += [dpt_conv.magnitude[0]]

    row_data += [rh.magnitude[0]]

    u = units.Quantity(point_data.U10.values, 'm/s')
    v = units.Quantity(point_data.V10.values, 'm/s')
    ws_ms = mpcalc.wind_speed(u, v)
    ws_conv = ws_ms.to(ureg(windunit_str))
    row_data += [ws_conv.magnitude[0]]

    wd_deg = mpcalc.wind_direction(u, v)
    row_data += [wd_deg.magnitude[0]]

    '''
    Below is the post-processing algorithm implementation of 
    wind gust speeds, used by the HRRR model. This algorithm
    is sourced from this document:
    https://rapidrefresh.noaa.gov/Diag-vars-NOAA-TechMemo.pdf

    This implementation uses the following equation to calculate
    wind gust (technically 'gust potential') speed:

            gust = wsfc + max(weight_coef * [w(lvl) - wsfc])

    where wsfc is the surface wind speed (in this case @ 10m)
    in m/s, w(lvl) is the wind speed at that level in m/s, and
    weight_coef is the weight coefficent calculated for that
    level. The difference between the wind speed at that level
    and the surface wind speed is calculated (known as the 
    "excess wind speed"), which is then multiplied by the weight
    coefficent for that level. The weight coefficent ranges from
    1.0 at the surface to 0.5 at 1 km (in increments of 1m here),
    and is 0.5 at > 1 km. This is done for all levels below the
    top of the PBL, and the max of all the results is taken and 
    added to the surface wind speed to get the "wind gust potential".
    '''

    #U/V column winds must be destaggered off of mass grid
    wrf_uwinds = wrf.destagger(point_data.U, 2, meta=True)
    point_uwinds = wrf_uwinds.sel(west_east=idx_x)

    wrf_vwinds = wrf.destagger(point_data.V, 2, meta=True)
    point_vwinds = wrf_vwinds.sel(south_north=idx_y)

    col_winds_ms = mpcalc.wind_speed(units.Quantity(point_uwinds.data, 'm/s'), 
                                     units.Quantity(point_vwinds.data, 'm/s'))

    gph = units.Quantity(point_data.PHP.values, 'm**2/s**2')
    hgt_m = mpcalc.geopotential_to_height(gph)
    hgt_m = hgt_m[0]

    pbl_hgt_m = units.Quantity(point_data.PBLH.values[0], 'm')

    closest_lev_to_pbl = np.abs(hgt_m[hgt_m < pbl_hgt_m] - pbl_hgt_m) #This only looks at levels below PBL top

    weight_coefs_sfcto1km_per1m = np.linspace(1, 0.5, 1000)

    if closest_lev_to_pbl.size > 0:
        closest_lev_to_pbl = closest_lev_to_pbl.argmin()
        closest_lev_hgt_to_pbl = hgt_m[closest_lev_to_pbl]

        winds_below_pbl = col_winds_ms[0][0:closest_lev_to_pbl]

        excess_wind_speeds = winds_below_pbl - ws_ms

        weight_coefs_list = []
        for item in winds_below_pbl:
            item_idx = np.where(winds_below_pbl == item)[0][0]
            ws_hgt = hgt_m[item_idx] #use the index of the wind speed to get the height of that level

            if ws_hgt.magnitude > 1000:
                weight_coefs_list.append(0.5)
            else:
                weight_coefs_list.append(weight_coefs_sfcto1km_per1m[int(round(ws_hgt.magnitude))])

        weight_coefs = np.asarray(weight_coefs_list)
        
        excess_ws_weighted = excess_wind_speeds * weight_coefs
        max_excess_ws = max(excess_ws_weighted)

        wg_ms = units.Quantity(ws_ms.magnitude + max_excess_ws.magnitude, 'm/s')
        wg_conv = wg_ms.to(ureg(windunit_str))
        row_data += [wg_conv.magnitude[0]]

    else:
        if USE_MISSING_FOR_WRF_GUSTS:
            row_data += [MISSING_FLAG]
        else:
            wg_ms = ws_ms
            wg_conv = wg_ms.to(ureg(windunit_str))
            row_data += [wg_conv.magnitude[0]]

    precip_rate = units.Quantity(point_data.INST_PRATE.values, "mm s**-1")
    row_data += [precip_rate.magnitude[0]]

    pct_frz_precip = units.Quantity(point_data.SR.values, "%")
    row_data += [pct_frz_precip.magnitude[0]]

    if kwargs.get('precip_type'):
        raise NotImplementedError("Precip type not yet implemented for WRF.")

        # cat_ra = point_sfc_inst.crain.values
        # cat_fzrn = point_sfc_inst.cfrzr.values
        # cat_ip = point_sfc_inst.cicep.values
        # cat_sn = point_sfc_inst.csnow.values

        # model_ptypes_name = ["Rain", "Freezing Rain", "Ice Pellets", "Snow"]
        # model_ptypes_present = [cat_ra, cat_fzrn, cat_ip, cat_sn]

        # precip_type = "None"
        # for pval, pname in zip(model_ptypes_present, model_ptypes_name):
        #     # print(f"{pname}: {pval}")
        #     if pval == 1:
        #         precip_type += f"{pname}, "

        # if precip_type != "None":
        #     incl_precip_type = precip_type[:-1].replace("None", "")
        #     incl_pct_frz_precip = pct_frz_precip.magnitude
        #     incl_precip_rate = precip_rate.magnitude
        # else:
        #     incl_precip_type = precip_type
        #     incl_pct_frz_precip = 0
        #     incl_precip_rate = 0

        # # print(incl_precip_type, incl_precip_rate, incl_pct_frz_precip)

        # row_data += [incl_precip_type]

    else:
        row_data += [MISSING_FLAG]

    #print(row_data)
    
    return row_data


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

    if args.loc:
        meteogram_lat = args.loc[0]
        meteogram_lon = args.loc[1]
        if (meteogram_lat > 90) or (meteogram_lat < -90):
            raise ValueError("Invalid latitude for meteogram point, must +90 =< x =< -90")
        if (meteogram_lon > 180) or (meteogram_lon < -180):
            raise ValueError("Invalid longitude for meteogram point, must +180 =< x =< -180")
    else:
        raise ValueError("A location must be specified using the -l flag.")

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






