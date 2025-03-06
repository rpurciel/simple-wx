## The file at the following link lists units that can be used in this context:
## https://github.com/hgrecco/pint/blob/master/pint/default_en.txt
## Note: A dimensionality check is automatically performed on any specified units

## Default units for pressure values
## (dimensionality: '[mass] / [length] / [time] ** 2')
DEFAULT_PRESSURE_UNITS = 'hectopascal'

## Optional unit 
## Options are: degC, deg
DEFAULT_TEMP_UNITS = 'degC'

## Default units for wind speed/wind gust values
## (dimensionality: '[length] / [time]')
DEFAULT_WIND_UNITS = 'kt'

DEFAULT_WRF_INPUT_FREQUENCY = '5min'






'''


            STATIC CODE BELOW HERE


'''

import os
import glob
from datetime import datetime
import argparse
import csv
import warnings
warnings.filterwarnings("ignore")

import xarray as xr
import numpy as np
import pandas as pd
import metpy
from metpy.units import units
import metpy.calc as mpcalc
from tqdm.auto import tqdm
import pint
from wrf import getvar
from netCDF4 import Dataset

from __internal_funcs import (plot_towns, generate_verification_plot, remove_files)

np.seterr(divide = 'ignore', invalid='ignore')

MISSING_FLAG = -999

NETCDF_MODELS = ['wrf']

GRIB2_MODELS = ['hrrr', 
                'gfs', 
                'era5']

def destag_variable(stag_var: xr.DataArray) -> xr.DataArray:

    from wrf import destagger

    dims = stag_var.dims
    for dim in dims:
        if "stag" in dim:
            stag_dim_idx = dims.index(dim)
            
    destag = destagger(stag_var, 
                       stag_dim_idx,
                       meta=True)

    destag = destag.assign_coords(stag_var.coords)

    return destag

def get_file_time(input_file: str,
                  model: str) -> datetime:

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

    if model in NETCDF_MODELS:
        data = xr.open_dataset(input_file, engine="netcdf4")
    else:
        data = xr.open_dataset(input_file, engine="cfgrib", filter_by_keys={'typeOfLevel': 'isobaricInhPa'})

    if (model == "hrrr") or (model == "era5"):
        return pd.Timestamp(data.time.values).to_pydatetime()
    elif model == "wrf":
        return pd.Timestamp(data.XTIME.values[0]).to_pydatetime()

def create_file_times_dict(possible_input_files: list[str, ...],
                           model: str) -> dict:

    file_times_lookup_table = {}

    for file in possible_input_files:
        if ".idx" in file:
            continue

        if model in NETCDF_MODELS:
            data = xr.open_dataset(file, engine="netcdf4")
        else:
            data = xr.open_dataset(file, engine="cfgrib", filter_by_keys={'typeOfLevel': 'isobaricInhPa'})

        if (model == "hrrr") or (model == "era5"):
            file_time = pd.Timestamp(data.time.values).to_pydatetime()
        elif model == "wrf":
            file_time = pd.Timestamp(data.XTIME.values[0]).to_pydatetime()

        file_times_lookup_table.update({file_time.strftime("%Y-%m-%d %H:%M:%S"): file})

    return file_times_lookup_table

def parse_file_times(file_times_lookup_table: dict,
                     point_time: datetime,
                     data_file_freq: str) -> str:

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

    nearest_time_pd = pd.Timestamp(point_time).round(freq=data_file_freq).to_pydatetime()
    nearest_time = datetime.strftime(nearest_time_pd, "%Y-%m-%d %H:%M:%S")

    return file_times_lookup_table[nearest_time]

def parse_hrrr_data(file_path: str,
                    save_dir: str, 
                    lat: float, 
                    lon: float, 
                    point_time: datetime,
                    point_index: int = 1,
                    previous_points: list[tuple[float, float], ...] = None,
                    ureg: pint.UnitRegistry = pint.UnitRegistry(),
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

    col_data = xr.open_dataset(file_path, engine="cfgrib", filter_by_keys={'typeOfLevel': 'isobaricInhPa'})
    sfc_data = xr.open_dataset(file_path, engine="cfgrib", filter_by_keys={'typeOfLevel': 'surface', 'stepType': 'instant'})
    m2_data = xr.open_dataset(file_path, engine="cfgrib", filter_by_keys={'typeOfLevel': 'heightAboveGround', 'level': 2})
    m10_data = xr.open_dataset(file_path, engine="cfgrib", filter_by_keys={'typeOfLevel': 'heightAboveGround', 'level': 10})
    
    # Becuase the grib reader reads longitude from 0-360 and not -180-180
    # we have to adjust the `lon`.
    orig_lon = lon
    if lon < 0:
        lon += 360
    
    # Find the nearest neightbor from the maximum absolute differences...
    abslat = np.abs(col_data.latitude-lat)
    abslon = np.abs(col_data.longitude-lon)

    c = np.maximum(abslon, abslat)
    
    # Get the index of the nearest point
    # (Again, the values are backwards you might expect because the 
    # coordinates are in (y, x) order.)
    ([idx_y], [idx_x]) = np.where(c == np.min(c))
    
    #Now that we have the index location of the nearest point to the requested lat/lon, 
    #we can use the select function sel to get all the data in the height coordinate at a single point in the dataset.
    
    point_col = col_data.sel(y=idx_y, x=idx_x)
    point_sfc = sfc_data.sel(y=idx_y, x=idx_x)
    point_2m = m2_data.sel(y=idx_y, x=idx_x)
    point_10m = m10_data.sel(y=idx_y, x=idx_x)

    sel_pt_lat = float(point_col.latitude.data)
    sel_pt_lon = abs(float(point_col.longitude.data - 360))

    if previous_points is not None:
        if (sel_pt_lat, sel_pt_lon) in previous_points:
            return None, previous_points
        else:
            previous_points.append((sel_pt_lat, sel_pt_lon))

    time_utc = pd.to_datetime(point_col.time.values)
    date_str = datetime.strftime(time_utc, "%Y-%m-%d %H:%M:%S")

    col_press_hpa = units.Quantity(point_col.isobaricInhPa.data, 'hectopascal')

    col_temp_k = units.Quantity(point_col.t.data, 'K')
    col_dpt_k = units.Quantity(point_col.dpt.data, 'K')
    col_temp_c = col_temp_k.to(ureg('degC'))
    col_dpt_c = col_dpt_k.to(ureg('degC'))

    col_u_ms = units.Quantity(point_col.u.data, 'm/s')
    col_v_ms = units.Quantity(point_col.v.data, 'm/s')
    col_u_kts = col_u_ms.to(ureg('kts'))
    col_v_kts = col_v_ms.to(ureg('kts'))

    print(point_col.gh)

    col_geopot = units.Quantity(point_col.gh.data, 'm**2/s**2')
    col_layer_height = mpcalc.geopotential_to_height(col_geopot)

    sfc_press_pa = units.Quantity(point_sfc.sp.data, 'pascal')
    sfc_press_hpa = sfc_press_pa.to(ureg('hectopascal'))

    sfc_temp_k = units.Quantity(point_2m.t2m.data, 'K')
    sfc_dpt_k = units.Quantity(point_2m.d2m.data, 'K')
    sfc_temp_c = sfc_temp_k.to(ureg('degC'))
    sfc_dpt_c = sfc_dpt_k.to(ureg('degC'))

    sfc_u_ms = units.Quantity(point_10m.u10.data, 'm/s')
    sfc_v_ms = units.Quantity(point_10m.v10.data, 'm/s')
    sfc_u_kts = sfc_u_ms.to(ureg('kts'))
    sfc_v_kts = sfc_v_ms.to(ureg('kts'))

    sfc_geopot = units.Quantity(point_sfc.orog.data, 'm**2/s**2')
    print(point_sfc.orog)
    sfc_height = mpcalc.geopotential_to_height(sfc_geopot)

    surface_data = pd.DataFrame(data=[sfc_press_hpa.magnitude,
                                      sfc_temp_c.magnitude,
                                      sfc_dpt_c.magnitude,
                                      sfc_u_kts.magnitude,
                                      sfc_v_kts.magnitude,
                                      sfc_geopot.magnitude])

    column_data = pd.DataFrame(data=[col_press_hpa.magnitude,
                                     col_temp_c.magnitude,
                                     col_dpt_c.magnitude,
                                     col_u_kts.magnitude,
                                     col_v_kts.magnitude,
                                     col_geopot.magnitude])

    #remove data that is below the surface level
    column_data = column_data.T
    column_data = column_data[column_data[0] < sfc_press_hpa.magnitude]

    sounding_data = pd.concat([surface_data.T, column_data],
                              axis=0,
                              ignore_index=True)

    if sel_pt_lat >= 0:
        lat_hemisphere = 'N'
    else:
        lat_hemisphere = 'S'

    if sel_pt_lon >= 0:
        lon_hemisphere = 'W'
    else:
        lon_hemisphere = 'E'

    internal_name = f'Pt: {round(lat, 4)} {round(orig_lon, 4)} Gd: {round(sel_pt_lat, 4)} {round(sel_pt_lon, 4)} Md: HRRR (No. {point_index})'

    header = {0:['RAOB/CSV','DTG','LAT','LON','ELEV','MOISTURE','WIND','GPM','MISSING','RAOB/DATA','PRES'],
              1:[internal_name,date_str,np.abs(sel_pt_lat),np.abs(sel_pt_lon),sfc_geopot.magnitude,'TD','kts','MSL',-999,'','TEMP'],
              2:['','',lat_hemisphere,lon_hemisphere,'m','','U/V','','','','TD'],
              3:['','','','','','','','','','','UU'],
              4:['','','','','','','','','','','VV'],
              5:['','','','','','','','','','','GPM']}

    header_data = pd.DataFrame(data=header)

    final_data = pd.concat([header_data, sounding_data],
                           axis=0,
                           ignore_index=True)

    file_name = f'HRRR_RAOB_pt{point_index}_{datetime.strftime(point_time, "%Y%m%d_%H%M%S")}_frm{datetime.strftime(time_utc, "%Y%m%d_%H%M%S")}.csv'

    if not os.path.exists(save_dir):
        os.makedirs(save_dir)

    final_data.to_csv(os.path.join(save_dir, file_name), index=False, header=False)

    if previous_points is not None:
        return (sel_pt_lat, sel_pt_lon, point_index), previous_points
    else:
        return (sel_pt_lat, sel_pt_lon, point_index)

def parse_gfs_data(file_path, **kwargs):
    raise NotImplementedError("Parsing of GFS data not yet implemented.")

def parse_era5_data(file_path: str,
                    save_dir: str, 
                    lat: float, 
                    lon: float, 
                    point_time: datetime,
                    point_index: int = 1,
                    previous_points: list[tuple[float, float], ...] = None,
                    ureg: pint.UnitRegistry = pint.UnitRegistry(),
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

    col_data = xr.open_dataset(file_path, engine="cfgrib", drop_variables='i10fg', filter_by_keys={'typeOfLevel': 'isobaricInhPa'})

    sfc_file_path = file_path.replace("UA", "SFC")
    if os.path.isfile(sfc_file_path):
        sfc_data = xr.open_dataset(sfc_file_path, drop_variables='i10fg', engine="cfgrib")
        use_sfc_data = True
    else:
        warnings.warn("Surface data was not found for ERA5. Absent surface data degrades sounding quality.")
        use_sfc_data = False
    
    # Becuase the grib reader reads longitude from 0-360 and not -180-180
    # we have to adjust the `lon`.
    orig_lon = lon

    point_col = col_data.sel(longitude=lon, latitude=lat, method='nearest')

    sel_pt_lat = float(point_col.latitude.data)
    sel_pt_lon = float(point_col.longitude.data)

    if previous_points is not None:
        if (sel_pt_lat, sel_pt_lon) in previous_points:
            return None, previous_points
        else:
            previous_points.append((sel_pt_lat, sel_pt_lon))

    time_utc = pd.to_datetime(point_col.time.values)
    date_str = datetime.strftime(time_utc, "%Y-%m-%d %H:%M:%S")

    col_press_hpa = units.Quantity(point_col.isobaricInhPa.data, 'hectopascal')

    col_temp_k = units.Quantity(point_col.t.data, 'K')
    col_temp_c = col_temp_k.to(ureg('degC'))

    col_rh_pct = units.Quantity(point_col.r.data, '%')

    col_dpt_c = mpcalc.dewpoint_from_relative_humidity(col_temp_c,
                                                       col_rh_pct)

    col_u_ms = units.Quantity(point_col.u.data, 'm/s')
    col_v_ms = units.Quantity(point_col.v.data, 'm/s')
    col_u_kts = col_u_ms.to(ureg('kts'))
    col_v_kts = col_v_ms.to(ureg('kts'))

    col_geopot = units.Quantity(point_col.z.data, 'm**2/s**2')
    col_layer_height = mpcalc.geopotential_to_height(col_geopot)

    column_data = pd.DataFrame(data=[col_press_hpa.magnitude,
                                     col_temp_c.magnitude,
                                     col_dpt_c.magnitude,
                                     col_u_kts.magnitude,
                                     col_v_kts.magnitude,
                                     col_layer_height.magnitude])

    if use_sfc_data:

        point_sfc = sfc_data.sel(longitude=lon, latitude=lat, method='nearest')

        sfc_press_pa = units.Quantity(point_sfc.sp.data, 'pascal')
        sfc_press_hpa = sfc_press_pa.to(ureg('hectopascal'))

        sfc_temp_k = units.Quantity(point_sfc.t2m.data, 'K')
        sfc_dpt_k = units.Quantity(point_sfc.d2m.data, 'K')
        sfc_temp_c = sfc_temp_k.to(ureg('degC'))
        sfc_dpt_c = sfc_dpt_k.to(ureg('degC'))

        sfc_u_ms = units.Quantity(point_sfc.u10.data, 'm/s')
        sfc_v_ms = units.Quantity(point_sfc.v10.data, 'm/s')
        sfc_u_kts = sfc_u_ms.to(ureg('kts'))
        sfc_v_kts = sfc_v_ms.to(ureg('kts'))

        sfc_geopot = units.Quantity(point_sfc.z.data, 'm**2/s**2')
        sfc_layer_height = mpcalc.geopotential_to_height(sfc_geopot)

        surface_data = pd.DataFrame(data=[sfc_press_hpa.magnitude,
                                          sfc_temp_c.magnitude,
                                          sfc_dpt_c.magnitude,
                                          sfc_u_kts.magnitude,
                                          sfc_v_kts.magnitude,
                                          sfc_layer_height.magnitude])

        bot_lvl = sfc_layer_height.magnitude

        column_data = column_data.T
        #remove data that is below the surface level
        column_data = column_data[column_data[0] < sfc_press_hpa.magnitude]

        sounding_data = pd.concat([surface_data.T, column_data],
                                   axis=0,
                                   ignore_index=True)
    else:
        bot_lvl = col_layer_height[0].magnitude

        sounding_data = column_data.T

    if sel_pt_lat >= 0:
        lat_hemisphere = 'N'
    else:
        lat_hemisphere = 'S'

    if sel_pt_lon >= 0:
        lon_hemisphere = 'W'
    else:
        lon_hemisphere = 'E'

    internal_name = f'Pt: {round(lat, 4)} {round(orig_lon, 4)} Gd: {round(sel_pt_lat, 4)} {round(sel_pt_lon, 4)} Md: ERA5 (No. {point_index})'

    header = {0:['RAOB/CSV','DTG','LAT','LON','ELEV','MOISTURE','WIND','GPM','MISSING','RAOB/DATA','PRES'],
              1:[internal_name,date_str,np.abs(sel_pt_lat),np.abs(sel_pt_lon),bot_lvl,'TD','kts','MSL',-999,'','TEMP'],
              2:['','',lat_hemisphere,lon_hemisphere,'m','','U/V','','','','TD'],
              3:['','','','','','','','','','','UU'],
              4:['','','','','','','','','','','VV'],
              5:['','','','','','','','','','','GPM']}

    header_data = pd.DataFrame(data=header)

    final_data = pd.concat([header_data, sounding_data],
                           axis=0,
                           ignore_index=True)

    file_name = f'ERA5_RAOB_pt{point_index}_{datetime.strftime(point_time, "%Y%m%d_%H%M%S")}_frm{datetime.strftime(time_utc, "%Y%m%d_%H%M%S")}.csv'

    if not os.path.exists(save_dir):
        os.makedirs(save_dir)

    final_data.to_csv(os.path.join(save_dir, file_name), index=False, header=False)

    if previous_points is not None:
        return (sel_pt_lat, sel_pt_lon, point_index), previous_points
    else:
        return (sel_pt_lat, sel_pt_lon, point_index)

def parse_wrf_data(file_path: str,
                   save_dir: str, 
                   lat: float, 
                   lon: float, 
                   point_time: datetime,
                   point_index: int = 1,
                   previous_points: list[tuple[float, float], ...] = None,
                   ureg: pint.UnitRegistry = pint.UnitRegistry(),
                   **kwargs) -> list:

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

    # print(np.where(c == np.min(c)))
    (idx_y, idx_x) = np.where(c == np.min(c))

    idx_y = idx_y[0] if len(idx_y) > 1 else idx_y
    idx_x = idx_x[0] if len(idx_x) > 1 else idx_x

    wrf_data['U'] = destag_variable(wrf_data.U)
    wrf_data['V'] = destag_variable(wrf_data.V)
    wrf_data['PH'] = destag_variable(wrf_data.PH)
    wrf_data['PHB'] = destag_variable(wrf_data.PHB)

    point_data = wrf_data.sel(south_north=idx_y.flatten(), west_east=idx_x.flatten())

    sel_pt_lat = float(point_data.XLAT.values[0])
    sel_pt_lon = float(point_data.XLONG.values[0])

    if previous_points is not None:
        if (sel_pt_lat, sel_pt_lon) in previous_points:
            return None, previous_points
        else:
            previous_points.append((sel_pt_lat, sel_pt_lon))

    time_utc = pd.to_datetime(point_data.XTIME.values[0])
    date_str = datetime.strftime(time_utc, "%Y-%m-%d %H:%M:%S")

    col_press = point_data.P.data + point_data.PB.data
    col_press = col_press.flatten()

    col_press_pa = units.Quantity(col_press, 'pascal')
    col_press_hpa = col_press_pa.to(ureg('hectopascal'))

    col_theta_k = units.Quantity(point_data.T.data.flatten() + 300, 'K')
    col_theta_c = col_theta_k.to(ureg('degC'))

    col_temp_k = mpcalc.temperature_from_potential_temperature(col_press_hpa,
                                                             col_theta_k)

    col_temp_c = col_temp_k.to(ureg('degC'))

    # MOLECULAR_WEIGHT_RATIO_OF_WATER_VAPOR_TO_DRY_AIR = 0.622

    col_wv_mixing_ratio_kg_kg = units.Quantity(point_data.QVAPOR.data.flatten(), 'kg/kg')
    # col_vapor_pressure = col_press_pa * col_wv_mixing_ratio_kg_kg / (MOLECULAR_WEIGHT_RATIO_OF_WATER_VAPOR_TO_DRY_AIR + col_water_vapor_mixing_ratio_kg_kg)
    # print(col_vapor_pressure)

    col_rh = mpcalc.relative_humidity_from_mixing_ratio(col_press_hpa,
                                                        col_temp_c,
                                                        col_wv_mixing_ratio_kg_kg)

    col_rh_pct = units.Quantity(col_rh, '%')

    col_u_ms = units.Quantity(point_data.U.data.flatten(), 'm/s')
    col_v_ms = units.Quantity(point_data.V.data.flatten(), 'm/s')
    col_u_kts = col_u_ms.to(ureg('kts'))
    col_v_kts = col_v_ms.to(ureg('kts'))

    # col_perturb_geopot = units.Quantity(point_col.PH.data, 'm**2/s**2')
    # col_base_geopot = units.Quantity(point_col.PHB.data, 'm**2/s**2')

    col_geopot = units.Quantity(point_data.PHP.data.flatten(), 'm**2/s**2')
    col_layer_height_msl = mpcalc.geopotential_to_height(col_geopot)

    sfc_press_pa = units.Quantity(point_data.PSFC.data.flatten(), 'pascal')
    sfc_press_hpa = sfc_press_pa.to(ureg('hectopascal'))

    sfc_temp_k = units.Quantity(point_data.T2.data.flatten(), 'K')
    sfc_temp_c = sfc_temp_k.to(ureg('degC'))
   
    sfc_wv_mixing_ratio_kg_kg = units.Quantity(point_data.Q2.data.flatten(), 'kg/kg')
    sfc_rh = mpcalc.relative_humidity_from_mixing_ratio(sfc_press_hpa,
                                                        sfc_temp_c,
                                                        sfc_wv_mixing_ratio_kg_kg)
    sfc_rh_pct = units.Quantity(sfc_rh, '%')

    sfc_u_ms = units.Quantity(point_data.U10.data.flatten(), 'm/s')
    sfc_v_ms = units.Quantity(point_data.V10.data.flatten(), 'm/s')
    sfc_u_kts = sfc_u_ms.to(ureg('kts'))
    sfc_v_kts = sfc_v_ms.to(ureg('kts'))

    sfc_hgt = units.Quantity(point_data.HGT.data.flatten(), 'm')

    surface_data = pd.DataFrame(data=[sfc_press_hpa.magnitude,
                                      sfc_temp_c.magnitude,
                                      sfc_rh_pct.magnitude,
                                      sfc_u_kts.magnitude,
                                      sfc_v_kts.magnitude,
                                      sfc_hgt.magnitude])

    column_data = pd.DataFrame(data=[col_press_hpa.magnitude,
                                     col_temp_c.magnitude,
                                     col_rh_pct.magnitude,
                                     col_u_kts.magnitude,
                                     col_v_kts.magnitude,
                                     col_layer_height_msl.magnitude])

    #remove data that is below the surface level
    column_data = column_data.T
    column_data = column_data[column_data[0] < sfc_press_hpa.magnitude[0]]

    sounding_data = pd.concat([surface_data.T, column_data],
                              axis=0,
                              ignore_index=True)

    sounding_data = sounding_data.round(4)
    print(sounding_data)

    if sel_pt_lat >= 0:
        lat_hemisphere = 'N'
    else:
        lat_hemisphere = 'S'

    if sel_pt_lon >= 0:
        lon_hemisphere = 'W'
    else:
        lon_hemisphere = 'E'

    internal_name = f'Pt: {round(lat, 4)} {round(lon, 4)} Gd: {round(sel_pt_lat, 4)} {round(sel_pt_lon, 4)} Md: WRF (No. {point_index})'

    header = {0:['RAOB/CSV','DTG','LAT','LON','ELEV','MOISTURE','WIND','GPM','MISSING','RAOB/DATA','PRES'],
              1:[internal_name,date_str,np.abs(sel_pt_lat),np.abs(sel_pt_lon),sfc_hgt.magnitude[0],'RH','kts','MSL',-999,'','TEMP'],
              2:['','',lat_hemisphere,lon_hemisphere,'m','','U/V','','','','RH'],
              3:['','','','','','','','','','','UU'],
              4:['','','','','','','','','','','VV'],
              5:['','','','','','','','','','','GPM']}

    header_data = pd.DataFrame(data=header)

    final_data = pd.concat([header_data, sounding_data],
                           axis=0,
                           ignore_index=True)

    file_name = f'WRF_RAOB_pt{point_index}_{datetime.strftime(point_time, "%Y%m%d_%H%M%S")}_frm{datetime.strftime(time_utc, "%Y%m%d_%H%M%S")}.csv'

    if not os.path.exists(save_dir):
        os.makedirs(save_dir)

    final_data.to_csv(os.path.join(save_dir, file_name), index=False, header=False)

    if previous_points is not None:
        return (sel_pt_lat, sel_pt_lon, point_index), previous_points
    else:
        return (sel_pt_lat, sel_pt_lon, point_index)

if __name__ == "__main__":

    parser = argparse.ArgumentParser(description=f'Create a meteogram data file from model data')
    parser.add_argument('-pt', '--point-sounding', 
                        help='specify for the script to create a single point sounding, which takes a point and time as input and produces a single sounding', 
                        nargs=7, 
                        metavar=('lat', 'lon', 'year', 'month', 'day', 'hour', 'minute'),
                        dest='point_sounding',
                        default=None)
    parser.add_argument('-cs', '--cross-section', 
                        help='specify for the script to create soundings for a cross-section, which takes input points and times and takes soundings corresponding to (the best match for) each of those points',
                        nargs=1, 
                        type=str,
                        metavar=('point_file'),
                        dest='cross_section',
                        default=None)
    parser.add_argument('-th', '--time-height', 
                        help='specify for the script to create soundings for a time-height plot, by iterating through all the files and taking a sounding from each at the same point',
                        nargs=2, 
                        type=float,
                        metavar=('lat', 'lon'),
                        dest='time_height',
                        default=None)
    parser.add_argument('-d', '--skip-duplicates',
                        help='if you would like duplicate points to be skipped (only valid for cross-section)',
                        action='store_true',
                        dest='skip_duplicate_pts')
    parser.add_argument('-vp', '--verification-plot',
                        help='generate a plot to verify location of sounding(s)',
                        dest='ver_plot',
                        action='store_true')
    parser.add_argument('-gp', '--grid-points',
                        help='save a csv to show location of sounding grid point(s)',
                        dest='save_gp',
                        action='store_true')
    parser.add_argument('--wrf-input-freq', 
                        help='specify the frequency of WRF input files, as a pandas frequency string. will produce an undefined result if frequency does not match input files',
                        nargs=1, 
                        type=str,
                        metavar=('freq_str'),
                        dest='wrf_freq',
                        default=None)
    parser.add_argument('--force-time', 
                        help='specify for the script to create soundings for a cross-section, which takes input points and times and takes soundings corresponding to (the best match for) each of those points', 
                        type=str,
                        metavar=('point_file'),
                        dest='forced_time',
                        default=None)
    parser.add_argument('model',  
                        help='specify data from what model is to be used. will error or produce undefined result if data format differs from specifed model',
                        choices=['hrrr', 'gfs', 'era5', 'wrf'],
                        type=str)
    parser.add_argument('input_file_directory',  
                        help='directory to read input files from',
                        type=str)
    parser.add_argument('save_directory',  
                        help='directory to save sounding(s) to',
                        type=str)

    args = parser.parse_args()
    print(args)

    ## Quality Checks/Sanitizing of input arguments

    if args.time_height or args.point_sounding:
        if args.time_height:
            sounding_lat = float(args.time_height[0])
            sounding_lon = float(args.time_height[1])
        if args.point_sounding:
            sounding_lat = float(args.point_sounding[0])
            sounding_lon = float(args.point_sounding[1])

        if (sounding_lat > 90) or (sounding_lat < -90):
            raise ValueError("Invalid latitude for sounding point, must +90 =< x =< -90")
        if (sounding_lon > 180) or (sounding_lon < -180):
            raise ValueError("Invalid longitude for sounding point, must +180 =< x =< -180")
    elif args.cross_section:
        cs_points = []
        try:
            with open(args.cross_section[0], newline='', mode='r', encoding='utf-8-sig') as cspointcsv:
                reader = csv.reader(cspointcsv, delimiter=',')
                for row in reader:
                    cs_point = (float(row[0]), float(row[1]), datetime.strptime(str(row[2]), "%Y-%m-%d %H:%M:%S"))
                    cs_points += [cs_point]
        except Exception as err:
            raise err
    elif (args.time_height and args.cross_section) or (args.time_height and args.point_sounding) or (args.cross_section and args.point_sounding):
        raise ValueError("Please only select one sounding creation methodology.")
    else:
        raise ValueError("Please select a sounding methodology.")

    if args.input_file_directory[-1:] != "/":
        input_dir = args.input_file_directory + "/"
    else:
        input_dir = args.input_file_directory

    print(f"Reading files from '{input_dir}'")

    if args.save_directory[-1:] != "/":
        save_dir = args.save_directory + "/"
    else:
        save_dir = args.save_directory

    print(f"Saving sounding(s) to '{save_dir}'")

    if args.model == 'era5':
        input_files = sorted(glob.glob(f'{input_dir}*UA_ERA5.grib'))
    elif args.model == 'hrrr':
        input_files = sorted(glob.glob(f'{input_dir}*.grib2'))
    elif args.model == 'wrf':
        input_files = sorted(glob.glob(f'{input_dir}wrfout*'))

    if not input_files:
        raise FileNotFoundError("No input files found in specified directory.")
    else:
        num_input_files = len(input_files)
        print(f"Successfully found {num_input_files} files.\nIndexing valid times...")
        file_times = create_file_times_dict(input_files, args.model)
        print("Done.")

    unit_registry = pint.UnitRegistry()

    if args.ver_plot:
        print("Verification plot requested.")

    if args.wrf_freq:
        wrf_freq = args.wrf_freq
    else:
        wrf_freq = DEFAULT_WRF_INPUT_FREQUENCY

    if args.forced_time:
        forced_time = args.forced_time
    else:
        forced_time = None

    grid_pts = []
    point_idx = 0

    if args.cross_section:
        with tqdm(miniters=0, total=len(cs_points), desc=f'Making RAOB cross-section from {args.model.upper()} data...', ascii=" ✎✏✐█") as progress:
            if args.skip_duplicate_pts:
                dup_pts = []
                if args.model == "hrrr":
                    for cs_point in cs_points:
                        grid_pt, dup_pts = parse_hrrr_data(parse_file_times(file_times, cs_point[2], '1h'),
                                                            save_dir,
                                                            cs_point[0],
                                                            cs_point[1],
                                                            cs_point[2],
                                                            point_idx,
                                                            dup_pts)
                        if grid_pt:
                            grid_pts.append(grid_pt)
                        progress.update()
                        point_idx += 1
                elif args.model == "era5":
                    for cs_point in cs_points:
                        file_time = parse_file_times(file_times, forced_time, wrf_freq) if forced_time else parse_file_times(file_times, cs_point[2], wrf_freq)
                        grid_pt, dup_pts = parse_era5_data(file_time,
                                                          save_dir,
                                                          cs_point[0],
                                                          cs_point[1],
                                                          cs_point[2],
                                                          point_idx,
                                                          dup_pts)
                        if grid_pt:
                            grid_pts.append(grid_pt)
                        progress.update()
                        point_idx += 1
                elif args.model == "wrf":
                    for cs_point in cs_points:
                        grid_pt, dup_pts = parse_wrf_data(parse_file_times(file_times, cs_point[2], wrf_freq),
                                                          save_dir,
                                                          cs_point[0],
                                                          cs_point[1],
                                                          cs_point[2],
                                                          point_idx,
                                                          dup_pts)
                        if grid_pt:
                            grid_pts.append(grid_pt)
                        progress.update()
                        point_idx += 1
            else:
                if args.model == "hrrr":
                    for cs_point in cs_points:
                        grid_pt = parse_hrrr_data(parse_file_times(file_times, cs_point[2], '1h'),
                                                   save_dir,
                                                   cs_point[0],
                                                   cs_point[1],
                                                   cs_point[2],
                                                   point_idx)
                        if grid_pt:
                            grid_pts.append(grid_pt)
                        progress.update()
                        point_idx += 1
                elif args.model == "wrf":
                    for cs_point in cs_points:
                        grid_pt = parse_wrf_data(parse_file_times(file_times, cs_point[2], wrf_freq),
                                                 save_dir,
                                                 cs_point[0],
                                                 cs_point[1],
                                                 cs_point[2],
                                                 point_idx)
                        if grid_pt:
                            grid_pts.append(grid_pt)
                        progress.update()
                        point_idx += 1


    if args.time_height:
        with tqdm(miniters=0, total=num_input_files, desc=f'Making RAOB time-height from {args.model.upper()} data...', ascii=" ✎✏✐█") as progress:
            if args.model == "hrrr":
                for file in input_files:
                    grid_pt = parse_hrrr_data(file,
                                              save_dir,
                                              sounding_lat,
                                              sounding_lon,
                                              get_file_time(file, 'hrrr'))
                    grid_pts.append(grid_pt)
                    progress.update()
            elif args.model == "era5":
                for file in input_files:
                    grid_pt = parse_era5_data(file,
                                              save_dir,
                                              sounding_lat,
                                              sounding_lon,
                                              get_file_time(file, 'era5'))
                    grid_pts.append(grid_pt)
                    progress.update()
            elif args.model == "wrf":
                for file in input_files:
                    grid_pt = parse_wrf_data(file,
                                             save_dir,
                                             sounding_lat,
                                             sounding_lon,
                                             get_file_time(file, 'wrf'))
                    grid_pts.append(grid_pt)
                    progress.update()

    if args.point_sounding:
        print(f"Making RAOB point sounding from {args.model.upper()} data...")
        pt_time = datetime(int(args.point_sounding[2]), int(args.point_sounding[3]), int(args.point_sounding[4]), int(args.point_sounding[5]), int(args.point_sounding[6]))
        if args.model == "hrrr":
            grid_pt = parse_hrrr_data(parse_file_times(file_times, pt_time, '1h'),
                                      save_dir,
                                      sounding_lat,
                                      sounding_lon,
                                      pt_time)
            grid_pts.append(grid_pt)
        elif args.model == "era5":
            grid_pt = parse_era5_data(parse_file_times(file_times, pt_time, '1h'),
                                      save_dir,
                                      sounding_lat,
                                      sounding_lon,
                                      pt_time)
            grid_pts.append(grid_pt)
        elif args.model == "wrf":
            grid_pt = parse_wrf_data(parse_file_times(file_times, pt_time, wrf_freq),
                                     save_dir,
                                     sounding_lat,
                                     sounding_lon,
                                     pt_time)
            grid_pts.append(grid_pt)


    if args.ver_plot:
        if args.time_height:
            lats = args.time_height[0]
            lons = args.time_height[1]
        elif args.point_sounding:
            lats = args.point_sounding[0]
            lons = args.point_sounding[1]
        elif args.cross_section:
            lats = []
            lons = []
            for cs_point in cs_points:
                lats.append(cs_point[0])
                lons.append(cs_point[1])

        generate_verification_plot(lats,
                                   lons,
                                   save_dir)

        print("Verification plot generated.")

    if args.save_gp:
        grid_df = pd.DataFrame(data=grid_pts, columns=['lat', 'lon', 'num'])
        grid_df.to_csv(os.path.join(save_dir, f"RAOB_ModelSoundings_Gridpoints_{args.model.upper()}.csv"), index=None)
        print("Used grid points saved to file.")

    remove_files(input_dir)

    print("Done!")






