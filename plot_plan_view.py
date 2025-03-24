'''
IMPLEMENTED VARIABLES AND DISPLAY TYPES:

    Shaded/Filled contours: 
    Cloud Cover ['cldcvr_cf'], Relative Humidity (%) ['rh_cf'], Wind Speed (kts) ['wind_cf'],
    Vertical Velocity (ft/min) ['vv_cf'], Temperature (°C) ['temp_cf']

    Contour Lines:
    Temperature (°C) ['temp_c'], Geopotential Height (m) ['gpm_c']

    Vectors:
    Wind (kts) ['wind_vectors']

    Barbs:
    Wind (kts) ['wind_barbs']

    Grid Point Values:
    Wind (kts) ['wind_values']
    
'''

## Turn on N-point variable smoothing. Only applies to shaded/filled contour variables.
## (Default: False)
VARIABLE_SMOOTHING = False

## Amount of points to use for N-point variable smoothing. Must be 5 or 9.
## (Default: 9)
SMOOTHING_POINTS = 9

## Amount of passes to use for N-point variable smoothing.
## (Default: 1)
SMOOTHING_PASSES = 1

DEFAULT_LEVELS_TO_PLOT = [9999]#[1000, 950, 900, 850, 700, 500, 300, 250]

## Sort output plots into indivdual folders by level
## (Default: True)
SORT_PLOTS_BY_LEVEL = True

## Set the DPI of the saved output file
## (Default: 300)
FILE_DPI = 200












'''


            STATIC CODE BELOW HERE


'''





import os
import glob
from datetime import datetime, timedelta
import argparse
import csv
import json

import pandas as pd
import xarray as xr
import numpy as np
import matplotlib as mpl
from matplotlib import cm
import matplotlib.pyplot as plt
import matplotlib.patheffects as PathEffects
import cartopy
import cartopy.crs as crs
import cartopy.feature as cfeature
import cartopy.io.shapereader as shpreader
from cartopy.feature import NaturalEarthFeature
import metpy
from metpy.calc import (smooth_n_point, wind_speed, 
                        geopotential_to_height, vertical_velocity)
import metpy.calc as mpcalc
from metpy.units import units
import netCDF4
from tqdm.auto import tqdm
from adjustText import adjust_text 
import pint

from __internal_funcs import (plot_towns, draw_logo, 
                              plot_points, define_hi_res_fig,
                              define_gearth_compat_fig, save_figs_to_kml)

HRRR_VARIABLE_TABLE = {
    #internal id: (grib id, base unit, data file (if diff from pressure levels))
    'temp': ('t', 'K'),
    'dpt': ('dpt', 'K'),
    'rh': ('r', '%'),
    'gpm': ('gh', 'm^2/s^2'),
    'vv': ('w', 'pascal/s'),
    'cldcvr': ('cc', '%'),
    'u': ('u', 'm/s'),
    'v': ('v', 'm/s'),
    'mxgrat' : ('q', 'kg/kg'),
    'pressure': ('isobaricInhPa', 'hectopascal'),
    'SFC_t': ('t2m', 'K', 'm2_data'),
    'SFC_dpt': ('d2m', 'K', 'm2_data'),
    'SFC_rh': ('r2', '%', 'm2_data'),
    'SFC_landt': ('t', 'K', 'sfc_data'),
    'SFC_u': ('u10', 'm/s', 'm10_data'),
    'SFC_v': ('v10', 'm/s', 'm10_data')
}

ERA5_VARIABLE_TABLE = {
    'temp': 't',
    'dpt': 'dpt',
    'rh': 'r',
    'gpm': 'z',
    'vv': 'w',
    'cldcvr': 'cc',
    'u': 'u',
    'v': 'v',
    'mxgrat' : 'q',
    'pressure': 'isobaricInhPa',
}

CONTOURF_ZORDER = 1
CONTOUR_ZORDER = 3
SHAPE_ZORDER = 6

np.seterr(divide = 'ignore', invalid='ignore')
mpl.use('agg')

def plot_plan_view_hrrr(file_path: str,
                        save_dir: str,
                        sel_level: int,
                        products: list[str, ...],
                        points: list[tuple[float, float, str], ...],
                        bbox: list[float, float, float, float],
                        unit_reg: pint.UnitRegistry,
                        **kwargs) -> None:

    if level == 9999: #Surface flag
        sfc_flag = True
    else:
        sfc_flag = False

    if sfc_flag:
        sfc_data = xr.open_dataset(file_path, engine="cfgrib", filter_by_keys={'typeOfLevel': 'surface', 'stepType': 'instant'})
        m2_data = xr.open_dataset(file_path, engine="cfgrib", filter_by_keys={'typeOfLevel': 'heightAboveGround', 'level': 2})
        m10_data = xr.open_dataset(file_path, engine="cfgrib", filter_by_keys={'typeOfLevel': 'heightAboveGround', 'level': 10})

        file_time_str = str(sfc_data['time'].values)
    else:
        col_data = xr.open_dataset(file_path, engine="cfgrib", filter_by_keys={'typeOfLevel': 'isobaricInhPa'})
        data = col_data.sel(isobaricInhPa=level)
        file_time_str = str(data['time'].values)

    
    file_day = file_time_str[:10].replace("-", "_")
    file_time = file_time_str[11:16].replace(":", "") + "UTC"

    if kwargs.get('save_to_kmz'):  
        fig, ax, cbfig, cbax = define_gearth_compat_fig((bbox[2], bbox[3]),
                                                        (bbox[0], bbox[1]))

    else:
        fig, ax = define_hi_res_fig((bbox[2], bbox[3]),
                                    (bbox[0], bbox[1]))

        cbfig = None
        cbax = None

    prodstr = ""
    filestr = ""

    if sfc_flag:
        prodstr, filestr = draw_surface_products(fig, ax,
                                                 products,
                                                 'hrrr',
                                                 prodstr,
                                                 filestr,
                                                 unit_reg,
                                                 fig_for_cb=cbfig,
                                                 ax_for_cb=cbax,
                                                 sfc_data=sfc_data,
                                                 m2_data=m2_data,
                                                 m10_data=m10_data,
                                                 **kwargs)

        prodstr, filestr = draw_surface_wind_display(fig, ax,
                                                     m10_data,
                                                     products,
                                                     'hrrr',
                                                     prodstr,
                                                     filestr,
                                                     unit_reg,
                                                     **kwargs)


    else:
        prodstr, filestr = draw_contourf_lines(fig, ax,
                                               data,
                                               sel_level,
                                               products,
                                               'hrrr',
                                               prodstr,
                                               filestr,
                                               unit_reg,
                                               fig_for_cb=cbfig,
                                               ax_for_cb=cbax,
                                               **kwargs)

        prodstr, filestr = draw_contour_lines(fig, ax,
                                              data,
                                              sel_level,
                                              products,
                                              'hrrr',
                                              prodstr,
                                              filestr,
                                              unit_reg,
                                              **kwargs)

        prodstr, filestr = draw_wind_display(fig, ax,
                                             data,
                                             sel_level,
                                             products,
                                             'hrrr',
                                             prodstr,
                                             filestr,
                                             unit_reg,
                                             **kwargs)

    if points:
        plot_points(plt, ax,
                    points,
                    **kwargs)

    prod_titles = prodstr.split('#')
    if len(prod_titles[1:]) > 1:
        prod_titles[-1] = "and " + prod_titles[-1]
        if len(prodstr) >= 40:
            prod_titles[round(len(prod_titles)/2)] = '\n' + prod_titles[round(len(prod_titles)/2)]
    
    titlestr = ", ".join(prod_titles)
    titlestr = titlestr[2:]

    if SORT_PLOTS_BY_LEVEL:
        save_dir = os.path.join(save_dir, f"{sel_level} Level")

    if not os.path.exists(save_dir):
        os.makedirs(save_dir)

    if sfc_flag:
        file_name = "PlanView_SFC" + filestr + "_" + file_day + "_" + file_time +"_HRRR"
        descstr = f"Surface Level\nValid at {np.datetime_as_string(sfc_data['valid_time'].values, timezone='UTC')[:-11].replace('T', ' ')} UTC"
    else:
        file_name = "PlanView_" + str(sel_level) + filestr + "_" + file_day + "_" + file_time +"_HRRR"
        descstr = f"{str(data['isobaricInhPa'].values)} hPa Level\nValid at {np.datetime_as_string(data['valid_time'].values, timezone='UTC')[:-11].replace('T', ' ')} UTC"

    if kwargs.get('save_to_kmz'):
        
        layer_name = f"HRRR Reanalysis at {np.datetime_as_string(data['valid_time'].values, timezone='UTC')[:-11].replace('T', ' ')} UTC"
        layer_desc = titlestr.replace('/\n', ' ') + f', at {sel_level} hPa level'

        dest_path = os.path.join(save_dir, file_name + ".kmz")

        save_figs_to_kml(dest_path,
                         (bbox[2], bbox[3]),
                         (bbox[0], bbox[1]),
                         fig,
                         [layer_name],
                         [layer_desc],
                         colorbar_fig=cbfig)

        plt.close(fig)
        if cbfig:
            plt.close(cbfig)

    else:
        plot_towns(ax, 
                   (bbox[2], bbox[3]),
                   (bbox[0], bbox[1]), 
                   scale_rank=kwargs.pop('plot_towns_scale_rank', 5))

        draw_logo(ax)
        
        ax.set_title(f'HRRR Reanalysis {titlestr}', loc='left', fontweight='bold', fontsize=15)
        plt.title(descstr, loc='right')

        dest_path = os.path.join(save_dir, file_name + ".png")

        fig.savefig(dest_path, bbox_inches="tight", dpi=FILE_DPI)
        plt.close(fig)

def plot_plan_view_era5(file_path: str,
                        save_dir: str,
                        sel_level: int,
                        products: list[str, ...],
                        points: list[tuple[float, float, str], ...],
                        bbox: list[float, float, float, float],
                        **kwargs) -> None:

    if level == "9999": #Surface flag
        sfc_flag = True
    else:
        sfc_flag = False

    if sfc_flag:
        sfc_data = xr.open_dataset(file_path, engine="cfgrib", filter_by_keys={'typeOfLevel': 'surface', 'stepType': 'instant'})
        m2_data = xr.open_dataset(file_path, engine="cfgrib", filter_by_keys={'typeOfLevel': 'heightAboveGround', 'level': 2})
        m10_data = xr.open_dataset(file_path, engine="cfgrib", filter_by_keys={'typeOfLevel': 'heightAboveGround', 'level': 10})
    else:
        col_data = xr.open_dataset(file_path, engine="cfgrib")
        data = col_data.sel(isobaricInhPa=level)

    file_time_str = str(data['time'].values)
    file_day = file_time_str[:10].replace("-", "_")
    file_time = file_time_str[11:16].replace(":", "") + "UTC"

    if kwargs.get('save_to_kmz'):  
        fig, ax, cbfig, cbax = define_gearth_compat_fig((bbox[2], bbox[3]),
                                                        (bbox[0], bbox[1]))

    else:
        fig, ax = define_hi_res_fig((bbox[2], bbox[3]),
                                    (bbox[0], bbox[1]))

        cbfig = None
        cbax = None

    prodstr = ""
    filestr = ""

    if sfc_flag:
        raise NotImplementedError("Not implemented!!!1!")
    else:
        prodstr, filestr = draw_contourf_lines(fig, ax,
                                               data,
                                               sel_level,
                                               products,
                                               'era5',
                                               prodstr,
                                               filestr,
                                               fig_for_cb=cbfig,
                                               ax_for_cb=cbax,
                                               **kwargs)

        prodstr, filestr = draw_contour_lines(fig, ax,
                                              data,
                                              sel_level,
                                              products,
                                              'era5',
                                              prodstr,
                                              filestr,
                                              **kwargs)

        prodstr, filestr = draw_wind_display(fig, ax,
                                             data,
                                             sel_level,
                                             products,
                                             'era5',
                                             prodstr,
                                             filestr,
                                             **kwargs)

    if points:
        plot_points(plt, ax,
                    points,
                    **kwargs)
    
    prod_titles = prodstr.split('#')
    if len(prod_titles[1:]) > 1:
        prod_titles[-1] = "and " + prod_titles[-1]
        if len(prodstr) >= 40:
            prod_titles[round(len(prod_titles)/2)] = '\n' + prod_titles[round(len(prod_titles)/2)]
    
    titlestr = ", ".join(prod_titles)
    titlestr = titlestr[2:]

    if kwargs.get('save_to_kmz'):
        file_name = "PlanView_" + str(sel_level) + filestr + "_" + file_day + "_" + file_time +"_ERA5"
        layer_name = f"ERA5 Reanalysis at {np.datetime_as_string(data['valid_time'].values, timezone='UTC')[:-11].replace('T', ' ')} UTC"
        layer_desc = titlestr.replace('/\n', ' ') + f', at {sel_level} hPa level'
        
        if not os.path.exists(save_dir):
            os.makedirs(save_dir)

        dest_path = os.path.join(save_dir, file_name + ".kmz")

        save_figs_to_kml(dest_path,
                         (bbox[2], bbox[3]),
                         (bbox[0], bbox[1]),
                         fig,
                         [layer_name],
                         [layer_desc],
                         colorbar_fig=cbfig)

        plt.close(fig)
        if cbfig:
            plt.close(cbfig)

    else:
        plot_towns(ax, 
               (bbox[2], bbox[3]),
               (bbox[0], bbox[1]), 
               scale_rank=kwargs.pop('plot_towns_scale_rank', 5))

        draw_logo(ax)
    
        descstr = f"{str(data['isobaricInhPa'].values)} hPa Level\nValid at {np.datetime_as_string(data['valid_time'].values, timezone='UTC')[:-11].replace('T', ' ')} UTC"
        ax.set_title(f'ERA5 Reanalysis {titlestr}', loc='left', fontweight='bold', fontsize=15)
        ax.set_title(descstr, loc='right')

        file_name = "PlanView_" + str(sel_level) + filestr + "_" + file_day + "_" + file_time +"_ERA5"
        
        if not os.path.exists(save_dir):
            os.makedirs(save_dir)

        dest_path = os.path.join(save_dir, file_name + ".png")

        fig.savefig(dest_path, bbox_inches="tight", dpi=FILE_DPI)
        plt.close(fig)

def plot_plan_view_wrf(file_path: str,
                        save_dir: str,
                        sel_level: int,
                        products: list[str, ...],
                        points: list[tuple[float, float, str], ...],
                        bbox: list[float, float, float, float],
                        **kwargs) -> None:

    if level == "9999": #Surface flag
        sfc_flag = True
    else:
        sfc_flag = False

    if sfc_flag:
        raise NotImplementedError("Surface level plotting for WRF not implemented")

    else:
        col_data = xr.open_dataset(file_path, engine="cfgrib", filter_by_keys={'typeOfLevel': 'isobaricInhPa'})
        data = col_data.sel(isobaricInhPa=level)

    file_time_str = str(data['time'].values)
    file_day = file_time_str[:10].replace("-", "_")
    file_time = file_time_str[11:16].replace(":", "") + "UTC"

    if kwargs.get('save_to_kmz'):  
        fig, ax, cbfig, cbax = define_gearth_compat_fig((bbox[2], bbox[3]),
                                                        (bbox[0], bbox[1]))

    else:
        fig, ax = define_hi_res_fig((bbox[2], bbox[3]),
                                    (bbox[0], bbox[1]))

        cbfig = None
        cbax = None

    prodstr = ""
    filestr = ""

    prodstr, filestr = draw_contourf_lines(fig, ax,
                                           data,
                                           sel_level,
                                           products,
                                           'hrrr',
                                           prodstr,
                                           filestr,
                                           fig_for_cb=cbfig,
                                           ax_for_cb=cbax,
                                           **kwargs)

    prodstr, filestr = draw_contour_lines(fig, ax,
                                          data,
                                          sel_level,
                                          products,
                                          'hrrr',
                                          prodstr,
                                          filestr,
                                          **kwargs)

    prodstr, filestr = draw_wind_display(fig, ax,
                                         data,
                                         sel_level,
                                         products,
                                         'hrrr',
                                         prodstr,
                                         filestr,
                                         **kwargs)

    if points:
        plot_points(plt, ax,
                    points,
                    **kwargs)

    prod_titles = prodstr.split('#')
    if len(prod_titles[1:]) > 1:
        prod_titles[-1] = "and " + prod_titles[-1]
        if len(prodstr) >= 40:
            prod_titles[round(len(prod_titles)/2)] = '\n' + prod_titles[round(len(prod_titles)/2)]
    
    titlestr = ", ".join(prod_titles)
    titlestr = titlestr[2:]

    if kwargs.get('save_to_kmz'):
        file_name = "PlanView_" + str(sel_level) + filestr + "_" + file_day + "_" + file_time +"_HRRR"
        layer_name = f"HRRR Reanalysis at {np.datetime_as_string(data['valid_time'].values, timezone='UTC')[:-11].replace('T', ' ')} UTC"
        layer_desc = titlestr.replace('/\n', ' ') + f', at {sel_level} hPa level'
        
        if not os.path.exists(save_dir):
            os.makedirs(save_dir)

        dest_path = os.path.join(save_dir, file_name + ".kmz")

        save_figs_to_kml(dest_path,
                         (bbox[2], bbox[3]),
                         (bbox[0], bbox[1]),
                         fig,
                         [layer_name],
                         [layer_desc],
                         colorbar_fig=cbfig)

        plt.close(fig)
        if cbfig:
            plt.close(cbfig)

    else:
        plot_towns(ax, 
                   (bbox[2], bbox[3]),
                   (bbox[0], bbox[1]), 
                   scale_rank=kwargs.pop('plot_towns_scale_rank', 5))

        draw_logo(ax)
        
        descstr = f"{str(data['isobaricInhPa'].values)} hPa Level\nValid at {np.datetime_as_string(data['valid_time'].values, timezone='UTC')[:-11].replace('T', ' ')} UTC"
        ax.set_title(f'HRRR Reanalysis {titlestr}', loc='left', fontweight='bold', fontsize=15)
        plt.title(descstr, loc='right')

        file_name = "PlanView_" + str(sel_level) + filestr + "_" + file_day + "_" + file_time +"_HRRR"
        
        if not os.path.exists(save_dir):
            os.makedirs(save_dir)

        dest_path = os.path.join(save_dir, file_name + ".png")

        fig.savefig(dest_path, bbox_inches="tight", dpi=FILE_DPI)
        plt.close(fig)



def draw_contourf_lines(fig: mpl.figure.Figure,
                        ax: mpl.axes.Axes,
                        data: xr.Dataset | netCDF4.Dataset,
                        sel_level: int,
                        products: list[str, ...],
                        model: str,
                        prodstr: str,
                        prodabbr: str,
                        unit_reg: pint.UnitRegistry,
                        fig_for_cb: mpl.figure.Figure = None,
                        ax_for_cb: mpl.axes.Axes = None,
                        **kwargs) -> tuple[str, str]:

    if fig_for_cb and ax_for_cb:
        cb_fig = fig_for_cb
        cb_ax = ax_for_cb

        cbar_opts = {'rotation': -90, 'color': 'k', 'labelpad': 20}
    else:
        cb_fig = fig
        cb_ax = None

        cbar_opts = {}


    if model == 'hrrr':
        vtable = HRRR_VARIABLE_TABLE
        lat = data.variables['latitude'][:]
        lon = data.variables['longitude'][:]
    elif model == 'era5':
        vtable = ERA5_VARIABLE_TABLE
        lat = data.variables['latitude'][:]
        lon = data.variables['longitude'][:]

    if "temp_cf" in products:
        levels = np.arange(kwargs.pop('temp_level_min', -30 if sel_level > 500 else -70),
                           kwargs.pop('temp_level_max', 45 if sel_level > 500 else 20),
                           kwargs.pop('temp_level_step', 2))

        var = data.variables[vtable['temp'][0]][:]
        var_baseunits = units.Quantity(var, vtable['temp'][1])

        default_unit = 'degC'

        if kwargs.get('temp_units'):
            unit = unit_reg.Unit(kwargs.pop('temp_units', default_unit))
            var_units = var_baseunits.to(unit_reg(unit))
        else:
            unit = unit_reg.Unit(default_unit)
            var_units = var_baseunits

        if VARIABLE_SMOOTHING == True:

            var_contourf = ax.contourf(lon, lat, 
                                       smooth_n_point(var_units.magnitude,
                                                      SMOOTHING_POINTS,
                                                      SMOOTHING_PASSES), 
                                       levels=levels, 
                                       transform=crs.PlateCarree(), 
                                       cmap=kwargs.pop('temp_pallete', cm.Spectral_r), 
                                       vmin=levels.min(),
                                       vmax=levels.max(),
                                       extend='both',
                                       zorder=kwargs.pop('cf_zorder', CONTOURF_ZORDER))
        else:
            var_contourf = ax.contourf(lon, lat, 
                                       var_units.magnitude, 
                                       levels=levels, 
                                       transform=crs.PlateCarree(), 
                                       cmap=kwargs.pop('temp_pallete', cm.Spectral_r), 
                                       vmin=levels.min(),
                                       vmax=levels.max(),
                                       extend='both',
                                       zorder=kwargs.pop('cf_zorder', CONTOURF_ZORDER))
        cb = cb_fig.colorbar(var_contourf, 
                          ax=ax, 
                          cax=cb_ax,
                          ticks=levels[::2], 
                          orientation=kwargs.pop('colorbar_orientation', 'vertical'), 
                          shrink=0.77)
        cb.set_label(f'Temperature ({unit:~C})', size='x-large', **cbar_opts)

        prodstr += f"#Temperature ({unit:~C}, shaded contours)"
        prodabbr += "_T"

    if "cldcvr_cf" in products:
        levels = np.arange(kwargs.pop('cldcvr_level_min', 0),
                           kwargs.pop('cldcvr_level_max', 1.1),
                           kwargs.pop('cldcvr_level_step', 0.1))

        var = data.variables[vtable['cldcvr'][0]][:]
        var_units = units.Quantity(var, vtable['cldcvr'][1])

        if VARIABLE_SMOOTHING == True:
            var_contourf = ax.contourf(lon, lat, 
                                       smooth_n_point(var_units.magnitude,
                                                      SMOOTHING_POINTS,
                                                      SMOOTHING_PASSES), 
                                       levels=levels, 
                                       transform=crs.PlateCarree(), 
                                       cmap=kwargs.pop('cldcvr_pallete', cm.Blues), 
                                       zorder=kwargs.pop('cf_zorder', CONTOURF_ZORDER))
        else:
            var_contourf = ax.contourf(lon, lat, 
                                       var_units.magnitude, 
                                       levels=levels, 
                                       transform=crs.PlateCarree(), 
                                       cmap=kwargs.pop('cldcvr_pallete', cm.Blues), 
                                       vmin=levels.min(),
                                       vmax=levels.max(),
                                       zorder=kwargs.pop('cf_zorder', CONTOURF_ZORDER))
        cb = cb_fig.colorbar(var_contourf, 
                          ax=ax, 
                          cax=cb_ax,
                          ticks=levels[::2], 
                          orientation=kwargs.pop('colorbar_orientation', 'vertical'), 
                          shrink=0.77)
        cb.set_label(f'Fraction of Cloud Cover ({unit:~C})', size='x-large', **cbar_opts)

        prodstr += f"#Fraction of Cloud Cover (({unit:~C}), shaded contours)"
        prodabbr += "_CC"

    if "rh_cf" in products:
        levels = np.arange(kwargs.pop('rh_level_min', 0),
                           kwargs.pop('rh_level_max', 1.1),
                           kwargs.pop('rh_level_step', 0.1))

        var = data.variables[vtable['rh'][0]][:]
        var_baseunits = units.Quantity(var, vtable['rh'][1])

        default_unit = '%'

        if kwargs.get('rh_units'):
            unit = unit_reg.Unit(kwargs.pop('rh_units', default_unit))
            var_units = var_baseunits.to(unit_reg(unit))
        else:
            unit = unit_reg.Unit(default_unit)
            var_units = var_baseunits


        if VARIABLE_SMOOTHING == True:
            var_contourf = ax.contourf(lon, lat, 
                                       smooth_n_point(var_units.magnitude, 
                                                      SMOOTHING_POINTS,
                                                      SMOOTHING_PASSES), 
                                       levels=levels, 
                                       transform=crs.PlateCarree(), 
                                       cmap=kwargs.pop('rh_pallete', cm.terrain_r), 
                                       extend='max',
                                       zorder=kwargs.pop('cf_zorder', CONTOURF_ZORDER))
        else:
            var_contourf = ax.contourf(lon, lat, 
                                       var_units.magnitude, 
                                       levels=levels, 
                                       transform=crs.PlateCarree(), 
                                       cmap=kwargs.pop('rh_pallete', cm.terrain_r), 
                                       vmin=levels.min(),
                                       vmax=levels.max(),
                                       extend='max',
                                       zorder=kwargs.pop('cf_zorder', CONTOURF_ZORDER))
        cb = cb_fig.colorbar(var_contourf, 
                          ax=ax, 
                          cax=cb_ax,
                          ticks=levels[::2], 
                          orientation=kwargs.pop('colorbar_orientation', 'vertical'), 
                          shrink=0.77,
                          extendrect=True)
        cb.set_label(f'Relative Humidty ({unit:~C})', size='x-large', **cbar_opts)

        prodstr += f"#Relative Humidity ({unit:~C}, shaded)"
        prodabbr += "_RH"

    if "wind_cf" in products:
        levels = np.arange(kwargs.pop('wind_level_min', 0),
                           kwargs.pop('wind_level_max', 40 if sel_level > 500 else 100),
                           kwargs.pop('wind_level_step', 2 if sel_level > 500 else 5))

        wind_speed = mpcalc.wind_speed(data[vtable['u'][0]], data[vtable['v'][0]])

        default_unit = 'knots'

        if kwargs.get('wind_units'):
            unit = unit_reg.Unit(kwargs.get('wind_units', default_unit))
            var_units = wind_speed.metpy.convert_units(unit)
        else:
            unit = unit_reg.Unit(default_unit)
            var_units = wind_speed.metpy.convert_units(unit)

        if VARIABLE_SMOOTHING == True:
            var_contourf = ax.contourf(lon, lat, 
                                       smooth_n_point(var_units, 
                                                      SMOOTHING_POINTS,
                                                      SMOOTHING_PASSES), 
                                       levels=levels, 
                                       transform=crs.PlateCarree(), 
                                       cmap=kwargs.pop('wind_pallete', cm.jet), 
                                       extend='max',
                                       zorder=kwargs.pop('cf_zorder', CONTOURF_ZORDER))
        else:
            var_contourf = ax.contourf(lon, lat, 
                                       var_units, 
                                       levels=levels, 
                                       transform=crs.PlateCarree(), 
                                       cmap=kwargs.pop('wind_pallete', cm.jet), 
                                       vmin=levels.min(),
                                       vmax=levels.max(),
                                       extend='max',
                                       zorder=kwargs.pop('cf_zorder', CONTOURF_ZORDER))
        cb = cb_fig.colorbar(var_contourf, 
                          ax=ax, 
                          cax=cb_ax,
                          ticks=levels[::2], 
                          orientation=kwargs.pop('colorbar_orientation', 'vertical'), 
                          shrink=0.77)
        cb.set_label(f'Wind Speed ({unit:~C})', size='x-large', **cbar_opts)

        prodstr += f"#Wind Speed ({unit:~C}, shaded)"
        prodabbr += "_WS"

    if "wd_cf" in products:
        levels = np.arange(kwargs.pop('wd_level_min', 0),
                           kwargs.pop('wd_level_max', 365),
                           kwargs.pop('wd_level_step', 15))

        wind_direction = mpcalc.wind_direction(data[vtable['u'][0]], data[vtable['v'][0]])

        default_unit = 'degrees'
        
        if kwargs.get('wd_units'):
            unit = unit_reg.Unit(kwargs.pop('wd_units', default_unit))
            var_units = wind_direction.metpy.convert_units(unit)
        else:
            unit = unit_reg.Unit(default_unit)
            var_units = wind_direction.metpy.convert_units(unit)

        if VARIABLE_SMOOTHING == True:
            var_contourf = ax.contourf(lon, lat, 
                                       smooth_n_point(var_units, 
                                                      SMOOTHING_POINTS,
                                                      SMOOTHING_PASSES), 
                                       levels=levels, 
                                       transform=crs.PlateCarree(), 
                                       cmap=kwargs.pop('wd_pallete', cm.twilight), 
                                       zorder=kwargs.pop('cf_zorder', CONTOURF_ZORDER))
        else:
            var_contourf = ax.contourf(lon, lat, 
                                       var_units, 
                                       levels=levels, 
                                       transform=crs.PlateCarree(), 
                                       cmap=kwargs.pop('wd_pallete', cm.twilight), 
                                       vmin=levels.min(),
                                       vmax=levels.max(),
                                       zorder=kwargs.pop('cf_zorder', CONTOURF_ZORDER))
        cb = cb_fig.colorbar(var_contourf, 
                          ax=ax, 
                          cax=cb_ax,
                          ticks=levels[::2], 
                          orientation=kwargs.pop('colorbar_orientation', 'vertical'), 
                          shrink=0.77)
        cb.set_label(f'Wind Direction ({unit:~C})', size='x-large', **cbar_opts)

        prodstr += f"#Wind Direction ({unit:~C}, shaded)"
        prodabbr += "_WD"

    if "vv_cf" in products:
        levels = np.arange(kwargs.pop('vv_level_min', -450),
                           kwargs.pop('vv_level_max', 450),
                           kwargs.pop('vv_level_step', 10))

        if model == "hrrr":
            vv = vertical_velocity(data[vtable['vv'][0]],
                                   data[vtable['pressure'][0]],
                                   data[vtable['temp'][0]],
                                   data[vtable['mxgrat'][0]])

            vv = vv.values * 196.9 #m/s to ft/min

        if VARIABLE_SMOOTHING == True:
            var_contourf = ax.contourf(lon, lat, 
                                       smooth_n_point(vv, 
                                                      SMOOTHING_POINTS,
                                                      SMOOTHING_PASSES), 
                                       levels=levels, 
                                       transform=crs.PlateCarree(), 
                                       cmap=kwargs.pop('vv_pallete', cm.bwr), 
                                       extend='both',
                                       zorder=kwargs.pop('cf_zorder', CONTOURF_ZORDER))
        else:
            var_contourf = ax.contourf(lon, lat, 
                                       vv, 
                                       levels=levels, 
                                       transform=crs.PlateCarree(), 
                                       cmap=kwargs.pop('vv_pallete', cm.bwr), 
                                       vmin=levels.min(),
                                       vmax=levels.max(),
                                       extend='both',
                                       zorder=kwargs.pop('cf_zorder', CONTOURF_ZORDER))
        cb = cb_fig.colorbar(var_contourf, 
                          ax=ax, 
                          cax=cb_ax,
                          ticks=levels[::2], 
                          orientation=kwargs.pop('colorbar_orientation', 'vertical'), 
                          shrink=0.77)
        cb.set_label('Vertical Velocity (ft/s)', size='x-large', **cbar_opts)

        prodstr += "#Vertical Velocity (ft/s, shaded)"
        prodabbr += "_VV"

    return prodstr, prodabbr

def draw_contour_lines(fig: mpl.figure.Figure,
                       ax: mpl.axes.Axes,
                       data: xr.Dataset | netCDF4.Dataset,
                       sel_level: int,
                       products: list[str, ...],
                       model: str,
                       prodstr: str,
                       prodabbr: str,
                       unit_reg: pint.UnitRegistry,
                       **kwargs) -> tuple[str, str]:

    num_prods = 0

    if model == 'hrrr':
        vtable = HRRR_VARIABLE_TABLE
        lat = data.variables['latitude'][:]
        lon = data.variables['longitude'][:]
    elif model == 'era5':
        vtable = ERA5_VARIABLE_TABLE
        lat = data.variables['latitude'][:]
        lon = data.variables['longitude'][:]


    if "gpm_c" in products:
        num_prods += 1
        levels = np.arange(kwargs.pop('gpm_level_min', sel_level - 500),
                           kwargs.pop('gpm_level_max', sel_level + 500),
                           kwargs.pop('gpm_level_step', 50))

        gpm = data[vtable['gpm'][0]]

        contours = ax.contour(lon, lat, 
                               gpm, 
                               levels=levels, 
                               cmap=kwargs.pop('gpm_c_colors', "gray"), 
                               alpha=kwargs.pop('gpm_c_alpha', 1), 
                               transform=crs.PlateCarree(), 
                               zorder=kwargs.pop('gpm_c_zorder', CONTOUR_ZORDER + (0.5*num_prods)), 
                               linewidths=kwargs.pop('gpm_c_linewidth', 1.8))
        ax.clabel(contours, 
                  inline=1, 
                  fontsize=kwargs.pop('gpm_c_label_fontsize', 10), 
                  fmt="%i", 
                  zorder=kwargs.pop('gpm_c_label_zorder', 6))

        prodstr += "#Geopotential Height (m, contours)"
        prodabbr += "_GPM"

    if "temp_c" in products:
        num_prods += 1
        levels = np.arange(kwargs.pop('temp_level_min', -30 if sel_level > 500 else -70),
                           kwargs.pop('temp_level_max', 45 if sel_level > 500 else 20),
                           kwargs.pop('temp_level_step', 1))

        contours = ax.contour(lon, lat, 
                              data[vtable['temp']] - 273.15, 
                              levels=levels, 
                              cmap=kwargs.pop('temp_c_colors', cm.RdYlBu_r), 
                              alpha=kwargs.pop('temp_c_alpha', 0.85), 
                              transform=crs.PlateCarree(), 
                              zorder=kwargs.pop('temp_c_zorder', CONTOUR_ZORDER + (0.5*num_prods)), 
                              linestyles=kwargs.pop('temp_c_linestyle', 'dashed'),
                              linewidths=kwargs.pop('temp_c_linewidth', 1.2))
        ax.clabel(contours, 
                  inline=1, 
                  fontsize=kwargs.pop('temp_c_label_fontsize', 10), 
                  fmt="%i", 
                  zorder=kwargs.pop('temp_c_label_zorder', 6))

        prodstr += "#Temperature (°C, colored contours)"
        prodabbr += "_T"

    return prodstr, prodabbr

def draw_wind_display(fig: mpl.figure.Figure,
                      ax: mpl.axes.Axes,
                      data: xr.Dataset | netCDF4.Dataset,
                      sel_level: int,
                      products: list[str, ...],
                      model: str,
                      prodstr: str,
                      prodabbr: str,
                      unit_reg: pint.UnitRegistry,
                      **kwargs) -> tuple[str, str]:

    slicer_idx = kwargs.pop('wind_display_indexing', 4)

    if model == 'hrrr':
        vtable = HRRR_VARIABLE_TABLE
        lat = data.variables['latitude'][:]
        lon = data.variables['longitude'][:]
    elif model == 'era5':
        vtable = ERA5_VARIABLE_TABLE
        lat = data.variables['latitude'][:]
        lon = data.variables['longitude'][:]

    u = data[vtable['u'][0]]
    v = data[vtable['v'][0]]

    wind_speed = mpcalc.wind_speed(data[vtable['u'][0]], data[vtable['v'][0]])

    default_unit = 'knots'

    if kwargs.get('wind_display_units'):
        unit = unit_reg.Unit(kwargs.pop('wind_display_units', default_unit))

        u_units = u.metpy.convert_units(unit)
        v_units = v.metpy.convert_units(unit)

        wind_speed_units = wind_speed.metpy.convert_units(unit)

    else:
        unit = unit_reg.Unit(default_unit)

        u_units = u.metpy.convert_units(default_unit)

        u_units = u.metpy.convert_units(unit)
        v_units = v.metpy.convert_units(unit)

        wind_speed_units = wind_speed.metpy.convert_units(unit)

    if "wind_vectors" in products:
        if slicer_idx > 0:
            if model == 'era5':
                obj = ax.quiver(lon[::slicer_idx], lat[::slicer_idx], 
                            u_units[::slicer_idx, ::slicer_idx], 
                            v_units[::slicer_idx, ::slicer_idx], 
                            color=kwargs.pop('wind_display_color', 'greenyellow'), 
                            transform=crs.PlateCarree(),
                            zorder=kwargs.pop('wind_display_zorder', SHAPE_ZORDER))
            else:
                obj = ax.quiver(lon[::slicer_idx, ::slicer_idx], lat[::slicer_idx, ::slicer_idx], 
                            u_units[::slicer_idx, ::slicer_idx], 
                            v_units[::slicer_idx, ::slicer_idx], 
                            color=kwargs.pop('wind_display_color', 'greenyellow'), 
                            transform=crs.PlateCarree(),
                            scale=0.5,
                            scale_units='dots',
                            zorder=kwargs.pop('wind_display_zorder', SHAPE_ZORDER))
        else:
            obj = ax.quiver(lon, lat, 
                            u_units, 
                            v_units, 
                            color=kwargs.pop('wind_display_color', 'greenyellow'),
                            transform=crs.PlateCarree(),
                            scale=0.5,
                            scale_units='dots',
                            zorder=kwargs.pop('wind_display_zorder', SHAPE_ZORDER))
        prodabbr += "_Wv"

    if "wind_barbs" in products:
        if slicer_idx > 0:
            if model == 'era5':
                obj = ax.barbs(lon[::slicer_idx], lat[::slicer_idx], 
                               u_units[::slicer_idx, ::slicer_idx], 
                               v_units[::slicer_idx, ::slicer_idx], 
                               color=kwargs.pop('wind_display_color', 'greenyellow'), 
                               length=kwargs.pop('wind_vector_length', 6), 
                               transform=crs.PlateCarree(),
                               zorder=kwargs.pop('wind_display_zorder', SHAPE_ZORDER))
            else:
                obj = ax.barbs(lon[::slicer_idx, ::slicer_idx], lat[::slicer_idx, ::slicer_idx], 
                               u_units[::slicer_idx, ::slicer_idx], 
                               v_units[::slicer_idx, ::slicer_idx], 
                               color=kwargs.pop('wind_display_color', 'greenyellow'), 
                               length=kwargs.pop('wind_vector_length', 6), 
                               transform=crs.PlateCarree(),
                               zorder=kwargs.pop('wind_display_zorder', SHAPE_ZORDER))
        else:
            obj = ax.barbs(lon, lat, 
                           u_units, 
                           v_units, 
                           color=kwargs.pop('wind_display_color', 'greenyellow'), 
                           length=kwargs.pop('wind_vector_length', 6), 
                           transform=crs.PlateCarree(),
                           zorder=kwargs.pop('wind_display_zorder', SHAPE_ZORDER))
        prodabbr += "_Wb"

    if "wind_values" in products:
        lat_val, lon_val = np.meshgrid(lat.values, lon.values)
        if slicer_idx > 0:
            for ws_row, lat_row, lon_row in zip(wind_speed_units, lat_val, lon_val):
                for this_ws, this_lat, this_lon in zip(ws_row[::slicer_idx], lat_row[::slicer_idx], lon_row[::slicer_idx]):
                    obj = ax.annotate(str(round(this_ws.magnitude, 1)), 
                                      (this_lon, this_lat), 
                                      horizontalalignment='right', 
                                      verticalalignment='top', 
                                      color=kwargs.pop('wind_display_color', 'white'), 
                                      clip_box=ax.bbox, 
                                      fontsize=kwargs.pop('wind_values_fontsize', 8), 
                                      transform=crs.PlateCarree(), 
                                      annotation_clip=False, 
                                      zorder=kwargs.pop('wind_display_zorder', SHAPE_ZORDER))
        else:
            for ws_row, lat_row, lon_row in zip(wind_speed_units, lat_val, lon_val):
                for this_ws, this_lat, this_lon in zip(ws_row, lat_row, lon_row,):
                    obj = ax.annotate(str(round(this_ws.magnitude, 1)), 
                                      (this_lon, this_lat), 
                                      horizontalalignment='right', 
                                      verticalalignment='top', 
                                      color=kwargs.pop('wind_display_color', 'white'), 
                                      clip_box=ax.bbox, 
                                      fontsize=kwargs.pop('wind_values_fontsize', 8), 
                                      transform=crs.PlateCarree(), 
                                      annotation_clip=False, 
                                      zorder=kwargs.pop('wind_display_zorder', SHAPE_ZORDER))
        prodabbr += "_Wd"

    # obj.set_path_effects([PathEffects.withStroke(linewidth=4, 
    #                                              foreground=kwargs.pop('wind_display_outline_coor', 'black'))])

    return prodstr, prodabbr

def draw_surface_products(fig: mpl.figure.Figure,
                          ax: mpl.axes.Axes,
                          products: list[str, ...],
                          model: str,
                          prodstr: str,
                          prodabbr: str,
                          unit_reg: pint.UnitRegistry,
                          fig_for_cb: mpl.figure.Figure = None,
                          ax_for_cb: mpl.axes.Axes = None,
                          **kwargs) -> tuple[str, str]:

    if fig_for_cb and ax_for_cb:
        cb_fig = fig_for_cb
        cb_ax = ax_for_cb

        cbar_opts = {'rotation': -90, 'color': 'k', 'labelpad': 20}
    else:
        cb_fig = fig
        cb_ax = None

        cbar_opts = {}


    meta_df = kwargs.get('sfc_data')
    if model == 'hrrr':
        vtable = HRRR_VARIABLE_TABLE
        lat = meta_df.variables['latitude'][:]
        lon = meta_df.variables['longitude'][:]
    elif model == 'era5':
        vtable = ERA5_VARIABLE_TABLE
        lat = meta_df.variables['latitude'][:]
        lon = meta_df.variables['longitude'][:]

    if "temp_cf" in products:
        levels = np.arange(kwargs.pop('temp_level_min', -30),
                           kwargs.pop('temp_level_max', 45),
                           kwargs.pop('temp_level_step', 2))

        data_file = kwargs.get(vtable['SFC_t'][2])
        var = data_file.variables[vtable['SFC_t'][0]][:]
        var_baseunits = units.Quantity(var, vtable['SFC_t'][1])

        default_unit = 'degC'

        if kwargs.get('temp_units'):
            unit = unit_reg.Unit(kwargs.pop('temp_units', default_unit))
            var_units = var_baseunits.to(unit_reg(unit))
        else:
            unit = unit_reg.Unit(default_unit)
            var_units = var_baseunits


        if VARIABLE_SMOOTHING == True:

            var_contourf = ax.contourf(lon, lat, 
                                       smooth_n_point(var_units.magnitude,
                                                      SMOOTHING_POINTS,
                                                      SMOOTHING_PASSES), 
                                       levels=levels, 
                                       transform=crs.PlateCarree(), 
                                       cmap=kwargs.pop('temp_pallete', cm.Spectral_r), 
                                       vmin=levels.min(),
                                       vmax=levels.max(),
                                       extend='both',
                                       zorder=kwargs.pop('cf_zorder', CONTOURF_ZORDER))
        else:
            var_contourf = ax.contourf(lon, lat, 
                                       var_units.magnitude, 
                                       levels=levels, 
                                       transform=crs.PlateCarree(), 
                                       cmap=kwargs.pop('temp_pallete', cm.Spectral_r), 
                                       vmin=levels.min(),
                                       vmax=levels.max(),
                                       extend='both',
                                       zorder=kwargs.pop('cf_zorder', CONTOURF_ZORDER))
        cb = cb_fig.colorbar(var_contourf, 
                          ax=ax, 
                          cax=cb_ax,
                          ticks=levels[::2], 
                          orientation=kwargs.pop('colorbar_orientation', 'vertical'), 
                          shrink=0.77)
        cb.set_label(f'2m Temperature ({unit:~C})', size='x-large', **cbar_opts)

        prodstr += f"#2m Temperature ({unit:~C}, shaded)"
        prodabbr += "_T"

    if "rh_cf" in products:
        levels = np.arange(kwargs.pop('rh_level_min', 0),
                           kwargs.pop('rh_level_max', 1.1),
                           kwargs.pop('rh_level_step', 0.1))

        data_file = kwargs.get(vtable['SFC_rh'][2])
        var = data_file.variables[vtable['SFC_rh'][0]][:]

        if VARIABLE_SMOOTHING == True:
            var_contourf = ax.contourf(lon, lat, 
                                       smooth_n_point(var, 
                                                      SMOOTHING_POINTS,
                                                      SMOOTHING_PASSES), 
                                       levels=levels, 
                                       transform=crs.PlateCarree(), 
                                       cmap=kwargs.pop('rh_pallete', cm.terrain_r), 
                                       extend='max',
                                       zorder=kwargs.pop('cf_zorder', CONTOURF_ZORDER))
        else:
            var_contourf = ax.contourf(lon, lat, 
                                       var, 
                                       levels=levels, 
                                       transform=crs.PlateCarree(), 
                                       cmap=kwargs.pop('rh_pallete', cm.terrain_r), 
                                       vmin=levels.min(),
                                       vmax=levels.max(),
                                       extend='max',
                                       zorder=kwargs.pop('cf_zorder', CONTOURF_ZORDER))
        cb = cb_fig.colorbar(var_contourf, 
                          ax=ax, 
                          cax=cb_ax,
                          ticks=levels[::2], 
                          orientation=kwargs.pop('colorbar_orientation', 'vertical'), 
                          shrink=0.77,
                          extendrect=True)
        cb.set_label('2m Relative Humidty (%)', size='x-large', **cbar_opts)

        prodstr += "#2m Relative Humidity (%, shaded)"
        prodabbr += "_RH"

    if "wind_cf" in products:
        levels = np.arange(kwargs.pop('wind_level_min', 0),
                           kwargs.pop('wind_level_max', 40),
                           kwargs.pop('wind_level_step', 2))

        data_file = kwargs.get(vtable['SFC_u'][2])

        base_unit = vtable['SFC_u'][1]
        default_unit = 'knots'

        u_baseunits = data_file[vtable['SFC_u'][0]].values
        v_baseunits = data_file[vtable['SFC_v'][0]].values

        var_baseunits = mpcalc.wind_speed(units.Quantity(u_baseunits, base_unit),
                                          units.Quantity(v_baseunits, base_unit))

        if kwargs.get('wind_units'):
            unit = unit_reg.Unit(kwargs.get('wind_units'))
            var_units = var_baseunits.to(unit)
        else:
            unit = unit_reg.Unit(default_unit)
            var_units = var_baseunits.to(unit)

        if VARIABLE_SMOOTHING == True:
            var_contourf = ax.contourf(lon, lat, 
                                       smooth_n_point(var_units, 
                                                      SMOOTHING_POINTS,
                                                      SMOOTHING_PASSES), 
                                       levels=levels, 
                                       transform=crs.PlateCarree(), 
                                       cmap=kwargs.pop('wind_pallete', cm.jet), 
                                       extend='max',
                                       zorder=kwargs.pop('cf_zorder', CONTOURF_ZORDER))
        else:
            var_contourf = ax.contourf(lon, lat, 
                                       var_units, 
                                       levels=levels, 
                                       transform=crs.PlateCarree(), 
                                       cmap=kwargs.pop('wind_pallete', cm.jet), 
                                       vmin=levels.min(),
                                       vmax=levels.max(),
                                       extend='max',
                                       zorder=kwargs.pop('cf_zorder', CONTOURF_ZORDER))
        cb = cb_fig.colorbar(var_contourf, 
                          ax=ax, 
                          cax=cb_ax,
                          ticks=levels[::2], 
                          orientation=kwargs.pop('colorbar_orientation', 'vertical'), 
                          shrink=0.77)
        cb.set_label(f'10m Wind Speed ({unit:~C})', size='x-large', **cbar_opts)

        prodstr += f"#10m Wind Speed ({unit:~C}, shaded)"
        prodabbr += "_WS"

    if "wd_cf" in products:
        levels = np.arange(kwargs.pop('wd_level_min', 0),
                           kwargs.pop('wd_level_max', 365),
                           kwargs.pop('wd_level_step', 15))

        data_file = kwargs.get(vtable['SFC_u'][2])
        wind_direction = mpcalc.wind_direction(data_file[vtable['SFC_u'][0]], data_file[vtable['SFC_v'][0]])

        default_unit = 'degrees'
        
        if kwargs.get('wd_units'):
            unit = unit_reg.Unit(kwargs.pop('wd_units', default_unit))
            var_units = wind_direction.metpy.convert_units(unit)
        else:
            unit = unit_reg.Unit(default_unit)
            var_units = wind_direction.metpy.convert_units(unit)

        if VARIABLE_SMOOTHING == True:
            var_contourf = ax.contourf(lon, lat, 
                                       smooth_n_point(var_units, 
                                                      SMOOTHING_POINTS,
                                                      SMOOTHING_PASSES), 
                                       levels=levels, 
                                       transform=crs.PlateCarree(), 
                                       cmap=kwargs.pop('wd_pallete', cm.twilight), 
                                       zorder=kwargs.pop('cf_zorder', CONTOURF_ZORDER))
        else:
            var_contourf = ax.contourf(lon, lat, 
                                       var_units, 
                                       levels=levels, 
                                       transform=crs.PlateCarree(), 
                                       cmap=kwargs.pop('wd_pallete', cm.twilight), 
                                       vmin=levels.min(),
                                       vmax=levels.max(),
                                       zorder=kwargs.pop('cf_zorder', CONTOURF_ZORDER))
        cb = cb_fig.colorbar(var_contourf, 
                          ax=ax, 
                          cax=cb_ax,
                          ticks=levels[::2], 
                          orientation=kwargs.pop('colorbar_orientation', 'vertical'), 
                          shrink=0.77)
        cb.set_label(f'10m Wind Direction ({unit:~C})', size='x-large', **cbar_opts)

        prodstr += f"#10m Wind Direction ({unit:~C}, shaded)"
        prodabbr += "_WD"

    return prodstr, prodabbr

def draw_surface_wind_display(fig: mpl.figure.Figure,
                      ax: mpl.axes.Axes,
                      m10_data: xr.Dataset | netCDF4.Dataset,
                      products: list[str, ...],
                      model: str,
                      prodstr: str,
                      prodabbr: str,
                      unit_reg: pint.UnitRegistry,
                      **kwargs) -> tuple[str, str]:

    slicer_idx = kwargs.pop('wind_display_indexing', 4)

    if model == 'hrrr':
        vtable = HRRR_VARIABLE_TABLE
        lat = m10_data.variables['latitude'][:]
        lon = m10_data.variables['longitude'][:]
    elif model == 'era5':
        vtable = ERA5_VARIABLE_TABLE
        lat = m10_data.variables['latitude'][:]
        lon = m10_data.variables['longitude'][:]

    u = m10_data[vtable['SFC_u'][0]]
    v = m10_data[vtable['SFC_v'][0]]

    wind_speed = mpcalc.wind_speed(m10_data[vtable['SFC_u'][0]], m10_data[vtable['SFC_v'][0]])

    default_unit = 'knots'

    if kwargs.get('wind_display_units'):
        unit = unit_reg.Unit(kwargs.pop('wind_display_units', default_unit))

        u_units = u.metpy.convert_units(unit)
        v_units = v.metpy.convert_units(unit)

        wind_speed_units = wind_speed.metpy.convert_units(unit)

    else:
        unit = unit_reg.Unit(default_unit)

        u_units = u.metpy.convert_units(default_unit)
        v_units = v.metpy.convert_units(default_unit)

        wind_speed_units = wind_speed.metpy.convert_units(default_unit)

    if "wind_vectors" in products:
        if slicer_idx > 0:
            if model == 'era5':
                obj = ax.quiver(lon[::slicer_idx], lat[::slicer_idx], 
                            u_units[::slicer_idx, ::slicer_idx], 
                            v_units[::slicer_idx, ::slicer_idx], 
                            color=kwargs.pop('wind_display_color', 'greenyellow'), 
                            transform=crs.PlateCarree(),
                            zorder=kwargs.pop('wind_display_zorder', SHAPE_ZORDER))
            else:
                obj = ax.quiver(lon[::slicer_idx, ::slicer_idx], lat[::slicer_idx, ::slicer_idx], 
                            u_units[::slicer_idx, ::slicer_idx], 
                            v_units[::slicer_idx, ::slicer_idx], 
                            color=kwargs.pop('wind_display_color', 'greenyellow'), 
                            transform=crs.PlateCarree(),
                            scale=0.5,
                            scale_units='dots',
                            zorder=kwargs.pop('wind_display_zorder', SHAPE_ZORDER))
        else:
            obj = ax.quiver(lon, lat, 
                            u_units, 
                            v_units, 
                            color=kwargs.pop('wind_display_color', 'greenyellow'),
                            transform=crs.PlateCarree(),
                            scale=0.5,
                            scale_units='dots',
                            zorder=kwargs.pop('wind_display_zorder', SHAPE_ZORDER))
        prodabbr += "_Wv"

    if "wind_barbs" in products:
        if slicer_idx > 0:
            if model == 'era5':
                obj = ax.barbs(lon[::slicer_idx], lat[::slicer_idx], 
                               u_units[::slicer_idx, ::slicer_idx], 
                               v_units[::slicer_idx, ::slicer_idx], 
                               color=kwargs.pop('wind_display_color', 'greenyellow'), 
                               length=kwargs.pop('wind_vector_length', 6), 
                               transform=crs.PlateCarree(),
                               zorder=kwargs.pop('wind_display_zorder', SHAPE_ZORDER))
            else:
                obj = ax.barbs(lon[::slicer_idx, ::slicer_idx], lat[::slicer_idx, ::slicer_idx], 
                               u_units[::slicer_idx, ::slicer_idx], 
                               v_units[::slicer_idx, ::slicer_idx], 
                               color=kwargs.pop('wind_display_color', 'greenyellow'), 
                               length=kwargs.pop('wind_vector_length', 6), 
                               transform=crs.PlateCarree(),
                               zorder=kwargs.pop('wind_display_zorder', SHAPE_ZORDER))
        else:
            obj = ax.barbs(lon, lat, 
                           u_units, 
                           v_units, 
                           color=kwargs.pop('wind_display_color', 'greenyellow'), 
                           length=kwargs.pop('wind_vector_length', 6), 
                           transform=crs.PlateCarree(),
                           zorder=kwargs.pop('wind_display_zorder', SHAPE_ZORDER))
        prodabbr += "_Wb"

    if "wind_values" in products:
        lat_val, lon_val = np.meshgrid(lat.values, lon.values)
        if slicer_idx > 0:
            for ws_row, lat_row, lon_row in zip(wind_speed_units, lat_val, lon_val):
                for this_ws, this_lat, this_lon in zip(ws_row[::slicer_idx], lat_row[::slicer_idx], lon_row[::slicer_idx]):
                    obj = ax.annotate(str(round(this_ws.magnitude, 1)), 
                                      (this_lon, this_lat), 
                                      horizontalalignment='right', 
                                      verticalalignment='top', 
                                      color=kwargs.pop('wind_display_color', 'white'), 
                                      clip_box=ax.bbox, 
                                      fontsize=kwargs.pop('wind_values_fontsize', 8), 
                                      transform=crs.PlateCarree(), 
                                      annotation_clip=False, 
                                      zorder=kwargs.pop('wind_display_zorder', SHAPE_ZORDER))
        else:
            for ws_row, lat_row, lon_row in zip(wind_speed_units, lat_val, lon_val):
                for this_ws, this_lat, this_lon in zip(ws_row, lat_row, lon_row,):
                    obj = ax.annotate(str(round(this_ws.magnitude, 1)), 
                                      (this_lon, this_lat), 
                                      horizontalalignment='right', 
                                      verticalalignment='top', 
                                      color=kwargs.pop('wind_display_color', 'white'), 
                                      clip_box=ax.bbox, 
                                      fontsize=kwargs.pop('wind_values_fontsize', 8), 
                                      transform=crs.PlateCarree(), 
                                      annotation_clip=False, 
                                      zorder=kwargs.pop('wind_display_zorder', SHAPE_ZORDER))
        prodabbr += "_Wd"

    # obj.set_path_effects([PathEffects.withStroke(linewidth=4, 
    #                                              foreground=kwargs.pop('wind_display_outline_coor', 'black'))])

    return prodstr, prodabbr

if __name__ == "__main__":

    parser = argparse.ArgumentParser(description=f'Plot GOES 16/17/18 GLM Data, either for single time periods or cumulative over all input data.')
    parser.add_argument('--level-range', 
                        help=f'specify range of levels to plot. default: {DEFAULT_LEVELS_TO_PLOT}', 
                        nargs=3, 
                        type=int,
                        metavar=('max', 'min', 'step'),
                        dest='levels_range',
                        default=None)
    parser.add_argument('--bbox', 
                        help='specify bounding box of the image (N-S-E-W, decimal lat/lon)', 
                        nargs=4, 
                        type=float,
                        metavar=('N', 'S', 'E', 'W'),
                        default=None)
    parser.add_argument('--center', 
                        help='specifiy center of the image and a radius to plot', 
                        nargs=3, 
                        type=float,
                        metavar=('lat', 'lon', 'rad'),
                        default=None)
    parser.add_argument('--points-from-file',
                        help='specify a CSV file to read in points from. CSV file must have no header and be in format: lat,lon,label',
                        type=str,
                        metavar='file_path',
                        dest='points_file',
                        default=None)
    parser.add_argument('--settings-from-file',
                        help='specify a JSON file to read in plot settings from.',
                        type=str,
                        metavar='file_path',
                        dest='settings_file',
                        default=None)
    parser.add_argument('--save-as-kmz', 
                        help='save outputs as georeferenced KMZ files instead of images', 
                        action='store_true',
                        dest='save_to_kmz',
                        default=False)
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
    parser.add_argument('variables',
                        help='a list of variables to plot. variable names are specified in the script',
                        type=str,
                        nargs='+',
                        metavar='var_name')

    args = parser.parse_args()
    print(args)

    ## Quality Checks/Sanitizing

    if args.bbox:
        bbox = args.bbox
    else:
        if args.center:
            center_lat = args.center[0]
            center_lon = args.center[1]
            rad = args.center[2]

            bbox = [center_lat+rad, center_lat-rad, center_lon+rad, center_lon-rad] #NSEW
        else: #no center or bbox params defined
            bbox = [51, 23.5, -64.4, -129.1] #default, CONUS NSEW

    if args.input_file_directory[-1:] != "/":
        input_dir = args.input_file_directory + "/"
    else:
        input_dir = args.input_file_directory

    print(input_dir)

    if args.save_directory[-1:] != "/":
        save_dir = args.save_directory + "/"
    else:
        save_dir = args.save_directory

    print(save_dir)

    if args.levels_range:
        plot_levels = range(args.levels_range[1], args.levels_range[0]+args.levels_range[2], args.levels_range[2])
    else:
        plot_levels = DEFAULT_LEVELS_TO_PLOT
    num_levels = len(plot_levels)

    if args.settings_file:
        try:
            with open(args.settings_file, newline='') as settingsjson:
                user_settings = json.load(settingsjson)
        except Exception as err:
            raise err
    else:
        user_settings = dict()

    if args.save_to_kmz:
        user_settings.update({'save_to_kmz': True})

    points = []
    if args.points_file:
        try:
            with open(args.points_file, newline='') as pointcsv:
                reader = csv.reader(pointcsv, delimiter=',')
                for row in reader:

                    point = (float(row[0]),
                             float(row[1]),
                             row[2] if len(row) > 2 else None,
                             row[3] if len(row) > 3 else None,
                             row[4] if len(row) > 4 else None,
                             row[5] if len(row) > 5 else None)
                    points += [point]
        except Exception as err:
            raise err

    if args.model == 'era5':
        input_files = sorted(glob.glob(f'{input_dir}*UA_ERA5.grib'))
    elif args.model == 'era5_sfc':
        input_files = sorted(glob.glob(f'{input_dir}*SFC_ERA5.grib'))
    elif args.model == 'era5_land':
        input_files = sorted(glob.glob(f'{input_dir}*ERA5Land.grib'))
    elif args.model == 'hrrr':
        input_files = sorted(glob.glob(f'{input_dir}*.grib2'))
    elif args.model == 'wrf':
        input_files = sorted(glob.glob(f'{input_dir}wrfout*'))
    elif args.model == 'gfs':
        input_files = sorted(glob.glob(f'{input_dir}*pgrb2*'))

    unit_registry = pint.UnitRegistry()

    num_input_files = len(input_files)
    num_plot_steps = num_input_files * num_levels

    plot_times = []
    start_of_period = True

    with tqdm(miniters=0, total=num_plot_steps, desc=f'Plotting {len(args.variables)} variables on {num_levels} levels...', ascii=" ░▒▓█") as progress:
        for input_file_path in input_files:
            for level in plot_levels:
                if args.model == "hrrr":
                    plot_plan_view_hrrr(input_file_path, 
                                        save_dir, 
                                        level, 
                                        args.variables, 
                                        points, 
                                        bbox, 
                                        unit_registry,
                                        **user_settings)
                elif args.model == "era5":
                    plot_plan_view_era5(input_file_path, 
                                        save_dir, 
                                        level, 
                                        args.variables, 
                                        points, 
                                        bbox, 
                                        **user_settings)
                elif args.model == "wrf":
                    pass
                progress.update()

    print("Done!")






