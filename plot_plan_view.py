'''
IMPLEMENTED VARIABLES AND DISPLAY TYPES:

    Shaded/Filled contours: 
    Cloud Cover ['cldcvr_cf'], Relative Humidity (%) ['rh_cf'], Wind Speed (kts) ['wind_cf'],
    Vertical Velocity (ft/min) ['vv_cf']

    Contours:
    Temperature (°C) ['temp_c'], Geopotential Height (m) ['gpm_c']

    Vectors:
    Wind (kts) ['wind_vectors']

    Barbs:
    Wind (kts) ['wind_barbs']

    Grid Point Values:
    Wind (kts) ['wind_values']
    
'''

## Turn on N-point variable VARIABLE_SMOOTHING. Only applies to shaded/filled contour variables.
## (Default: False)
VARIABLE_SMOOTHING = False

## Amount of points to use for N-point VARIABLE_SMOOTHING. Must be 5 or 9.
## (Default: 9)
SMOOTHING_POINTS = 9

## Amount of passes to use for N-point VARIABLE_SMOOTHING.
## (Default: 1)
SMOOTHING_PASSES = 1

DEFAULT_LEVELS_TO_PLOT = [1000, 950, 900, 850, 700, 500, 300, 250]

## Set the DPI of the saved output file
## (Default: 300)
FILE_DPI = 300

TOWN_SCALE_RANK = 14

'''
Used to set the style of any user input points shown on the map.
POINT_STYLE is used for the point(s), and POINT_LABEL_STYLE is used for any
point label(s) given. No labels are shown if POINT_LABEL_VISIBLE is False.
'''
POINT_STYLE = {'color': 'black', 'markersize': 8,}
POINT_LABEL_VISIBLE = True
DRAW_LABEL_ARROWS = True
POINT_LABEL_STYLE = {'color': 'red', 'fontsize': 14, 'fontweight': 'bold'}

## Set the positioning of the labels relative to the points being plotted
## (Default: X=0.2, Y=0.1)
Y_LABEL_OFFSET = -0.2
X_LABEL_OFFSET = 0.2














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
                        geopotential_to_height, vertical_velocity_pressure)
import metpy.calc as mpcalc
import netCDF4
from tqdm.auto import tqdm
from adjustText import adjust_text 

from __internal_funcs import (plot_towns, draw_logo, 
                              plot_points, define_hi_res_fig,
                              define_gearth_compat_fig, save_figs_to_kml)

HRRR_VARIABLE_TABLE = {
    'temp': 't',
    'dpt': 'dpt',
    'rh': 'r',
    'gpm': 'gh',
    'vv': 'w',
    'cldcvr': 'cc',
    'u': 'u',
    'v': 'v',
    'mxgrat' : 'q',
    'pressure': 'isobaricInhPa',
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
                    x_label_offset=X_LABEL_OFFSET,
                    y_label_offset=Y_LABEL_OFFSET,
                    draw_labels=POINT_LABEL_VISIBLE,
                    draw_arrows=DRAW_LABEL_ARROWS,
                    point_style=POINT_STYLE,
                    point_label_style=POINT_LABEL_STYLE)

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
                   scale_rank=TOWN_SCALE_RANK)

        draw_logo(ax)
        
        descstr = f"{str(data['isobaricInhPa'].values)} hPa Level\nValid at {np.datetime_as_string(data['valid_time'].values, timezone='UTC')[:-11].replace('T', ' ')} UTC"
        ax.set_title(f'HRRR Reanalysis {titlestr}', loc='left', fontweight='bold', fontsize=15)
        plt.title(descstr, loc='right')

        file_name = "PlanView_" + str(sel_level) + filestr + "_" + file_day + "_" + file_time +"_ERA5"
        
        if not os.path.exists(save_dir):
            os.makedirs(save_dir)

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
                    x_label_offset=X_LABEL_OFFSET,
                    y_label_offset=Y_LABEL_OFFSET,
                    draw_labels=POINT_LABEL_VISIBLE,
                    draw_arrows=DRAW_LABEL_ARROWS,
                    point_style=POINT_STYLE,
                    point_label_style=POINT_LABEL_STYLE)
    
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
               scale_rank=TOWN_SCALE_RANK)

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

def draw_contourf_lines(fig: mpl.figure.Figure,
                        ax: mpl.axes.Axes,
                        data: xr.Dataset | netCDF4.Dataset,
                        sel_level: int,
                        products: list[str, ...],
                        model: str,
                        prodstr: str,
                        prodabbr: str,
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


    if "cldcvr_cf" in products:
        levels = np.arange(kwargs.pop('cldcvr_level_min', 0),
                           kwargs.pop('cldcvr_level_max', 1.1),
                           kwargs.pop('cldcvr_level_step', 0.1))

        if VARIABLE_SMOOTHING == True:
            var_contourf = ax.contourf(lon, lat, 
                                       smooth_n_point(data.variables[vtable['cldcvr']][:],
                                                      SMOOTHING_POINTS,
                                                      SMOOTHING_PASSES), 
                                       levels=levels, 
                                       transform=crs.PlateCarree(), 
                                       cmap=kwargs.pop('cldcvr_pallete', cm.Blues), 
                                       zorder=kwargs.pop('cf_zorder', CONTOURF_ZORDER))
        else:
            var_contourf = ax.contourf(lon, lat, 
                                       data.variables[vtable['cldcvr']][:], 
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
        cb.set_label('Fraction of Cloud Cover', size='x-large', **cbar_opts)

        prodstr += "#Fraction of Cloud Cover (shaded)"
        prodabbr += "_CC"

    if "rh_cf" in products:
        levels = np.arange(kwargs.pop('rh_level_min', 0),
                           kwargs.pop('rh_level_max', 1.1),
                           kwargs.pop('rh_level_step', 0.1))

        if VARIABLE_SMOOTHING == True:
            var_contourf = ax.contourf(lon, lat, 
                                       smooth_n_point(data.variables[vtable['rh']][:], 
                                                      SMOOTHING_POINTS,
                                                      SMOOTHING_PASSES), 
                                       levels=levels, 
                                       transform=crs.PlateCarree(), 
                                       cmap=kwargs.pop('rh_pallete', cm.terrain_r), 
                                       extend='max',
                                       zorder=kwargs.pop('cf_zorder', CONTOURF_ZORDER))
        else:
            var_contourf = ax.contourf(lon, lat, 
                                       data.variables[vtable['rh']][:], 
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
        cb.set_label('Relative Humidty (%)', size='x-large', **cbar_opts)

        prodstr += "#Relative Humidity (%, shaded)"
        prodabbr += "_RH"

    if "wind_cf" in products:
        levels = np.arange(kwargs.pop('wind_level_min', 0),
                           kwargs.pop('wind_level_max', 40 if sel_level > 500 else 100),
                           kwargs.pop('wind_level_step', 2 if sel_level > 500 else 5))

        wind_speed = mpcalc.wind_speed(data[vtable['u']], data[vtable['v']]).metpy.convert_units('knots')

        if VARIABLE_SMOOTHING == True:
            var_contourf = ax.contourf(lon, lat, 
                                       smooth_n_point(wind_speed, 
                                                      SMOOTHING_POINTS,
                                                      SMOOTHING_PASSES), 
                                       levels=levels, 
                                       transform=crs.PlateCarree(), 
                                       cmap=kwargs.pop('wind_pallete', cm.jet), 
                                       extend='max',
                                       zorder=kwargs.pop('cf_zorder', CONTOURF_ZORDER))
        else:
            var_contourf = ax.contourf(lon, lat, 
                                       wind_speed, 
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
        cb.set_label('Wind Speed (kts)', size='x-large', **cbar_opts)

        prodstr += "#Wind Speed (kts, shaded)"
        prodabbr += "_WS"

    if "vv_cf" in products:
        levels = np.arange(kwargs.pop('vv_level_min', -450),
                           kwargs.pop('vv_level_max', 450),
                           kwargs.pop('vv_level_step', 10))

        if model == "hrrr":
            vv = vertical_velocity(data[vtable['vv']],
                                   data[vtable['pres']],
                                   data[vtable['temp']],
                                   data[vtable['mxgrat']])

            vv = vv.magnitude * 196.9 #m/s to ft/min

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

        gpm = geopotential_to_height(data[vtable['gpm']])

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

    ukt = data[vtable['u']]
    vkt = data[vtable['v']]

    wind_speed = mpcalc.wind_speed(data[vtable['u']], data[vtable['v']]).metpy.convert_units('knots')

    if "wind_vectors" in products:
        if slicer_idx > 0:
            if model == 'era5':
                obj = ax.quiver(lon[::slicer_idx], lat[::slicer_idx], 
                            ukt[::slicer_idx, ::slicer_idx], 
                            vkt[::slicer_idx, ::slicer_idx], 
                            color=kwargs.pop('wind_display_color', 'greenyellow'), 
                            transform=crs.PlateCarree(),
                            zorder=kwargs.pop('wind_display_zorder', SHAPE_ZORDER))
            else:
                obj = ax.quiver(lon[::slicer_idx, ::slicer_idx], lat[::slicer_idx, ::slicer_idx], 
                            ukt[::slicer_idx, ::slicer_idx], 
                            vkt[::slicer_idx, ::slicer_idx], 
                            color=kwargs.pop('wind_display_color', 'greenyellow'), 
                            transform=crs.PlateCarree(),
                            zorder=kwargs.pop('wind_display_zorder', SHAPE_ZORDER))
        else:
            obj = ax.quiver(lon, lat, 
                            ukt, 
                            vkt, 
                            color=kwargs.pop('wind_display_color', 'greenyellow'),
                            transform=crs.PlateCarree(),
                            zorder=kwargs.pop('wind_display_zorder', SHAPE_ZORDER))
        prodabbr += "_Wv"

    if "wind_barbs" in products:
        if slicer_idx > 0:
            if model == 'era5':
                obj = ax.barbs(lon[::slicer_idx], lat[::slicer_idx], 
                               ukt[::slicer_idx, ::slicer_idx], 
                               vkt[::slicer_idx, ::slicer_idx], 
                               color=kwargs.pop('wind_display_color', 'greenyellow'), 
                               length=kwargs.pop('wind_vector_length', 6), 
                               transform=crs.PlateCarree(),
                               zorder=kwargs.pop('wind_display_zorder', SHAPE_ZORDER))
            else:
                obj = ax.barbs(lon[::slicer_idx, ::slicer_idx], lat[::slicer_idx, ::slicer_idx], 
                               ukt[::slicer_idx, ::slicer_idx], 
                               vkt[::slicer_idx, ::slicer_idx], 
                               color=kwargs.pop('wind_display_color', 'greenyellow'), 
                               length=kwargs.pop('wind_vector_length', 6), 
                               transform=crs.PlateCarree(),
                               zorder=kwargs.pop('wind_display_zorder', SHAPE_ZORDER))
        else:
            obj = ax.barbs(lon, lat, 
                           ukt, 
                           vkt, 
                           color=kwargs.pop('wind_display_color', 'greenyellow'), 
                           length=kwargs.pop('wind_vector_length', 6), 
                           transform=crs.PlateCarree(),
                           zorder=kwargs.pop('wind_display_zorder', SHAPE_ZORDER))
        prodabbr += "_Wb"

    if "wind_values" in products:
        lat_val, lon_val = np.meshgrid(lat.values, lon.values)
        if slicer_idx > 0:
            for ws_row, lat_row, lon_row in zip(ws, lat_val, lon_val):
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
            for ws_row, lat_row, lon_row in zip(ws, lat_val, lon_val):
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

    obj.set_path_effects([PathEffects.withStroke(linewidth=4, 
                                                 foreground=kwargs.pop('wind_display_outline_coor', 'black'))])

    return prodstr, prodabbr

if __name__ == "__main__":

    parser = argparse.ArgumentParser(description=f'Plot GOES 16/17/18 GLM Data, either for single time periods or cumulative over all input data.')
    parser.add_argument('--level-range', 
                        help=f'specify range of levels to plot. default: {DEFAULT_LEVELS_TO_PLOT}', 
                        nargs=2, 
                        type=int,
                        metavar=('max', 'min'),
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
        plot_levels = range(args.levels_range[1], args.levels_range[0]+25, 25)
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
                    if row[2] == "" or row[2].isspace():
                        marker = 'x'
                    else:
                        marker = row[2]

                    if row[3] == "" or row[3].isspace():
                        point = (float(row[0]), float(row[1]), marker, None)
                    else:
                        point = (float(row[0]), float(row[1]), marker, row[3])
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






