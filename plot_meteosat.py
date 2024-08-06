## Set to True if you want more towns plotted on the map
## (Default: False)
TOWN_SCALE_RANK = 14

## Set the DPI of the saved output file
FILE_DPI = 300

## Set the default color pallete to be used when plotting single-band imagery.
## (Default: "Greys_r")
DEFAULT_SB_PALLETE = "Greys_r"

DEFAULT_COMPOSITE_PRODUCT = "true_color"

'''
Used to set the style of any user input points shown on the map.
POINT_STYLE is used for the point(s), and POINT_LABEL_STYLE is used for any
point label(s) given. No labels are shown if POINT_LABEL_VISIBLE is False.
'''
POINT_STYLE = {'color': 'black', 'markersize': 8}
POINT_LABEL_VISIBLE = True
POINT_LABEL_STYLE = {'color': 'black', 'fontsize': 12}
DRAW_LABEL_ARROWS = True


## Set the positioning of the labels relative to the points being plotted
## (Default: X=0.2, Y=0.1)
Y_LABEL_OFFSET = 0.0
X_LABEL_OFFSET = 0.35

## Set the positioning of the pixel values relative to the upper-left corner
## of the pixel.
## (Default: X=5, Y=10)
X_PIX_VAL_OFFSET = 0
Y_PIX_VAL_OFFSET = 0











# 6/13/24 - Added cmasher (cmr) import for cmasher colormap support




import sys
import os
import glob
from datetime import datetime
import re
import glob
import argparse
import csv
import warnings
warnings.simplefilter("ignore")

import xarray as xr
import numpy as np
import numpy.ma as ma
import metpy
import matplotlib
from matplotlib import cm
import matplotlib.pyplot as plt
import matplotlib.patheffects as PathEffects
from matplotlib.offsetbox import OffsetImage, AnnotationBbox
import cartopy 
import cartopy.crs as crs
import cartopy.crs as ccrs
import cartopy.feature as cfeat
import cartopy.io.shapereader as shpreader
from cartopy.feature import NaturalEarthFeature
from tqdm.auto import tqdm
from adjustText import adjust_text
import cmasher as cmr
import satpy

from __internal_funcs import (plot_towns, draw_logo, plot_points,
                              define_gearth_compat_fig, save_figs_to_kml,
                              define_hi_res_fig)

matplotlib.use('agg')
np.seterr(divide = 'ignore', invalid='ignore')

#global params
DEFAULT_BBOX = [-129.1, -64.4, 23.5, 51] #CONUS - WESN

DEF_GAMMA = 2.2

GEOG_VISIBLE = True
GEOG_DRAW_STATES = True
GEOG_DRAW_COASTLINES = True

COLORBAR_VISIBLE = False
COLORBAR_LABEL = 'Radiance Temperature (K)'

IMPLEMENTED_COMPOSITES = {
    'dust_rgb' : 'Dust RGB',
    'nt_microphysics' : 'Nighttime Microphysics',
    'day_land_cloud_fire' : 'Day Land Cloud Fire',
    'day_cloud_phase' : 'Day Cloud Phase',
    'true_color' : 'True Color'
}

IMPLEMENTED_BAND_NAMES = {
    '1': 'Visible (Red) (0.6 µm, 3 km)',
    '2': 'Visible (Green) (0.8 µm, 3 km)',
    '3': 'Near-Infrared (Blue) (1.6 µm, 3 km)',
    '4': 'Infrared (3.9 µm, 3 km)',
    '5': 'High-Level Water Vapor (6.2 µm, 3 km)',
    '6': 'Mid-Level Water Vapor (7.3 µm, 3 km)',
    '7': 'Infrared (8.7 µm, 3 km)',
    '8': 'Infrared (9.7 µm, 3 km)',
    '9': 'Infrared (10.8 µm, 3 km)',
    '10': 'Infrared (12.0 µm, 3 km)',
    '11': 'Infrared (13.4 µm, 3 km)',
    '12': 'High-Resolution Visible (0.3 - 1 µm, 1 km)',
}

IMPLEMENTED_BAND_IDS = {
    '1': 'VIS006',
    '2': 'VIS008',
    '3': 'IR_016',
    '4': 'IR_039',
    '5': 'WV_062',
    '6': 'WV_073',
    '7': 'IR_087',
    '8': 'IR_097',
    '9': 'IR_108',
    '10': 'IR_120',
    '11': 'IR_134',
    '12': 'HRV',
}

def plot_single_band_meteosat(file_path: str, 
                              save_dir: str, 
                              band: int, 
                              points: list[tuple[float, float, str], ...],
                              bbox: list[float, float, float, float] = DEFAULT_BBOX, 
                              pallete: str = DEFAULT_SB_PALLETE, 
                              **kwargs):
    """
    Using a bounding box, plots a single satellite band and any points.

    Band is specified as a int between 1-16, corresponding to GOES
    band id.

    Bounding box is specified via a list of length 4:
    [ll corner lat, ll corner lon, ur corner lat, ur corner lon]

    Points are specified as a list of tuples, including label:
    [(lat1, lon1, label1), (lat2, lon2, label2), ...]

    Pallete is specified as a matplotlib compatable color pallete
    name string. Suffix '_r' to reverse the color pallete

    Returns: Success code, time (in seconds) for function to run,
             path to file
    """

    try:
        data = satpy.Scene(reader="seviri_l1b_native",
                           filenames=[file_path],
                           reader_kwargs={'fill_disk': True})
    except Exception as e:
        raise e

    sel_band_name = IMPLEMENTED_BAND_NAMES[str(band)]
    sel_band_id = IMPLEMENTED_BAND_IDS[str(band)]

    data.load([sel_band_id])
    data = data.crop(ll_bbox=(bbox[3], bbox[1], bbox[2], bbox[0]))
    data = data.resample(data[sel_band_id].attrs['area'], resampler='native')

    
    sel_band = data[sel_band_id]

    scan_start = sel_band.time_parameters['observation_start_time']
    scan_end = sel_band.time_parameters['observation_end_time']
    sat_sensor = sel_band.attrs['sensor']
    sat_name = sel_band.attrs['platform_name']
    sat_id = sat_name.replace("Meteosat-", "MSG")

    
    geog_data = sel_band.attrs["area"].to_cartopy_crs()

    x = sel_band.x
    y = sel_band.y

    img_data = sel_band.values

    img_extent = sel_band.attrs["area"].area_extent

    proj_info = {}
    for entry in sel_band.attrs['area'].proj4_string.split(" +"):
        items = entry.replace('+','').split('=')
        proj_info.update({items[0]: items[1] if len(items) > 1 else None})

    # pix_lats, pix_lons, pix_vals = _calculate_pixel_lat_lon(x.values, y.values, sel_band, sel_band.attrs["area"])

    pix_lats, pix_lons = sel_band.attrs["area"].get_lonlats()

    # x_range = range(0, len(proj_x))
    # y_range = range(0, len(proj_y))
    # x_idx, y_idx = np.meshgrid(x_range, y_range)

    # for x, y in zip(x_idx, y_idx):
    #     lon, lat = proj_area_def.get_lonlat_from_array_coordinates(x, y)
    #     lats += [lat]
    #     lons += [lon]

    #     x, y = proj_area_def.lonlat2colrow(lon, lat)
    #     val = proj_data.values[y, x]
    #     vals += [val]
    
    #x, y = np.meshgrid(x, y)

    #DATA CORRECTIONS

    #Bands 1-6 are output in "Reflectance Factor", from 0-1 and can be clip corrected
    #Bands 7-16 are output in "Brightness Temperature", and cannot be clip or gamma corrected
    if band >= 7: 
        correct_clip = False
        correct_gamma = False

        img_vmin = 210
        img_vmax = 230

        cbar_label = "Radiance Temperature (K)"

    else:
        correct_clip = True
        correct_gamma = True

        img_vmin = 0
        img_vmax = 1

        cbar_label = "Apparent Brightness"

    if kwargs.get('save_to_kmz'):  
        fig, ax, cbfig, cbax = define_gearth_compat_fig((bbox[2], bbox[3]),
                                                        (bbox[0], bbox[1]))

        cbar_opts = {'rotation': -90, 'color': 'k', 'labelpad': 20}

    else:
        fig, ax = define_hi_res_fig((bbox[2], bbox[3]),
                                    (bbox[0], bbox[1]))

        cbfig = fig
        cbax = None

        cbar_opts = {}


    
    img = ax.imshow(img_data, 
              origin='upper',
              extent=(img_extent[0], img_extent[2], img_extent[1], img_extent[3]),#(x.min().values, x.max().values, y.min().values, y.max().values),
              transform=geog_data,
              interpolation='none',
              # vmin = 180,
              # vmax = 330,
              cmap = pallete)

    
    if kwargs.get('pixel_value'):

        for lat, lon, val in zip(pix_lats.flatten(), pix_lons.flatten(), img_data.flatten()):
            ax.annotate(str(round(val, 1)),
                        (lat, lon),
                        xytext=(X_PIX_VAL_OFFSET,
                                Y_PIX_VAL_OFFSET),
                        textcoords='offset points',
                        horizontalalignment='center',
                        verticalalignment='center', 
                        color='black',
                        clip_box=ax.bbox,
                        fontsize=6,
                        transform=geog_data._as_mpl_transform(ax),
                        annotation_clip=False, 
                        zorder=30)

        # for lat_row, lon_row, val_row in zip(pix_lats, pix_lons, reversed(img_data)):
        #     for lat, lon, val in zip(lat_row, lon_row, reversed(val_row)):
        #         if not (np.isnan(lat) or np.isnan(lon)):

        #             #To compensate for projection, add an adjustment factor to the pixel values.
        #             #The adjustment factor starts large, in the upper left corner, and decreases
        #             #toward the right and the bottom.

        #             ax.annotate(str(round(val, 1)), 
        #                         (lon, lat),
        #                         xytext=(X_PIX_VAL_OFFSET,
        #                                 Y_PIX_VAL_OFFSET),
        #                         textcoords='offset points',
        #                         horizontalalignment='center',
        #                         verticalalignment='center', 
        #                         color='black',
        #                         clip_box=ax.bbox,
        #                         fontsize=6,
        #                         transform=crs.Geostationary()._as_mpl_transform(ax),
        #                         annotation_clip=False, 
        #                         zorder=30)

    if kwargs.get('colorbar_visible'):
        cb = cbfig.colorbar(img,
                        ax=ax,
                       cax=cbax,
                       orientation = "horizontal", 
                       pad=.05,
                       shrink=0.7,
                       use_gridspec=False)
        cb.set_label(cbar_label, size='x-large', **cbar_opts)

    #POINT DRAWING

    if points:
        plot_points(plt, ax,
                    points,
                    x_label_offset=X_LABEL_OFFSET,
                    y_label_offset=Y_LABEL_OFFSET,
                    draw_labels=POINT_LABEL_VISIBLE,
                    draw_arrows=DRAW_LABEL_ARROWS,
                    point_style=POINT_STYLE,
                    point_label_style=POINT_LABEL_STYLE)

    
    if kwargs.get('save_to_kmz'):
        
        file_name = sat_id + "_" + sel_band_id + "_" + scan_end.strftime('%Y%m%d_%H%M%S%Z')
        layer_name = f"{sat_name} {sat_sensor.upper()} image at " + scan_end.strftime('%d %B %Y %H:%M UTC')
        layer_desc = sel_band_name

        if not os.path.exists(save_dir):
            os.makedirs(save_dir)

        dest_path = os.path.join(save_dir, file_name + ".kmz")

        save_figs_to_kml(dest_path,
                         (bbox[2], bbox[3]),
                         (bbox[0], bbox[1]),
                         fig,
                         [layer_name],
                         [layer_desc],
                         colorbar_fig=cbfig if kwargs.get('colorbar_visible') else None)

        plt.close(fig)
        if cbfig:
            plt.close(cbfig)

    else:
        plot_towns(ax, 
                  (bbox[2], bbox[3]),
                  (bbox[0], bbox[1]), 
                  scale_rank=TOWN_SCALE_RANK)
        draw_logo(ax)

        title_info = f"{sat_name} {sat_sensor.upper()} imagery\n{sel_band_name}"
        
        ax.set_title(title_info, loc='left', fontweight='bold', fontsize=15)
        ax.set_title('{}'.format(scan_end.strftime('%d %B %Y %H:%M UTC ')), loc='right')
        
        file_name = sat_id + "_" + sel_band_id + "_" + scan_end.strftime('%Y%m%d_%H%M%S%Z')
                
        if not os.path.exists(save_dir):
            os.makedirs(save_dir)

        dest_path = os.path.join(save_dir, file_name + ".png")

        fig.savefig(dest_path, bbox_inches="tight", dpi=FILE_DPI)
        plt.close(fig)

def plot_composite_meteosat(file_path: str, 
                            save_dir: str,
                            points: list[tuple[float, float, str], ...], 
                            bbox: list[float, float, float, float] = DEFAULT_BBOX, 
                            product: str = DEFAULT_COMPOSITE_PRODUCT,  
                            **kwargs):
    """
    Using a bounding box, plots a single satellite band and any points.

    Band is specified as a int between 1-16, corresponding to GOES
    band id.

    Bounding box is specified via a list of length 4:
    [ll corner lat, ll corner lon, ur corner lat, ur corner lon]

    Points are specified as a list of tuples, including label:
    [(lat1, lon1, label1), (lat2, lon2, label2), ...]

    Pallete is specified as a matplotlib compatable color pallete
    name string. Suffix '_r' to reverse the color pallete

    Returns: Success code, time (in seconds) for function to run,
             path to file
    """

    try:
        data = xr.open_dataset(file_path, engine="netcdf4")
    except Exception as e:
        raise e

    scan_start = datetime.strptime(data.time_coverage_start, '%Y-%m-%dT%H:%M:%S.%fZ')
    scan_end = datetime.strptime(data.time_coverage_end, '%Y-%m-%dT%H:%M:%S.%fZ')
    file_created = datetime.strptime(data.date_created, '%Y-%m-%dT%H:%M:%S.%fZ')
    orbital_slot = data.orbital_slot #GOES-East, GOES-West, GOES-Test, etc.
    sat_id = data.platform_ID #G18, G17, G16, etc.

    red, green, blue, pallete, human_product_name = _calculate_composite_product_data(data, product)

    data = data.metpy.parse_cf('CMI_C02')
    geog_data = data.metpy.cartopy_crs
    x = data.x
    y = data.y

    correct_clip = True #I think this can be removed but needs more testing

    if correct_clip:
        red = np.clip(red, 0, 1)
        green = np.clip(green, 0, 1)
        blue = np.clip(blue, 0, 1)

    rgb_composite = np.stack([red, green, blue], axis=2) #Create the composite image data using calculated data above

    bbox_wesn = [bbox[3], bbox[2], bbox[1], bbox[0]]
    title_info = orbital_slot.replace("-Test", "") + " (" + sat_id.replace("G", "GOES-") + ")\n" + human_product_name + " Composite"

    if kwargs.get('save_to_kmz'):  
        fig, ax, cbfig, cbax = define_gearth_compat_fig((bbox[2], bbox[3]),
                                                        (bbox[0], bbox[1]))

        cbar_opts = {'rotation': -90, 'color': 'k', 'labelpad': 20}

    else:
        fig, ax = define_hi_res_fig((bbox[2], bbox[3]),
                                    (bbox[0], bbox[1]))

        cbfig = fig
        cbax = None

        cbar_opts = {}
        
    ax.imshow(rgb_composite, 
               origin='upper',
               extent=(x.min(), x.max(), y.min(), y.max()),
               transform=geog_data,
               interpolation='none')
    
    if kwargs.get('colorbar_visible'):
        cb = cbfig.colorbar(ax=ax,
                       cax=cbax,
                       orientation = "horizontal", 
                       pad=.05,
                       shrink=0.7,
                       use_gridspec=False)
        cb.set_label(cbar_label, size='x-large', **cbar_opts)

        #POINT DRAWING

    if points:
        plot_points(plt, ax,
            points,
            x_label_offset=X_LABEL_OFFSET,
            y_label_offset=Y_LABEL_OFFSET,
            draw_labels=POINT_LABEL_VISIBLE,
            draw_arrows=DRAW_LABEL_ARROWS,
            point_style=POINT_STYLE,
            point_label_style=POINT_LABEL_STYLE)

    if kwargs.get('save_to_kmz'):

        file_name = sat_id + "_" + sel_band_str + "_" + scan_end.strftime('%Y%m%d_%H%M%S%Z')
        layer_name = sat_id.replace("G", "GOES-") + " image at " + scan_end.strftime('%d %B %Y %H:%M UTC ')
        layer_desc = f'{human_product_name} Composite'

        if not os.path.exists(save_dir):
            os.makedirs(save_dir)

        dest_path = os.path.join(save_dir, file_name + ".kmz")

        save_figs_to_kml(dest_path,
                         (bbox[2], bbox[3]),
                         (bbox[0], bbox[1]),
                         fig,
                         [layer_name],
                         [layer_desc],
                         colorbar_fig=cbfig if kwargs.get('colorbar_visible') else None)

        plt.close(fig)
        if cbfig:
            plt.close(cbfig)

    else:
    
        plot_towns(ax, 
                  (bbox[2], bbox[3]),
                  (bbox[0], bbox[1]), 
                  scale_rank=TOWN_SCALE_RANK)
        draw_logo(ax)
    
        
        title_info = orbital_slot.replace("-Test", "") + " (" + sat_id.replace("G", "GOES-") + ")\n" + human_product_name + " Composite"
        ax.set_title(title_info, loc='left', fontweight='bold', fontsize=15)
        ax.set_title('{}'.format(scan_end.strftime('%d %B %Y %H:%M:%S UTC ')), loc='right')
        
        file_name = sat_id + "_" + product + "_" + scan_end.strftime('%Y%m%d_%H%M%S%Z')
        
        if not os.path.exists(save_dir):
            os.makedirs(save_dir)

        dest_path = os.path.join(save_dir, file_name + ".png")

        fig.savefig(dest_path, 
                    bbox_inches="tight", 
                    dpi=FILE_DPI)

        plt.close(fig)

# def _calculate_pixel_lat_lon(proj_x, proj_y, proj_data, proj_area_def):

#     lons, lats = proj_area_def.get_lonlat_from_projection_coordinates(proj_x, proj_y)


#     for lon, lat in zip(lons, lats):
#         x, y = proj_area_def.lonlat2colrow(lon, lat)
#         val = proj_data.values[y, x]
#         vals += [val]

#     return lats, lons, vals

def _calculate_pixel_lat_lon(proj_x, proj_y, proj_data, proj_area_def):

    lats = []
    lons = []
    vals = []

    x_range = range(0, len(proj_x))
    y_range = range(0, len(proj_y))
    x_idx, y_idx = np.meshgrid(x_range, y_range)

    for x, y in zip(x_idx, y_idx):
        lon, lat = proj_area_def.get_lonlat_from_array_coordinates(x, y)
        lats += [lat]
        lons += [lon]

        x, y = proj_area_def.lonlat2colrow(lon, lat)
        val = proj_data.values[y, x]
        vals += [val]

    return lats, lons, vals

def _calculate_composite_product_data(data, product_name):

    if product_name == 'day_land_cloud_fire':

        red = data['CMI_C06'].data
        green = data['CMI_C03'].data
        blue = data['CMI_C02'].data

        red_bounds = (0.0, 1.0) #in albedo (%)
        green_bounds = (0.0, 1.0) #in albedo (%)
        blue_bounds = (0.0, 1.0) #in albedo (%)

        red = ((red - red_bounds[0]) / (red_bounds[1] - red_bounds[0]))
        green = ((green - green_bounds[0]) / (green_bounds[1] - green_bounds[0]))
        blue = ((blue - blue_bounds[0]) / (blue_bounds[1] - blue_bounds[0]))

        pallete = None

        human_product_name = "Day Land Cloud/Fire"

    if product_name == 'day_cloud_phase':

        red = data['CMI_C13'].data
        green = data['CMI_C02'].data
        blue = data['CMI_C05'].data

        red = red - 273.15 #convert from kelvin to celsius

        red_bounds = (-53.5, 7.5) #in degrees C
        green_bounds = (0.0, 0.78) #in albedo (%)
        blue_bounds = (0.01, 0.59) #in albedo (%)

        #normalize
        red = ((red - red_bounds[1]) / (red_bounds[0] - red_bounds[1]))
        green = ((green - green_bounds[0]) / (green_bounds[1] - green_bounds[0]))
        blue = ((blue - blue_bounds[0]) / (blue_bounds[1] - blue_bounds[0]))

        pallete = None

        human_product_name = 'Day Cloud Phase'

    if product_name == 'nt_microphysics':

        c15 = data['CMI_C15'].data #12.4 micron
        c13 = data['CMI_C13'].data #10.4 micron
        c7 = data['CMI_C07'].data #3.9 micron

        red = c15 - c13
        green = c13 - c7
        blue = c13

        #red = red - 273.15
        #green = green - 273.15
        blue = blue - 273.15

        red_bounds = (-6.7, 2.6) #in degrees C
        green_bounds = (-3.1, 5.2) #in degrees C
        blue_bounds = (-29.6, 19.5) #in degrees C

        red = ((red - red_bounds[0]) / (red_bounds[1] - red_bounds[0]))
        green = ((green - green_bounds[0]) / (green_bounds[1] - green_bounds[0]))
        blue = ((blue - blue_bounds[0]) / (blue_bounds[1] - blue_bounds[0]))

        pallete = None

        human_product_name = "Nighttime Microphysics"

    if product_name == 'dust_rgb':

        c15 = data['CMI_C15'].data #12.4 micron
        c14 = data['CMI_C14'].data #11.2 micron
        c13 = data['CMI_C13'].data #10.4 micron
        c11 = data['CMI_C11'].data #8.4 micron

        red = c15 - c13
        green = c14 - c11
        blue = c13

        #red = red - 273.15
        # green = green - 273.15
        blue = blue - 273.15

        red_bounds = (-6.7, 2.6) #in degrees C
        green_bounds = (-0.5, 20.0) #in degrees C
        blue_bounds = (-11.95, 15.55) #in degrees C

        red = ((red - red_bounds[0]) / (red_bounds[1] - red_bounds[0]))
        green = ((green - green_bounds[0]) / (green_bounds[1] - green_bounds[0]))
        blue = ((blue - blue_bounds[0]) / (blue_bounds[1] - blue_bounds[0]))

        green = np.power(green, 1/2.5) #gamma correction

        pallete = None

        human_product_name = "Dust RGB"

    if product_name == 'true_color':

        c1 = data['CMI_C01'].data #12.4 micron
        c2 = data['CMI_C02'].data #11.2 micron
        c3 = data['CMI_C03'].data #10.4 micron

        red = c2
        veggie = c3
        blue = c1

        red = np.clip(red, 0, 1)
        veggie = np.clip(veggie, 0, 1)
        blue = np.clip(blue, 0, 1)

        gamma = 2.2
        red = np.power(red, 1/gamma)
        veggie = np.power(veggie, 1/gamma)
        blue = np.power(blue, 1/gamma)

        true_green = 0.45 * red + 0.1 * veggie + 0.45 * blue
        green = true_green

        pallete = None

        human_product_name = "True Color"

    if product_name == 'geocolor_day':

        c1 = data['CMI_C02'].data #12.4 micron
        c2 = data['CMI_C02'].data #11.2 micron
        c3 = data['CMI_C03'].data #10.4 micron

        red = c2
        green = c3
        blue = c1

        #red = red - 273.15
        # green = green - 273.15
        #blue = blue - 273.15

        # red_bounds = (-6.7, 2.6) #in degrees C
        # green_bounds = (-0.5, 20.0) #in degrees C
        # blue_bounds = (-11.95, 15.55) #in degrees C

        # red = ((red - red_bounds[0]) / (red_bounds[1] - red_bounds[0]))
        # green = ((green - green_bounds[0]) / (green_bounds[1] - green_bounds[0]))
        # blue = ((blue - blue_bounds[0]) / (blue_bounds[1] - blue_bounds[0]))

        true_green = 0.45 * red + 0.1 * green + 0.45 * blue #correct from "veggie" to true green
        green = true_green

        pallete = None

        human_product_name = "True Color"

    if product_name == 'airmass':

        c15 = data['CMI_C15'].data #12.4 micron
        c14 = data['CMI_C14'].data #11.2 micron
        c13 = data['CMI_C13'].data #10.4 micron
        c11 = data['CMI_C11'].data #8.4 micron

        red = c15 - c13
        green = c14 - c11
        blue = c13

        #red = red - 273.15
        # green = green - 273.15
        blue = blue - 273.15

        red_bounds = (-6.7, 2.6) #in degrees C
        green_bounds = (-0.5, 20.0) #in degrees C
        blue_bounds = (-11.95, 15.55) #in degrees C

        red = ((red - red_bounds[0]) / (red_bounds[1] - red_bounds[0]))
        green = ((green - green_bounds[0]) / (green_bounds[1] - green_bounds[0]))
        blue = ((blue - blue_bounds[0]) / (blue_bounds[1] - blue_bounds[0]))

        green = np.power(green, 1/2.5) #gamma correction

        pallete = None

        human_product_name = "Dust RGB"

 
    return red, green, blue, pallete, human_product_name

if __name__ == "__main__":

    parser = argparse.ArgumentParser(description=f'Plot GOES 16/17/18 Single-Band or Composite Data')
    parser.add_argument('--help-bands',
                        help='print details about available bands',
                        action='store_true')
    parser.add_argument('--help-composite',
                        help='print details about available composite products',
                        action='store_true')
    parser.add_argument('-b', '--band', 
                        help='plot band number # (1-16)', 
                        choices=range(1, 17), 
                        type=int,
                        default=None)
    parser.add_argument('-c', '--composite',
                        help='plot composite product (string)',
                        type=str,
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
                        help='specify a CSV file to read in points from. CSV file must have no header and be in format: lat,lon,marker,label',
                        type=str,
                        metavar='file_path',
                        dest='points_file',
                        default=None)
    parser.add_argument('-p', '--pallete',
                        help='specify color pallete to be used (single band only)',
                        type=str,
                        default=None)
    parser.add_argument('--save-as-kmz',
                        help='save as a georeferenced kmz rather than an image',
                        action='store_true',
                        dest='save_to_kmz')
    parser.add_argument('-cb', '--show-colorbar',
                        help='show a colorbar corresponding to pixel color (SB only)',
                        dest='show_colorbar',
                        action='store_true')
    parser.add_argument('-pv', '--pixel-value',
                        help='plot value of pixel (SB only)',
                        action='store_true',
                        dest='pixel_value')
    parser.add_argument('input_file_directory',  
                        help='directory to read input files from',
                        type=str)
    parser.add_argument('save_directory',  
                        help='directory to save output plots to',
                        type=str)

    args = parser.parse_args()
    print(args)

    ## Quality Checks/Sanitizing
    ##TODO: Add all CLI args that map to kwargs into a dict

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

    if args.composite:
        plot_composite = True
        product_desc = f"'{IMPLEMENTED_COMPOSITES[args.composite]}' composite"
    elif args.band:
        plot_composite = False
        product_desc = f"band {args.band} ({IMPLEMENTED_BAND_NAMES[str(args.band)]})"
    else:
        raise ValueError("A band or composite product to plot must be specified.")

    if args.pallete:
        pal = args.pallete

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

    kwargs = {}
    if args.save_to_kmz:
        kwargs.update({'save_to_kmz': True})

    if args.pixel_value:
        kwargs.update({'pixel_value': True})

    if args.show_colorbar:
        kwargs.update({'colorbar_visible': True})


    input_files = sorted(glob.glob(f'{input_dir}*.nat'))
    num_input_files = len(input_files)

    tot_files = 0

    with tqdm(miniters=0, total=len(input_files), desc=f'Plotting {product_desc}...', ascii=" ░▒▓█") as progress:
        for input_file_path in input_files:

            if plot_composite:
                plot_composite_goes(input_file_path, 
                                    save_dir, 
                                    points, 
                                    bbox, 
                                    args.composite, 
                                    **kwargs)
            else:
                plot_single_band_meteosat(input_file_path, 
                                          save_dir, 
                                          args.band, 
                                          points, 
                                          bbox, 
                                          pal, 
                                          **kwargs)
            progress.update()

    print("Done!")






