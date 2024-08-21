## Set to True if you want more towns plotted on the map
## (Default: False)
TOWN_SCALE_RANK = 5

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
POINT_STYLE = {'color': 'red', 'markersize': 6}
POINT_LABEL_VISIBLE = True
POINT_LABEL_STYLE = {'color': 'red', 'fontsize': 12}
DRAW_LABEL_ARROWS = True


## Set the positioning of the labels relative to the points being plotted
## (Default: X=0.2, Y=0.1)
Y_LABEL_OFFSET = -0.1
X_LABEL_OFFSET = 0.1

## Set the positioning of the pixel values relative to the upper-left corner
## of the pixel.
## (Default: X=5, Y=10)
X_PIX_VAL_OFFSET = -8
Y_PIX_VAL_OFFSET = -7











# 6/13/24 - Added cmasher (cmr) import for cmasher colormap support

#TODO:
#Add in ability to position points in CSV
#Add in ability to use settings file (see plot_plan_view.py)




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

from __internal_funcs import (plot_towns, draw_logo, plot_points,
                              define_hi_res_fig,
                              define_gearth_compat_fig, save_figs_to_kml)

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
    'nighttime_microphysics' : 'Nighttime Microphysics',
    'day_land_cloud_fire' : 'Day Land Cloud Fire',
    'day_land_cloud' : 'Day Land Cloud',
    'day_cloud_phase' : 'Day Cloud Phase',
    'day_snow_fog' : 'Day Snow/Fog',
    'true_color' : 'True Color',

}

IMPLEMENTED_BANDS = {
    '1': 'Blue (0.47 µm, 1 km)',
    '2': 'Red (0.64 µm, 0.5 km)',
    '3': 'Veggie (0.86 µm, 1 km)',
    '4': 'Cirrus (1.37 µm, 2 km)',
    '5': 'Snow/Ice (1.6 µm, 1 km)',
    '6': 'Cloud Particle Size (2.2 µm, 2 km)',
    '7': 'Shortwave IR (3.9 µm, 2 km)',
    '8': 'High-level Water Vapor (6.2 µm, 2 km)',
    '9': 'Mid-level Water Vapor (6.9 µm, 2 km)',
    '10': 'Low-level Water Vapor (7.3 µm, 2 km)',
    '11': 'Cloud-Top Phase (8.4 µm, 2 km)',
    '12': 'Ozone (9.6 µm, 2 km)',
    '13': '"Clean" Longwave IR (10.3 µm, 2 km)',
    '14': 'Longwave IR (11.2 µm, 2 km)',
    '15': '"Dirty" Longwave IR (12.3 µm, 2 km)',
    '16': '"CO2" Longwave IR" (13.3 µm, 2 km)',
}

def plot_single_band_goes(file_path: str, 
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
        data = xr.open_dataset(file_path, engine="netcdf4")
    except Exception as e:
        raise e

    scan_start = datetime.strptime(data.time_coverage_start, '%Y-%m-%dT%H:%M:%S.%fZ')
    scan_end = datetime.strptime(data.time_coverage_end, '%Y-%m-%dT%H:%M:%S.%fZ')
    file_created = datetime.strptime(data.date_created, '%Y-%m-%dT%H:%M:%S.%fZ')
    orbital_slot = data.orbital_slot #GOES-East, GOES-West, GOES-Test, etc.
    sat_id = data.platform_ID #G18, G17, G16, etc.

    proj_info = data.variables['goes_imager_projection']
    # print(data.variables)
    # radiance = data.variables['Rad'][:]

    # fk1 = data.variables['planck_fk1'][:]
    # fk2 = data.variables['planck_fk2'][:]
    # bc1 = data.variables['planck_bc1'][:]
    # bc2 = data.variables['planck_bc2'][:]
    
    #BAND SELECTION

    sel_band_str = 'CMI_C' + str(band).zfill(2)
    sel_band = data[sel_band_str].data
    sel_band_name = data['band_id_C' + str(band).zfill(2)].long_name

    x = data.x
    y = data.y
    
    pix_lats, pix_lons = _calculate_pixel_lat_lon(x, y, proj_info)

    data = data.metpy.parse_cf('CMI_C02')
    geog_data = data.metpy.cartopy_crs
    x = data.x
    y = data.y

    #x, y = np.meshgrid(x, y)

    #DATA CORRECTIONS

    #Bands 1-6 are output in "Reflectance Factor", from 0-1 and can be clip corrected
    #Bands 7-16 are output in "Brightness Temperature", and cannot be clip or gamma corrected
    if band >= 7: 
        correct_clip = False
        correct_gamma = False

        img_vmin = 180
        img_vmax = 330

        cbar_label = "Radiance Temperature (K)"

    else:
        correct_clip = True
        correct_gamma = True

        img_vmin = 0
        img_vmax = 1

        cbar_label = "Apparent Brightness"

    if correct_clip:
        sel_band = np.clip(sel_band, 0, 1)

    if correct_gamma:
        sel_band = np.power(sel_band, 1/DEF_GAMMA)

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
    
    img = ax.imshow(sel_band, 
                    origin='upper',
                    extent=(x.min(), x.max(), y.min(), y.max()),
                    transform=geog_data,
                    interpolation='none',
                    vmin = img_vmin,
                    vmax = img_vmax,
                    cmap = pallete)

    
    # ir = [fk2 / (np.log((fk1/radiance) + 1)) - bc1] / bc2
    
    # ir = ir[0,:,:]

    if kwargs.get('pixel_value'):

        pix_lats[pix_lats > bbox[0] + 1] = np.nan
        pix_lats[pix_lats < bbox[1] - 1] = np.nan

        pix_lons[pix_lons > bbox[2] + 1] = np.nan
        pix_lons[pix_lons < bbox[3] - 1] = np.nan
        
        for lat_row, lon_row, val_row in zip(pix_lats, pix_lons, sel_band[:]):
            for lat, lon, val in zip(lat_row, lon_row, val_row):
                if not (np.isnan(lat) or np.isnan(lon)):
                    ax.annotate(str(round(val, 1)), 
                                (lon, lat),
                                xytext=(X_PIX_VAL_OFFSET,Y_PIX_VAL_OFFSET),
                                textcoords='offset points',
                                horizontalalignment='center',
                                verticalalignment='center', 
                                color='black',
                                clip_box=ax.bbox,
                                fontsize=6,
                                transform=crs.PlateCarree(),
                                annotation_clip=False, 
                                zorder=30)

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

        file_name = sat_id + "_" + sel_band_str + "_" + scan_end.strftime('%Y%m%d_%H%M%S%Z')
        layer_name = sat_id.replace("G", "GOES-") + " image at " + scan_end.strftime('%d %B %Y %H:%M UTC ')
        layer_desc = IMPLEMENTED_BANDS[str(band)]

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

        title_info = orbital_slot.replace("-Test", "") + " (" + sat_id.replace("G", "GOES-") + ")\n" + IMPLEMENTED_BANDS[str(band)]
        
        ax.set_title(title_info, loc='left', fontweight='bold', fontsize=15)
        ax.set_title('{}'.format(scan_end.strftime('%d %B %Y %H:%M UTC ')), loc='right')
        
        file_name = sat_id + "_" + sel_band_str + "_" + scan_end.strftime('%Y%m%d_%H%M%S%Z')
        
        if not os.path.exists(save_dir):
            os.makedirs(save_dir)

        dest_path = os.path.join(save_dir, file_name + ".png")

        fig.savefig(dest_path, bbox_inches="tight", dpi=FILE_DPI)
        plt.close(fig)

def plot_composite_goes(file_path: str, 
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

def _calculate_pixel_lat_lon(proj_x, proj_y, proj_info):

    #print(proj_info)
    proj_x, proj_y = np.meshgrid(proj_x, proj_y)

    #Thanks to Mikhail Krotkin

    r_eq = 6378137 #Equitorial radius in m
    r_pol = 6356752.31414 #Polar radius in m
    
    ### Longitude of Origin ###
    lon_projection = proj_info.attrs['longitude_of_projection_origin']
    l_0 = lon_projection * (np.pi/180) #Longitude in Radians
    
    H = proj_info.attrs['perspective_point_height']+proj_info.attrs['semi_major_axis']
    
    ### Formula to Determine Lat/Longs from Satellite proj_x/proj_y coordinates that are given in Radians ###
    a = np.sin(proj_x)**2 + (np.cos(proj_x)**2 * (np.cos(proj_y)**2 + (r_eq**2 / r_pol**2) * np.sin(proj_y)**2))
    b = -2 * H * np.cos(proj_x) * np.cos(proj_y)
    c = (H**2) - (r_eq**2)
    
    r_s = (-b - np.sqrt(b**2 - 4*a*c))/(2*a)
    
    s_x = r_s * np.cos(proj_x) * np.cos(proj_y)
    s_y = -r_s * np.sin(proj_x)
    s_z = r_s * np.cos(proj_x) * np.sin(proj_y)
    
    lats = np.arctan((r_eq**2 / r_pol**2) * (s_z / np.sqrt((H-s_x)**2 +s_y**2))) * (180/np.pi)
    lons = (l_0 - np.arctan(s_y / (H-s_x))) * (180/np.pi)

    #print('{} N, {} W'.format(lats[318,1849],abs(lons[318,1849])))

    return lats, lons

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

    if product_name == 'day_land_cloud':

        red = data['CMI_C05'].data #1.6 micron
        green = data['CMI_C03'].data #0.86 micron
        blue = data['CMI_C02'].data #0.64 micron

        red_bounds = (0.0, 0.975) #in albedo (%)
        green_bounds = (0.0, 1.086) #in albedo (%)
        blue_bounds = (0.0, 1.0) #in albedo (%)

        red = ((red - red_bounds[0]) / (red_bounds[1] - red_bounds[0]))
        green = ((green - green_bounds[0]) / (green_bounds[1] - green_bounds[0]))
        blue = ((blue - blue_bounds[0]) / (blue_bounds[1] - blue_bounds[0]))

        pallete = None

        human_product_name = "Day Land Cloud/Fire"

    if product_name == 'day_snow_fog':

        c13 = data['CMI_C13'].data #10.3 micron
        c7 = data['CMI_C07'].data #3.9 micron
        c5 = data['CMI_C05'].data #1.6 micron
        c3 = data['CMI_C03'].data #0.86 micron

        red = c3
        green = c5
        blue = c7 - c13

        # blue = blue - 273.15 #convert from kelvin to celsius

        red_bounds = (0.0, 1.0) #in albedo (%)
        green_bounds = (0.0, 0.7) #in albedo (%)
        blue_bounds = (0.0, 30.0) #in degrees C

        red = ((red - red_bounds[0]) / (red_bounds[1] - red_bounds[0]))
        green = ((green - green_bounds[0]) / (green_bounds[1] - green_bounds[0]))
        blue = ((blue - blue_bounds[1]) / (blue_bounds[0] - blue_bounds[1]))

        red = np.power(red, 1/1.7)
        green = np.power(green, 1/1.7)
        blue = np.power(blue, 1/1.7)

        pallete = None

        human_product_name = "Day Snow/Fog"

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

    if product_name == 'nighttime_microphysics':

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
        product_desc = f"band {args.band} ({IMPLEMENTED_BANDS[str(args.band)]})"
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


    input_files = sorted(glob.glob(f'{input_dir}*.nc'))
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
                plot_single_band_goes(input_file_path, 
                                      save_dir, 
                                      args.band, 
                                      points, 
                                      bbox, 
                                      pal, 
                                      **kwargs)
            progress.update()

    print("Done!")






