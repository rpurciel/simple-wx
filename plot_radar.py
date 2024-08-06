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

DEFAULT_NEXRAD_PRODUCT = 'reflectivity'
DEFAULT_NEXRAD_SCAN_ANGLE = 0.5
## Set the DPI of the saved output file
## (Default: 300)
FILE_DPI = 300

TOWN_SCALE_RANK = 5

'''
Used to set the style of any user input points shown on the map.
POINT_STYLE is used for the point(s), and POINT_LABEL_STYLE is used for any
point label(s) given. No labels are shown if POINT_LABEL_VISIBLE is False.
'''
POINT_STYLE = {'color': 'black', 'markersize': 10, 'markeredgecolor': 'black'}
POINT_LABEL_VISIBLE = True
DRAW_LABEL_ARROWS = True
POINT_LABEL_STYLE = {'color': 'black', 'fontsize': 14, 'fontweight': 'bold'}

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
from matplotlib.colors import from_levels_and_colors
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
import pyart

from __internal_funcs import (plot_towns, draw_logo, 
                              plot_points, define_hi_res_fig,
                              define_gearth_compat_fig, save_figs_to_kml)

DBZ_LEVELS = np.arange(5., 75., 5.)

DBZ_COLORS = np.array([[4,233,231],
                      [1,159,244], 
                      [3,0,244],
                      [2,253,2], 
                      [1,197,1],
                      [0,142,0], 
                      [253,248,2],
                      [229,188,0], 
                      [253,149,0],
                      [253,0,0], 
                      [212,0,0],
                      [188,0,0],
                      [248,0,253],
                      [152,84,198]], 
                      np.float32) / 255.0

DBZ_CMAP, DBZ_NORM = from_levels_and_colors(DBZ_LEVELS, 
                                            DBZ_COLORS,
                                            extend="max")

np.seterr(divide = 'ignore', invalid='ignore')
mpl.use('agg')

def plot_nexrad_l2(file_path: str,
                   save_dir: str,
                   points: list[tuple[float, float, str, str], ...],
                   bbox: list[float, float, float, float],
                   product: str = DEFAULT_NEXRAD_PRODUCT,
                   sweep_idx: float = None,
                   sweep_angle: float = DEFAULT_NEXRAD_SCAN_ANGLE,
                   **kwargs) -> None:

    try:
        data = pyart.io.read_nexrad_archive(file_path)
        radar_display = pyart.graph.RadarMapDisplay(data)
    except Exception as e:
        raise e

    valid_products = [prod for prod in data.fields]
    if product not in valid_products:
        raise ValueError(f'Specified product not found, must be one of {valid_products}')

    scan_time = datetime.strptime(data.time['units'][14:], "%Y-%m-%dT%H:%M:%SZ")
    radar_id = data.metadata['instrument_name']
    sweep_angles = data.fixed_angle['data']
    radar_lat = data.latitude['data'][0]
    radar_lon = data.longitude['data'][0]
    print(radar_lat, radar_lon)
    
    if sweep_angle:
        sel_sweep_idx = (np.abs(sweep_angles-sweep_angle)).argmin()
        sel_sweep_ang = sweep_angles[sel_sweep_idx]
    else:
        sel_sweep_idx = sweep_idx
        sel_sweep_ang = sweep_angles[sel_sweep_idx]

    if kwargs.get('save_to_kmz'):  
        proj = crs.PlateCarree(central_longitude=radar_lon)

        fig, ax, cbfig, cbax = define_gearth_compat_fig((bbox[2], bbox[3]),
                                                        (bbox[0], bbox[1]),
                                                        projection=proj,
                                                        draw_earth_features=True,
                                                        dpi_pixels=2048)

        cbar_opts = {'rotation': -90, 'color': 'k', 'labelpad': 20}

    else:
        fig, ax = define_hi_res_fig((bbox[2], bbox[3]),
                                    (bbox[0], bbox[1]))

        cbfig = fig
        cbax = None

        cbar_opts = {}

    radar_display.plot_ppi_map(product, 
                               sel_sweep_idx, 
                               colorbar_flag=False, 
                               title_flag=False,
                               # projection=proj,
                               ticks=DBZ_LEVELS, 
                               cmap=DBZ_CMAP, 
                               norm=DBZ_NORM,
                               add_grid_lines=False,
                               embellish=False,
                               min_lat=bbox[1],
                               max_lat=bbox[0],
                               min_lon=bbox[3],
                               max_lon=bbox[2],
                               vmin=-12, 
                               vmax=64, 
                               ax=ax)

    if kwargs.get('show_colorbar'):
        cb = cbfig.colorbar(mpl.cm.ScalarMappable(norm=DBZ_NORM, cmap=DBZ_CMAP),
                            ax=ax,
                            cax=cbax,
                            orientation = "vertical", 
                            pad=.05,
                            shrink=0.7,
                            use_gridspec=False)
        cb.set_label('Reflectivity (dBZ)', 
                     size='x-large', 
                     **cbar_opts)

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

        file_name = file_name = f"{radar_id}_{str(round(sel_sweep_ang, 1)).replace('.','_')}degScan_{scan_time.strftime('%Y%m%d_%H%M%S%Z')}"
        layer_name = f"{radar_id} Radar {product.title()}"
        layer_desc = f"{round(sel_sweep_ang, 1)}° Beam Angle\n{scan_time.strftime('%d %B %Y %H:%M:%S UTC')}"

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
        
        ax.set_title(f"{radar_id} Radar {product.title()}", 
                     loc='left', 
                     fontweight='bold', 
                     fontsize=15)

        ax.set_title(f"{round(sel_sweep_ang, 1)}° Beam Angle\n{scan_time.strftime('%d %B %Y %H:%M:%S UTC')}", 
                     loc='right')
        
        file_name = f"{radar_id}_{str(round(sel_sweep_ang, 1)).replace('.','_')}degScan_{scan_time.strftime('%Y%m%d_%H%M%S%Z')}"
        
        if not os.path.exists(save_dir):
            os.makedirs(save_dir)

        dest_path = os.path.join(save_dir, file_name + ".png")

        fig.savefig(dest_path, bbox_inches="tight", dpi=FILE_DPI)
        plt.close(fig)

def plot_mrms(file_path: str,
                save_dir: str,
                points: list[tuple[float, float, str, str], ...],
                bbox: list[float, float, float, float],
                **kwargs) -> None:

    try:
        data = xr.open_dataset(file_path, engine="netcdf4")
    except Exception as e:
        raise e

    lat = data.variables['latitude'][:]
    lon = data.variables['longitude'][:]
    reflectivity = data.variables['unknown'][:]
    scan_time = pd.Timestamp(data.time.values).to_pydatetime()

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

    ax.contourf(lon, lat, 
                reflectivity, 
                transform = crs.PlateCarree(), 
                levels=DBZ_LEVELS, 
                cmap=DBZ_CMAP, 
                norm=DBZ_NORM, 
                extend='max')

    if kwargs.get('show_colorbar'):
        cb = cbfig.colorbar(ax=ax,
                            cax=cbax,
                            orientation = "horizontal", 
                            pad=.05,
                            shrink=0.7,
                            use_gridspec=False)
        cb.set_label('Reflectivity (dBZ)', size='x-large', **cbar_opts)

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

if __name__ == "__main__":

    parser = argparse.ArgumentParser(description=f'Plot NEXRAD archive data or MRMS (Multi-Radar Multi-Sensor) data.')
    parser.add_argument('--nexrad-product',
                        help=f'specify a NEXRAD product to be plotted (default: {DEFAULT_NEXRAD_PRODUCT})',
                        type=str,
                        metavar='product',
                        dest='rad_prod',
                        default=None)
    parser.add_argument('--nexrad-scan-angle',
                        help=f'specify a scan angle to be used when plotting NEXRAD radar (default: {DEFAULT_NEXRAD_SCAN_ANGLE})',
                        type=str,
                        metavar='angle',
                        dest='rad_sa',
                        default=None)
    parser.add_argument('--nexrad-scan-index',
                        help='specify scan angle (via index) to be used when plotting NEXRAD radar',
                        type=str,
                        metavar='index',
                        dest='rad_sa_idx',
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
                        help='specify a JSON file to read in plot settings from',
                        type=str,
                        metavar='file_path',
                        dest='settings_file',
                        default=None)
    parser.add_argument('--save-as-kmz',
                        help='save as a georeferenced kmz rather than an image',
                        action='store_true',
                        dest='save_to_kmz')
    parser.add_argument('--remove-colorbar',
                        help='prevent showing a colorbar corresponding to reflectivity color',
                        dest='remove_colorbar',
                        action='store_true')
    parser.add_argument('input_file_directory',  
                        help='directory to read input files from',
                        type=str)
    parser.add_argument('save_directory',  
                        help='directory to save sounding(s) to',
                        type=str)

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

    if args.rad_prod:
        product = args.rad_prod
    else:
        product = DEFAULT_NEXRAD_PRODUCT

    if args.rad_sa:
        scan_angle = args.rad_sa
    else:
        scan_angle = DEFAULT_NEXRAD_SCAN_ANGLE

    if args.rad_sa_idx:
        scan_idx = args.rad_sa_idx
    else:
        scan_idx = None

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

    if args.remove_colorbar:
        user_settings.update({'show_colorbar': False})
    else:
        user_settings.update({'show_colorbar': True})

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

    input_files = sorted(glob.glob(f'{input_dir}*V06'))
    data_is_mrms = False

    if len(input_files) == 0:
        print("No NEXRAD single radar files detected, searching for MRMS format files...")
        input_files = sorted(glob.glob(f'{input_dir}*.grib2'))
        if len(input_files) == 0:
            raise ValueError("No input files found in NEXRAD archive or MRMS format. Try reviewing input parameters")
        else:
            print("MRMS format files found, assuming they are MRMS data files.")
            data_is_mrms = True

    num_input_files = len(input_files)
    num_plot_steps = num_input_files

    with tqdm(miniters=0, total=num_plot_steps, desc=f"Plotting {'MRMS' if data_is_mrms else 'NEXRAD radar'} data...", ascii=" ░▒▓█") as progress:
        for input_file_path in input_files:
            if data_is_mrms:
                plot_mrms(input_file_path,
                          save_dir,
                          points,
                          bbox,
                          **user_settings)
            else:
                plot_nexrad_l2(input_file_path,
                               save_dir,
                               points,
                               bbox,
                               product,
                               scan_idx,
                               scan_angle,
                               **user_settings)
            progress.update()

    print("Done!")






