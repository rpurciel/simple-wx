## When plotting cumulative flashes, how many minutes to display historical flashes for
## (Default: 20)
DEFAULT_HISTORY_MINUTES = 20

## Set the DPI of the saved output file
FILE_DPI = 300

'''
Used to set the style of the flashes shown on the map.
New (within current 20s frame) flashes use the NEW_FLASH_STYLE,
and historical flashes use the OLD_FLASH_STYLE.
'''
NEW_FLASH_STYLE = {'markersize': 8, 'marker': 'x', 'color': 'red'}
OLD_FLASH_STYLE = {'markersize': 4, 'marker': 'x', 'color': 'darkgray'}











'''


            STATIC CODE BELOW HERE


'''





import os
import glob
from datetime import datetime, timedelta
import argparse
import csv
import warnings
import json
warnings.simplefilter("ignore")

import pandas as pd
import xarray as xr
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.patheffects as PathEffects
import cartopy
import cartopy.crs as crs
import cartopy.feature as cfeature
import cartopy.io.shapereader as shpreader
from cartopy.feature import NaturalEarthFeature
from tqdm.auto import tqdm
from adjustText import adjust_text 

from __internal_funcs import (plot_towns, draw_logo, plot_points,
                              define_hi_res_fig,
                              define_gearth_compat_fig, save_figs_to_kml)

np.seterr(divide = 'ignore', invalid='ignore')
mpl.use('agg')

def plot_single_time_glm(file_path: str, 
                         save_dir: str, 
                         points: list[tuple[float, float, str, str], ...], 
                         bbox: list[float, float, float, float], 
                         **kwargs) -> pd.DataFrame: 

    """
    This function will take a GLM data file in NetCDF format, and plot any flashes
    it finds in the file, along with any other specified points. Will only plot flashes
    found in this specific file (i.e. "time").

    Inputs
        - file_path: str:
                A path to a GLM NetCDF data file.
        - save_dir, str: 
                The desired destination to save the plot.
        - point, list[tuple[float, float, str, str]]:
                A list of points to add to the plot. Each point has the format
                '(lat [float], lon [float], marker [str], label [str])'
        - bbox, list[float, float, float, float]:
                A list of floats defining the bounding box of the image. The function
                pulls the boundaries in the following order: N, S, E, W.
    
    Returns
        A Pandas DataFrame, with the lats, lons, and times of all the flashes detected 
        in this file.
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
    
    flash_lats = data.variables['flash_lat'][:]
    flash_lons = data.variables['flash_lon'][:]

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
    
    flashes = ax.plot(flash_lats, 
                      flash_lons,
                      linestyle='None',
                      transform=crs.PlateCarree(), 
                      zorder=15)

    plt.setp(flashes, **NEW_FLASH_STYLE)

    if points:
        plot_points(plt, ax,
                    points,
                    **kwargs)

    if kwargs.get('save_to_kmz'):

        file_name = sat_id + "_GLM_" + scan_end.strftime('%Y%m%d_%H%M%S%Z')
        layer_name = f'{orbital_slot.replace("-Test", "")} ({sat_id.replace("G", "GOES-")})\nGLM Detected Lightning Flashes'
        layer_desc = f'Starting At {scan_start.strftime("%d %B %Y %H:%M:%S UTC")}\nThrough {scan_end.strftime("%d %B %Y %H:%M:%S UTC")}'

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
                  scale_rank=kwargs.pop('plot_towns_scale_rank', 5))
        draw_logo(ax)

        ax.set_title(f'{orbital_slot.replace("-Test", "")} ({sat_id.replace("G", "GOES-")})\nDetected Lightning Flashes', loc='left', fontweight='bold', fontsize=15)
        ax.set_title(f'Starting At {scan_start.strftime("%d %B %Y %H:%M:%S UTC")}\nThrough {scan_end.strftime("%d %B %Y %H:%M:%S UTC")}', loc='right')
        
        file_name = sat_id + "_GLM_" + scan_end.strftime('%Y%m%d_%H%M%S%Z')
        
        if not os.path.exists(save_dir):
            os.makedirs(save_dir)

        dest_path = os.path.join(save_dir, file_name + ".png")

        fig.savefig(dest_path, bbox_inches="tight", dpi=FILE_DPI)
        plt.close(fig)

def plot_cumulative_glm(file_path: str, 
                        save_dir: str, 
                        points: list[tuple[float, float, str], ...], 
                        bbox: list[float, float, float, float], 
                        hist_flashes_df: pd.DataFrame,
                        cum_minutes: int, 
                        **kwargs) -> pd.DataFrame:

    """
    This function will take a GLM data file in NetCDF format, and plot any flashes
    it finds in the file, along with any other specified points. This takes an 
    Inputs
        - file_path: str:
                A path to a GLM NetCDF data file.
        - save_dir, str: 
                The desired destination to save the plot.
        - point, list[tuple[float, float, str, str]]:
                A list of points to add to the plot. Each point has the format
                '(lat [float], lon [float], marker [str], label [str])'
        - bbox, list[float, float, float, float]:
                A list of floats defining the bounding box of the image. The function
                pulls the boundaries in the following order: N, S, E, W.
        - hist_flashes_df, pandas.DataFrame:
                A Pandas DataFrame, with the lats, lons, and times of all historical
                flashes to plot.
        - cum_minutes, int:
                The number of minutes to keep a rolling window of historical flashes for.
                Will remove historical flashes if t(flash) <= t(this_file) - cum_minutes 
    
    Returns
        A Pandas DataFrame, with the lats, lons, and times of all the historical flashes
        within the specified window.
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

    if 'start_of_period' in kwargs:
        start_of_period = kwargs.get('start_of_period')

    if start_of_period:
        period_st = scan_start
    else:
        hist_flash_filtered = hist_flashes_df[hist_flashes_df['flash_time'] > scan_start - timedelta(minutes=cum_minutes)]
        period_st = min(hist_flash_filtered['flash_time'])

    flash_lats = data.variables['flash_lat'].values
    flash_lons = data.variables['flash_lon'].values

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
    
    flashes = ax.plot(flash_lons, 
                      flash_lats,
                      linestyle='None',
                      transform=crs.PlateCarree(), 
                      zorder=15)



    mpl.artist.setp(flashes, **NEW_FLASH_STYLE)

    #create flash dataframe
    this_frame_data = {'flash_lat': flash_lats, 'flash_lon': flash_lons, 'flash_time': scan_start}
    this_frame_df = pd.DataFrame(data=this_frame_data)

    lon_bounds = (bbox[0], bbox[1])
    lat_bounds = (bbox[2], bbox[3])

    lat_max_E = max(lat_bounds)
    lat_min_W = min(lat_bounds)
    lon_max_N = max(lon_bounds)
    lon_min_S = min(lon_bounds)

    if kwargs.get('timestamp_flashes'):  
        for flash in this_frame_df.itertuples(name=None, index=False):
            flash_lat = float(flash[1])
            flash_lon = float(flash[0])

            if (flash_lat >= lat_min_W) and (flash_lat <= lat_max_E) and (flash_lon >= lon_min_S) and (flash_lon <= lon_max_N):
                pt_label = ax.annotate(scan_start,
                                       xy=(flash_lat, flash_lon),
                                       xytext=(0, -50), 
                                       xycoords='data',
                                       textcoords='offset pixels',
                                       horizontalalignment="center",
                                       verticalalignment="top",
                                       fontsize=8,
                                       color="red",
                                       transform=crs.PlateCarree(),
                                       annotation_clip=True, 
                                       zorder=14)

    if not start_of_period:
        hist_flashes = ax.plot(hist_flash_filtered['flash_lon'], 
                               hist_flash_filtered['flash_lat'],
                               linestyle='none',
                               transform=crs.PlateCarree(), 
                               zorder=10)

        if kwargs.get('timestamp_flashes'):
            for flash in hist_flash_filtered.itertuples(name=None, index=False):
                flash_lat = float(flash[1])
                flash_lon = float(flash[0])

                if (flash_lat >= lat_min_W) and (flash_lat <= lat_max_E) and (flash_lon >= lon_min_S) and (flash_lon <= lon_max_N):
                    pt_label = ax.annotate(flash[2].strftime("%Y-%m-%d %H:%M:%S"), 
                                           xy=(flash_lat, flash_lon),
                                           xytext=(0, -30), 
                                           xycoords='data',
                                           textcoords='offset pixels',
                                           horizontalalignment="center",
                                           verticalalignment="top",
                                           fontsize=6,
                                           fontstyle='italic',
                                           color="gray",
                                           transform=crs.PlateCarree(),
                                           annotation_clip=True, 
                                           zorder=9)


        mpl.artist.setp(hist_flashes, **OLD_FLASH_STYLE)
        hist_flashes_df = pd.concat([hist_flashes_df, this_frame_df], ignore_index=True)
    else:
        hist_flashes_df = this_frame_df

    if points:
        plot_points(plt, ax,
                    points,
                    **kwargs)

    if kwargs.get('save_to_kmz'):

        file_name = sat_id + "_GLM_" + scan_end.strftime('%Y%m%d_%H%M%S%Z')
        layer_name = f'{orbital_slot.replace("-Test", "")} ({sat_id.replace("G", "GOES-")})\nCumulative Detected Lightning Flashes'
        layer_desc = f'Starting At {period_st.strftime("%d %B %Y %H:%M:%S UTC")}\nThrough {scan_end.strftime("%d %B %Y %H:%M:%S UTC")}\n{cum_minutes} Min. Rolling Period'


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
                  scale_rank=kwargs.pop('plot_towns_scale_rank', 5))
        draw_logo(ax)

        ax.set_title(f'{orbital_slot.replace("-Test", "")} ({sat_id.replace("G", "GOES-")})\nCumulative Detected Lightning "Flashes"', loc='left', fontweight='bold', fontsize=15)
        ax.set_title(f'Starting At {period_st.strftime("%d %B %Y %H:%M:%S UTC")}\nThrough {scan_end.strftime("%d %B %Y %H:%M:%S UTC")}', loc='right')
        ax.set_title(f'{cum_minutes} Min. Rolling Period', loc='center', fontsize=8, fontstyle='italic')
        
        file_name = sat_id + "_GLM_" + period_st.strftime('%Y%m%d_%H%M%S%Z') + "_to" + scan_end.strftime('%Y%m%d_%H%M%S%Z')
        
        if not os.path.exists(save_dir):
            os.makedirs(save_dir)

        dest_path = os.path.join(save_dir, file_name + ".png")

        fig.savefig(dest_path, 
                    bbox_inches="tight", 
                    dpi=FILE_DPI)
        plt.close(fig)

    return hist_flashes_df

if __name__ == "__main__":

    parser = argparse.ArgumentParser(description=f'Plot GOES 16/17/18 GLM Data, either for single time periods or cumulative over all input data.')
    parser.add_argument('--single', 
                        help='plot flashes for single-time periods', 
                        action='store_true',
                        default=False)
    parser.add_argument('--cumulative', 
                        help=f'plot cumulative flashes over a rolling time period, optionally specified (default: {DEFAULT_HISTORY_MINUTES} min)', 
                        type=int,
                        nargs='?',
                        metavar='min',
                        const=DEFAULT_HISTORY_MINUTES,
                        default=DEFAULT_HISTORY_MINUTES)
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
    parser.add_argument('--timestamp-flashes', 
                        help='save all detected flashes in a kmz', 
                        action='store_true',
                        dest="timestamp_flashes",
                        default=False)
    parser.add_argument('--save-as-kmz', 
                        help='save all detected flashes in a kmz', 
                        action='store_true',
                        dest='save_to_kmz',
                        default=False)
    parser.add_argument('input_file_directory',  
                        help='directory to read input files from',
                        type=str)
    parser.add_argument('save_directory',  
                        help='directory to save output plots to',
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

    if args.single:
        plot_cumulative = False
        product_desc = f"GLM flashes"
    elif args.cumulative:
        plot_cumulative = True
        hist_data = pd.DataFrame(data=[], columns=['flash_lat', 'flash_lon', 'flash_time'])
        product_desc = "GLM flashes (cumulative)"
        cum_period = args.cumulative
    else:
        raise ValueError("An option for plotting single-time or cumulative flashes must be specified.")

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

    if args.settings_file:
        if ('default' in args.settings_file) and ('/' not in args.settings_file):
            pass
        try:
            with open(args.settings_file, newline='') as settingsjson:
                user_settings = json.load(settingsjson)
        except Exception as err:
            raise err
    else:
        user_settings = dict()

    if args.save_to_kmz:
        user_settings.update({'save_to_kmz': True})
    if args.timestamp_flashes:
        user_settings.update({'timestamp_flashes': True})

    input_files = sorted(glob.glob(f'{input_dir}*.nc'))
    num_input_files = len(input_files)

    plot_times = []
    start_of_period = True

    with tqdm(miniters=0, total=len(input_files), desc=f'Plotting {product_desc}...', ascii=" ░▒▓█") as progress:
        for input_file_path in input_files:
            if plot_cumulative:
                hist_data = plot_cumulative_glm(input_file_path, 
                                                save_dir, 
                                                points, 
                                                bbox, 
                                                hist_data, 
                                                cum_period,
                                                start_of_period=start_of_period,
                                                **user_settings)
                if start_of_period:
                    start_of_period = False
            else:
                plot_single_time_glm(input_file_path, 
                                     save_dir, 
                                     points, 
                                     bbox,
                                     **user_settings)
            progress.update()

    if args.save_to_kmz:
        hist_data.to_csv(os.path.join(save_dir, "GLM_all_flashes.csv"), index=None)

    print("Done!")






