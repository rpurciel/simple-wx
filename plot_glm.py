## When plotting cumulative flashes, how many minutes to display historical flashes for
## (Default: 20)
DEFAULT_HISTORY_MINUTES = 20

## Set the DPI of the saved output file
FILE_DPI = 300

'''
Used to set the style of the flashes shown on the map.
New (within current frame) flashes use the NEW_FLASH_STYLE,
and historical flashes use the OLD_FLASH_STYLE.
'''
NEW_FLASH_STYLE = {'markersize': 8, 'marker': 'x', 'color': 'red'}
OLD_FLASH_STYLE = {'markersize': 4, 'marker': 'x', 'color': 'darkgray'}

'''
Used to set the style of any user input points shown on the map.
POINT_STYLE is used for the point(s), and POINT_LABEL_STYLE is used for any
point label(s) given. No labels are shown if POINT_LABEL_VISIBLE is False.
'''
POINT_STYLE = {'color': 'black', 'markersize': 8, 'marker': 'o'}
POINT_LABEL_VISIBLE = True
POINT_LABEL_STYLE = {'color': 'black', 'fontsize': 10, 'fontweight': 'bold'}











'''


            STATIC CODE BELOW HERE


'''





import os
import glob
from datetime import datetime, timedelta
import argparse
import csv

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

from __internal_funcs import plot_towns, draw_logo

np.seterr(divide = 'ignore', invalid='ignore')
mpl.use('agg')

def plot_single_time_glm(file_path: str, 
                         save_dir: str, 
                         points: list[tuple[float, float, str], ...], 
                         bbox: list[float, float, float, float], 
                         **kwargs) -> pd.DataFrame: 

    '''
    A function to 
    '''

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

    fig = plt.figure(figsize=(22,16))
    ax = plt.axes(projection = crs.PlateCarree())
    
    bbox_wesn = [bbox[3], bbox[2], bbox[1], bbox[0]]
    ax.set_extent(bbox_wesn, crs=crs.PlateCarree())

    plot_logo(ax)
    
    states = NaturalEarthFeature(category="cultural", scale="50m",
                                          facecolor="none",
                                          name="admin_1_states_provinces")
    ax.add_feature(states, linewidth=1.0, edgecolor="black")
    ax.coastlines('50m', linewidth=1.5)
    ax.add_feature(cartopy.feature.LAKES.with_scale('10m'), linestyle='-', linewidth=0.5, alpha=1,edgecolor='blue',facecolor='none')
    ax.add_feature(cfeature.BORDERS, linewidth=1.5)
    plot_towns(ax, (bbox_wesn[0], bbox_wesn[1]), (bbox_wesn[2], bbox_wesn[3]), scale_rank=7)
    
    flashes = ax.plot(flash_lats, 
                      flash_lons,
                      linestyle='None',
                      transform=crs.PlateCarree(), 
                      zorder=15)

    plt.setp(flashes, **NEW_FLASH_STYLE)

    if points:
        input_pts = []
        input_pt_labels = []
        for point in points:
            x_axis = point[0]
            y_axis = point[1]
            label = point[2]
            
            pt = ax.plot([y_axis],[x_axis])
            input_pts += [pt]
           
            if POINT_LABEL_VISIBLE:
                pt_label = ax.annotate(label, 
                                      (y_axis, x_axis), 
                                      horizontalalignment='center', 
                                      verticalalignment='top', 
                                      transform=crs.PlateCarree(), 
                                      annotation_clip=True, 
                                      zorder=30)
                input_pt_labels = [pt_label]

        plt.setp(input_pts, **POINT_STYLE)
        plt.setp(input_pt_labels, **POINT_LABEL_STYLE)
        adjust_text(input_pt_labels)

    plt.title(f'{orbital_slot.replace("-Test", "")} ({sat_id.replace("G", "GOES-")})\nDetected Lightning "Flashes"', loc='left', fontweight='bold', fontsize=15)
    plt.title(f'Starting At {scan_start.strftime("%d %B %Y %H:%M:%S UTC")}\nThrough {scan_end.strftime("%d %B %Y %H:%M:%S UTC")}', loc='right')
    
    file_name = sat_id + "_GLM_" + scan_end.strftime('%Y%m%d_%H%M%S%Z')
    
    if not os.path.exists(save_dir):
        os.makedirs(save_dir)

    dest_path = os.path.join(save_dir, file_name + ".png")

    plt.savefig(dest_path, bbox_inches="tight", dpi=FILE_DPI)
    plt.close()

def plot_cumulative_glm(file_path: str, 
                        save_dir: str, 
                        points: list[tuple[float, float, str], ...], 
                        bbox: list[float, float, float, float], 
                        hist_flashes_df: pd.DataFrame,
                        cum_minutes: int = DEFAULT_HISTORY_MINUTES, 
                        **kwargs) -> pd.DataFrame:
    
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

    fig = plt.figure(figsize=(22,16))
    ax = plt.axes(projection = crs.PlateCarree())
    
    bbox_wesn = [bbox[3], bbox[2], bbox[1], bbox[0]]
    ax.set_extent(bbox_wesn, crs=crs.PlateCarree())

    draw_logo(ax)
    
    states = NaturalEarthFeature(category="cultural", scale="50m",
                                          facecolor="none",
                                          name="admin_1_states_provinces")
    ax.add_feature(states, linewidth=1.0, edgecolor="black")
    ax.coastlines('50m', linewidth=1.5)
    ax.add_feature(cartopy.feature.LAKES.with_scale('10m'), linestyle='-', linewidth=0.5, alpha=1,edgecolor='blue',facecolor='none')
    ax.add_feature(cfeature.BORDERS, linewidth=1.5)
    plot_towns(ax, (bbox_wesn[0], bbox_wesn[1]), (bbox_wesn[2], bbox_wesn[3]), scale_rank=7)
    
    flashes = ax.plot(flash_lons, 
                      flash_lats,
                      linestyle='None',
                      transform=crs.PlateCarree(), 
                      zorder=15)

    plt.setp(flashes, **NEW_FLASH_STYLE)

    #create flash tuples

    this_frame_data = {'flash_lat': flash_lats, 'flash_lon': flash_lons, 'flash_time': scan_start}
    this_frame_df = pd.DataFrame(data=this_frame_data)

    if not start_of_period:

        hist_flashes = ax.plot(hist_flash_filtered['flash_lon'], 
                               hist_flash_filtered['flash_lat'],
                               linestyle='none',
                               transform=crs.PlateCarree(), 
                               zorder=10)

        plt.setp(hist_flashes, **OLD_FLASH_STYLE)
        hist_flashes_df = pd.concat([hist_flashes_df, this_frame_df], ignore_index=True)
    else:
        hist_flashes_df = this_frame_df

    if points:
        input_pts = []
        input_pt_labels = []
        for point in points:
            x_axis = point[0]
            y_axis = point[1]
            label = point[2]
            
            pt = ax.plot([y_axis],[x_axis])
            input_pts += [pt]
           
            if POINT_LABEL_VISIBLE:
                pt_label = ax.annotate(label, 
                                       (y_axis, x_axis), 
                                       horizontalalignment='center', 
                                       verticalalignment='top', 
                                       transform=crs.PlateCarree(), 
                                       annotation_clip=True, 
                                       zorder=30)
                input_pt_labels = [pt_label]

        plt.setp(input_pts, **POINT_STYLE)
        plt.setp(input_pt_labels, **POINT_LABEL_STYLE)
        adjust_text(input_pt_labels)


    plt.title(f'{orbital_slot.replace("-Test", "")} ({sat_id.replace("G", "GOES-")})\nCumulative Detected Lightning "Flashes"', loc='left', fontweight='bold', fontsize=15)
    plt.title(f'Starting At {period_st.strftime("%d %B %Y %H:%M:%S UTC")}\nThrough {scan_end.strftime("%d %B %Y %H:%M:%S UTC")}', loc='right')
    
    file_name = sat_id + "_GLM_" + period_st.strftime('%Y%m%d_%H%M%S%Z') + "_to" + scan_end.strftime('%Y%m%d_%H%M%S%Z')
    
    if not os.path.exists(save_dir):
        os.makedirs(save_dir)

    dest_path = os.path.join(save_dir, file_name + ".png")

    plt.savefig(dest_path, bbox_inches="tight", dpi=FILE_DPI)
    plt.close()

    return hist_flashes_df

if __name__ == "__main__":

    parser = argparse.ArgumentParser(description=f'Plot GOES 16/17/18 GLM Data, either for single time periods or cumulative over all input data.')
    parser.add_argument('-s', '--single', 
                        help='plot flashes for single-time periods', 
                        action='store_true',
                        default=False)
    parser.add_argument('-c', '--cumulative', 
                        help='plot cumulative flashes over all time periods', 
                        action='store_true',
                        default=False)
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
    parser.add_argument('--save-kmz', 
                        help='save all detected flashes in a kmz', 
                        action='store_true',
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
    else:
        raise ValueError("An option for plotting single-time or cumulative flashes must be specified.")

    points = []
    if args.points_file:
        try:
            with open(args.points_file, newline='') as pointcsv:
                reader = csv.reader(pointcsv, delimiter=',')
                for row in reader:
                    point = (float(row[0]), float(row[1]), row[2])
                    points += [point]
        except Exception as err:
            raise err
    else:
        pass

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
                                                start_of_period=start_of_period)
                if start_of_period:
                    start_of_period = False
            else:
                plot_single_time_glm(input_file_path, 
                                     save_dir, 
                                     points, 
                                     bbox)
            progress.update()

    if args.save_kmz:
        hist_data.to_csv(os.path.join(save_dir, "GLM_all_flashes.csv"), index=None)

    print("Done!")






