import os
import sys
from typing import Any
from collections.abc import Callable
from uuid import uuid4
from glob import glob 
import requests

import numpy as np
import pandas as pd
import xarray as xr
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
import matplotlib.patheffects as patheffects
from matplotlib.offsetbox import OffsetImage, AnnotationBbox
import cartopy.mpl.geoaxes
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import cartopy.io.shapereader as shpreader
from cartopy.feature import NaturalEarthFeature
from adjustText import adjust_text 

GEONAMES_USERNAME = 'rpurciel'

def destag_variable(stag_var: xr.DataArray, 
                    *,
                    stag_dim_name: str = None) -> xr.DataArray:
    """
    Destagger a WRF variable. Can either automatically
    find the staggered dimensions or a staggered dimension
    can be specified via 'stag_dm_name'.

    Inputs
        - stag_var, xarray.DataArray: 
                An Xarray DataArray for a variable staggered
                on the WRF mass grid.
        - stag_dim_name, str [optional, def. None]:
                If the staggered dimension's name is known,
                it can be manually selected. Otherwise, 
                dimension is selected based on patterns.

    Returns
        An Xarray DataArray, destaggered if it was staggered
        on the WRF mass grid. If no dimensions are staggered
        or if wrfpython does not recognize the data, the 
        input DataArray is returned.
    """

    from wrf import destagger

    dims = stag_var.dims
    stag_dim_idx = None

    if stag_dim_name:
        stag_dim_idx = dims.index(stag_dim_name)
    else:
        for dim in dims:
            if "stag" in dim:
                stag_dim_idx = dims.index(dim)

    if not stag_dim_idx:
        return stag_var
    else:        
        try:
            destag_var = destagger(stag_var, 
                                   stag_dim_idx,
                                   meta=True)
        except:
            return stag_var

        destag_var = destag_var.assign_coords(stag_var.coords)

        return destag_var

def get_nearest_feature(lat: float,
                        lon: float,
                        feature_code: str = 'PPL',
                        show_distance: bool = True):

    keys = {'lat': lat, 'lng': lon, 'username': GEONAMES_USERNAME, 'featureCode': feature_code}

    r = requests.get('http://api.geonames.org/findNearbyJSON', params=keys)
    try:
        params = r.json()['geonames'][0]
    except:
        return "No distinguising features nearby"

    state_code = params['adminCode1']
    if show_distance:
        loc = f"{round(float(params['distance']), 2)} km from {params['name']}{f', {state_code}, ' if state_code is not None else ', '}{params['countryName']}"
    else:
        loc = f"Near {params['name']}{f', {state_code}, ' if state_code is not None else ', '}{params['countryName']}"
    return loc

def make_highly_visible(plt: mpl.figure.Figure,
                        obj: mpl.artist.Artist,
                        outline_width: int = 5,
                        outline_color: str = 'black',
                        object_color: str = 'greenyellow') -> None:

    """
    This function will take an internal matplotlib artist
    object, and will apply a 'highly visible' color and 
    ouline to it. Colors and outline width is customizable.

    Inputs
        - fig, matplotlib.figure.Figure: 
                A Figure object containing the plot.
        - obj, matplotlib.artist.Artist: 
                The object that the outline will be applied to,
                which is a subclass of mpl.artist.Artist. 
        - outline_width, int (optional: def 5):
                The width of the outline around the object.
        - outline_color, str (optional: def 'black')
                The color of the outline around the object.
        - object_color, str (optional: def 'greenyellow')
                The color to be applied to the object.

    Returns
        None
    """

    plt.setp(obj, color=object_color)
    obj.set_path_effects([patheffects.withStroke(linewidth=outline_width, 
                                                 foreground=outline_color)])

def plot_points(plt: mpl.figure.Figure,
                ax: mpl.axes.Axes | cartopy.mpl.geoaxes.GeoAxes, 
                points: list[tuple[float, float, str, str, float, float], ...],
                transform: cartopy.crs = ccrs.PlateCarree(), 
                default_marker: str = '.',
                default_label: str = None,
                default_x_offset: float = 0,
                default_y_offset: float = 0,
                zorder: int = 30,
                **kwargs) -> tuple[list, list]:

    input_pts = []
    input_pt_labels = []
    for point in points:
        x_axis = point[0]
        y_axis = point[1]
        marker = point[2] if point[2] is not None else default_marker
        label = point[3] if point[3] is not None else default_label
        x_offset = point[4] if point[4] is not None else default_x_offset
        y_offset = point[5] if point[5] is not None else default_y_offset
        
        pt = ax.plot([y_axis],[x_axis], 
                     marker=marker, 
                     linestyle='none',
                     markerfacecolor=kwargs.get("point_color", "black"),
                     markeredgecolor=kwargs.get("point_edge_color", "black"),
                     markersize=kwargs.get("point_size", 8),
                     zorder=zorder)

        input_pts += [pt]
       
        if (kwargs.get("point_draw_labels", True)) and (label is not None):
            if kwargs.get("point_draw_arrows", True):
                pt_label = ax.annotate(label, 
                                       xy=(y_axis, x_axis),
                                       xytext=(y_axis+float(x_offset), x_axis+float(y_offset)), 
                                       horizontalalignment=kwargs.get("point_label_horizontal_alignment", "right"), 
                                       verticalalignment=kwargs.get("point_label_vertical_alignment", "top"),
                                       fontsize=kwargs.get("point_label_fontsize", 12),
                                       color=kwargs.get("point_label_color", "black"),
                                       arrowprops=dict(arrowstyle='->', 
                                                       color=kwargs.get("point_label_color", "black")), 
                                       transform=transform,
                                       annotation_clip=True, 
                                       zorder=zorder)
                if kwargs.get("point_label_high_vis", False):
                    make_highly_visible(plt, pt_label)
            else:
                pt_label = ax.annotate(label, 
                                       xy=(y_axis, x_axis),
                                       xytext=(y_axis+float(y_offset), x_axis+float(x_offset)), 
                                       horizontalalignment=kwargs.get("point_label_horizontal_alignment", "right"), 
                                       verticalalignment=kwargs.get("point_label_vertical_alignment", "top"),
                                       fontsize=kwargs.get("point_label_fontsize", 12),
                                       color=kwargs.get("point_label_color", "black"),
                                       transform=transform,
                                       annotation_clip=True, 
                                       zorder=zorder)
                if kwargs.get("point_label_high_vis", False):
                    make_highly_visible(plt, pt_label)
            input_pt_labels = [pt_label]

    if kwargs.get("point_label_autoadjust_pos", False): #can create wierd behavior when specifying point label positions manually
        adjust_text(input_pt_labels)

    return (input_pts, input_pt_labels)

def plot_towns(ax: mpl.axes.Axes | cartopy.mpl.geoaxes.GeoAxes, 
               lon_bounds: tuple[float, float],
               lat_bounds: tuple[float, float], 
               resolution: str = '10m',
               scale_rank: int = 5, 
               transform: cartopy.crs = ccrs.PlateCarree(), 
               zorder: int = 1) -> None:
    """
    This function will download the 'populated_places' shapefile from
    NaturalEarth, trim the shapefile based on the limits of the provided
    lat & long coords, and then plot the locations and names of the towns
    on a given GeoAxes.

    Inputs
        - ax, cartopy.mpl.geoaxes.GeoAxes: 
                A GeoAxes object to plot onto.
        - lat_bounds, tuple: 
                A tuple of size 2 containing the latitude bounds of the plot area.
                Standard ordering is from lower to upper bound.
        - lon_bounds, tuple: 
                A tuple of size 2 containing the longitude bounds of the plot area.
                Standard ordering is from lower to upper bound.
        - resolution, str (optional, def '10m'):
                A string specifying what resolution of town data to download. Default
                is the finest resolution.
        - scale_rank, int (optional, def 5):
                Int representing the upper bound of the 'scale rank' of the towns to
                plot, from most-to-least culturally relevant. e.g. 1-2 will plot only
                the largest towns, >5 will plot tiny towns.
        - transform, cartopy.crs (optional, def ccrs.PlateCarree()):
                A transformation to use for the map. Should match transformation for
                GeoAxes object.
        - zorder, int (optional, def 3):
                Int specifying the zorder value to be fed into matplotlib.
    
    Returns
        None
    """

    shp = shpreader.Reader(shpreader.natural_earth(resolution=resolution, category='cultural', name='populated_places'))

    #get town names
    names = []
    x_pts = []
    y_pts = []

    for town in shp.records():
        if int(town.attributes['SCALERANK']) <= scale_rank:
            x_pts.append(town.attributes['LONGITUDE'])
            y_pts.append(town.attributes['LATITUDE'])
            name = town.attributes['NAME_EN']
            names.append(f'{name}')

    #create data frame and index by the region of the plot
    all_towns = pd.DataFrame({'names': names, 'x': x_pts, 'y': y_pts})
    region_towns = all_towns[(all_towns.y<max(lat_bounds)) & (all_towns.y>min(lat_bounds))
                           & (all_towns.x>min(lon_bounds)) & (all_towns.x<max(lon_bounds))]

    #plot the locations and labels of the towns in the region
    ax.scatter(region_towns.x.values, region_towns.y.values, 
               s=10,
               c ='black', 
               marker= '.', 
               transform=transform, 
               zorder=zorder)

    town_names = []
    for row in region_towns.itertuples():
        name = ax.text(float(row.x), 
                       float(row.y) * 0.9995, 
                       row.names,
                       fontsize=9, 
                       transform=transform,
                       style='italic',
                       horizontalalignment='left',
                       verticalalignment='top',
                       clip_box=ax.bbox)
        town_names.append(name)

    #use adjustText library to autoadjust town names to prevent ugly overlapping
    adjust_text(town_names)

def draw_logo(ax: mpl.axes.Axes | cartopy.mpl.geoaxes.GeoAxes, 
              logo_rel_path: str = 'logo.png',
              corner: str = 'lowerleft',
              scale: float = 1,
              zorder: int = 30) -> None:

    """
    This function will use a locally stored logo image, and 
    then plot the image on a corner of the plot.

    Inputs
        - ax, matplotlib.axes.Axes or cartopy.mpl.geoaxes.GeoAxes: 
                An Axes or GeoAxes object to plot onto.
        - logo_rel_path, str (optional, def 'logo.png'):
                A string specifying what file to use as a logo. Files are stored
                in the 'logos' subdirectory.
        - corner, str (optional, def 'lowerleft')
                A string specifying what corner (of the internal bounding box) to
                place the logo at. Default is the lower left corner.
        - scale, int (optional, def 1):
                Optional scaling factor of the logo's size.
        - zorder, int (optional, def 30):
                Int specifying the zorder value to be fed into matplotlib.

    Returns
        None
    """

    mpl.use('agg')


    try:
        with open(os.path.join(sys.path[0], f'logos/{logo_rel_path}'), 'rb') as file:
            logo_img = plt.imread(file)
    except Exception as OpenExcept:
        raise OpenExcept

    logobox = OffsetImage(logo_img, zoom=0.05 * scale)
    logobox.image.axes = ax

    if 'lower' in corner:
        logo_vpos = 0
    elif 'upper' in corner:
        logo_vpos = 1
    else:
        logo_vpos = 0

    if 'left' in corner:
        logo_hpos = 0
    elif 'right' in corner:
        logo_hpos = 1
    else:
        logo_hpos = 0

    logo = AnnotationBbox(logobox, 
                         (logo_hpos, logo_vpos),
                         xybox=(10., 10.),
                         xycoords='axes fraction',
                         boxcoords="offset points",
                         box_alignment=(0, 0),
                         pad=0.0,
                         frameon=False,
                         zorder=zorder)

    ax.add_artist(logo)

def define_hi_res_fig(lon_bounds: tuple[float, float],
                      lat_bounds: tuple[float, float],
                      projection: cartopy.crs = ccrs.PlateCarree(),
                      draw_earth_features: bool = True,
                      draw_gridlines: bool = False) -> tuple[mpl.figure.Figure, mpl.axes.Axes]:

    fig = plt.figure(figsize=(22,16))
    ax = plt.axes(projection = ccrs.PlateCarree())
    ax.set_extent([min(lon_bounds), max(lon_bounds), min(lat_bounds), max(lat_bounds)], 
                  crs=ccrs.PlateCarree())

    if draw_earth_features:
        states = NaturalEarthFeature(category="cultural", scale="50m",
                                              facecolor="none",
                                              name="admin_1_states_provinces")
        ax.add_feature(states, linewidth=1.0, edgecolor="black")
        ax.coastlines('50m', linewidth=1.5)
        ax.add_feature(cartopy.feature.LAKES.with_scale('10m'), linestyle='-', linewidth=0.5, alpha=1,edgecolor='blue',facecolor='none')
        ax.add_feature(cfeature.BORDERS, linewidth=1.5)

    if draw_gridlines:
        ax.set_xlabel('Latitude')
        ax.set_ylabel('Longitude')

        gl = ax.gridlines(crs=ccrs.PlateCarree(), 
                          linewidth=1, 
                          color='black', 
                          alpha=0.5, 
                          linestyle='--', 
                          draw_labels=True)

        gl.xlabels_top = False
        gl.ylabels_left = True
        gl.ylabels_right= False
        gl.xlines = True
        gl.ylines = True

        gl.xlocator = mticker.FixedLocator(np.arange(min(lon_bounds), 
                                                     max(lon_bounds), 
                                                     2))
        gl.ylocator = mticker.FixedLocator(np.arange(min(lat_bounds), 
                                                     max(lat_bounds), 
                                                     1))

    return fig, ax

def define_gearth_compat_fig(lon_bounds: tuple[float, float],
                             lat_bounds: tuple[float, float],
                             projection: cartopy.crs = ccrs.PlateCarree(),
                             draw_earth_features: bool = False,
                             draw_gridlines: bool = False,
                             dpi_pixels: int = 1024) -> tuple[mpl.figure.Figure, mpl.axes.Axes, mpl.figure.Figure, mpl.axes.Axes]:

    fig_x_size = np.ptp([max(lon_bounds), min(lon_bounds)]) * np.cos(np.mean([min(lat_bounds), max(lat_bounds)]) * np.pi/180.)
    fig_y_size = np.ptp([max(lat_bounds), min(lat_bounds)])
    fig_aspect_ratio = fig_x_size / fig_y_size

    if fig_aspect_ratio > 1.0:
        fig_size = (10.0 / fig_aspect_ratio, 10.0)
    else:
        fig_size = (10.0, 10.0 * fig_aspect_ratio)

    fig = plt.figure(figsize=fig_size,
                     frameon=False,
                     dpi=dpi_pixels // 10)

    ax = fig.add_axes([0, 0, 1, 1], projection=projection)
    ax.set_xlim(min(lon_bounds), max(lon_bounds))
    ax.set_ylim(min(lat_bounds), max(lat_bounds))

    if draw_earth_features:
        states = NaturalEarthFeature(category="cultural", scale="50m",
                                              facecolor="none",
                                              name="admin_1_states_provinces")
        ax.add_feature(states, linewidth=1.0, edgecolor="black")
        ax.coastlines('50m', linewidth=1.5)
        ax.add_feature(cartopy.feature.LAKES.with_scale('10m'), linestyle='-', linewidth=0.5, alpha=1,edgecolor='blue',facecolor='none')
        ax.add_feature(cfeature.BORDERS, linewidth=1.5)

    if draw_gridlines:
        ax.set_xlabel('Latitude')
        ax.set_ylabel('Longitude')

        gl = ax.gridlines(crs=projection, 
                          linewidth=1, 
                          color='black', 
                          alpha=0.5, 
                          linestyle='--', 
                          draw_labels=True)

        gl.xlabels_top = False
        gl.ylabels_left = True
        gl.ylabels_right= False
        gl.xlines = True
        gl.ylines = True

        gl.xlocator = mticker.FixedLocator(np.arange(min(lon_bounds), 
                                                     max(lon_bounds), 
                                                     2))
        gl.ylocator = mticker.FixedLocator(np.arange(min(lat_bounds), 
                                                     max(lat_bounds), 
                                                     1))

    ax.set_axis_off()

    ##Defining a figure to be used as a "colorbar" overlay,
    ##which can be pinned to the edge of the screen when saved

    cbfig = plt.figure(figsize=(1.0, 4.0), 
                       facecolor=None, 
                       frameon=True)
    cbax = cbfig.add_axes([0.0, 0.05, 0.2, 0.9])
    cbax.set_facecolor('white')

    return fig, ax, cbfig, cbax   

def save_figs_to_kml(save_path: str,
                     lon_bounds: tuple[float, float],
                     lat_bounds: tuple[float, float],
                     layer_figs: list[mpl.figure.Figure, ...] | list[str, ...],
                     layer_names: list[str, ...],
                     layer_descs: list[str, ...],
                     colorbar_fig: mpl.figure.Figure | str = None,
                     **kwargs) -> None:

    from simplekml import (Kml, OverlayXY, ScreenXY, Units, RotationXY,
                       AltitudeMode, Camera)

    paths_to_figs = []
    if type(layer_figs) is mpl.figure.Figure:
        for figure in [layer_figs]:
            fig_path = os.path.join(os.getcwd(), f'TEMP_{uuid4().hex}.png')
            figure.savefig(fig_path, 
                        bbox_inches="tight", 
                        dpi=600, 
                        transparent=False)
            paths_to_figs += [fig_path]
    elif type(layer_figs) is str:
        paths_to_figs = layer_figs

    if type(colorbar_fig) is mpl.figure.Figure:
        fig_path = os.path.join(os.getcwd(), f'TEMP_{uuid4().hex}.png')
        colorbar_fig.savefig(fig_path, 
                             bbox_inches="tight", 
                             dpi=150, 
                             transparent=False,
                             facecolor='White',
                             edgecolor='White')
        path_to_cb = fig_path
    elif type(colorbar_fig) is str:
        path_to_cb = colorbar_fig

    kml = Kml()
    camera = Camera(latitude=np.mean([max(lat_bounds), min(lat_bounds)]),
                    longitude=np.mean([max(lon_bounds), min(lon_bounds)]),
                    altitude=kwargs.pop('altitude', 2e7), 
                    roll=kwargs.pop('roll', 0), 
                    tilt=kwargs.pop('tilt', 0),
                    altitudemode=kwargs.pop('altitudemode', AltitudeMode.relativetoground))

    kml.document.camera = camera
    draworder = 0
    for fig, name, desc in zip(paths_to_figs, layer_names, layer_descs):  # NOTE: Overlays are limited to the same bbox.
        draworder += 1
        ground = kml.newgroundoverlay(name='GroundOverlay')
        ground.draworder = draworder
        ground.visibility = kwargs.pop('visibility', 1)
        ground.name = name
        ground.color = kwargs.pop('color', 'ffffffff')
        ground.atomauthor = kwargs.pop('author', 'WeatherExtreme')
        ground.latlonbox.rotation = kwargs.pop('rotation', 0)
        ground.description = desc
        ground.gxaltitudemode = kwargs.pop('gxaltitudemode', 'clampToSeaFloor')
        ground.opacity = kwargs.pop('opacity', 1)
        ground.icon.href = fig
        ground.latlonbox.east = max(lon_bounds)
        ground.latlonbox.south = min(lat_bounds)
        ground.latlonbox.north = max(lat_bounds)
        ground.latlonbox.west = min(lon_bounds)

    if colorbar_fig:  # Options for colorbar are hard-coded (to avoid a big mess).
        screen = kml.newscreenoverlay(name='Legend')
        screen.icon.href = path_to_cb
        screen.overlayxy = OverlayXY(x=0, y=0,
                                     xunits=Units.fraction,
                                     yunits=Units.fraction)
        screen.screenxy = ScreenXY(x=0.015, y=0.075,
                                   xunits=Units.fraction,
                                   yunits=Units.fraction)
        screen.rotationXY = RotationXY(x=0.5, y=0.5,
                                       xunits=Units.fraction,
                                       yunits=Units.fraction)
        screen.size.x = 0
        screen.size.y = 0
        screen.size.xunits = Units.fraction
        screen.size.yunits = Units.fraction
        screen.visibility = 1

    kml.savekmz(save_path)
    if type(layer_figs) is mpl.figure.Figure:
        for figure in paths_to_figs:
            os.remove(figure)

    if colorbar_fig:
        os.remove(path_to_cb)

def generate_verification_plot(lat: float | list[float, ...], 
                               lon: float | list[float, ...], 
                               save_dir: str,
                               extent_factor: int = None) -> None:
    """
    This function generates a plot of the input lat lon,
    to be used for verification of the input location.
    Plot is saved to [save_dir]/location_verification.png

    Inputs
        - lat, float: 
                Latitude of requested point. A list of latitudes is supported,
                if one is passed plot is centered around the last point.
        - lon, float: 
                Longitude of requested point. A list of latitudes is supported,
                if one is passed plot is centered around the last point.
        - save_dir, str: 
                The desired plot save destination.
        - extent_factor, int (optional: def 5):
                The factor (in degrees) to calculate the extent of the plot.
                e.g. 5 means 5 degrees out from the center on each side.

    Returns
        None
    """

    mpl.use('agg')

    #Generate figure and populate with features (borders, water, towns, etc.)
    fig = plt.figure(figsize=(22,16))
    ax = plt.axes(projection = ccrs.PlateCarree())

    if (type(lat) or type(lon)) is list:

        lon_dist = np.abs(max(lon) - min(lon))
        lat_dist = np.abs(max(lat) - min(lat))

        rad_max = lon_dist if lon_dist >= lat_dist else lat_dist
        fact = (rad_max + rad_max * 0.05) if extent_factor is None else extent_factor

        center_lat = np.mean(lat)
        center_lon = np.mean(lon)

        extent = [center_lon-fact,
                  center_lon+fact,
                  center_lat-fact,
                  center_lat+fact]
    else:
        fact = 0.5 if extent_factor is None else extent_factor

        extent = [float(lon)-fact,
                  float(lon)+fact,
                  float(lat)-fact,
                  float(lat)+fact]

    ax.set_extent(extent, crs=ccrs.PlateCarree())
    states = NaturalEarthFeature(category="cultural", scale="50m",
                                                  facecolor="none",
                                                  name="admin_1_states_provinces")
    ax.add_feature(states, linewidth=1.0, edgecolor="black")
    ax.coastlines('50m', linewidth=1.5)
    ax.add_feature(cartopy.feature.LAKES.with_scale('10m'), linestyle='-', linewidth=0.5, alpha=1,edgecolor='blue',facecolor='none')
    ax.add_feature(cfeature.BORDERS, linewidth=1.5)
    ax.add_feature(cfeature.LAND)
    plot_towns(ax, extent[0:2], extent[2:])

    #Add latitude/longitude gridlines and labels
    ax.set_xlabel('Latitude')
    ax.set_ylabel('Longitude')
    gl = ax.gridlines(crs=ccrs.PlateCarree(), linewidth=1, color='black', alpha=0.5, linestyle='--', draw_labels=True)
    gl.right_labels = True
    gl.left_labels = True
    gl.bottom_labels = True
    gl.top_labels = True
    gl.xlocator = mticker.FixedLocator(np.arange(extent[0], extent[1], 2))
    gl.ylocator = mticker.FixedLocator(np.arange(extent[2], extent[3], 2))

    #Plot requested point(s) and save
    if (type(lat) or type(lon)) is list:
        plt.plot(lon, lat, marker='x', color='red', transform=ccrs.PlateCarree(), zorder=15)
        # plt.annotate(f"Beginning\n{lat[len(lat)-1]}, {lon[len(lon)-1]}", 
        #              (lon[len(lon)-1], lat[len(lat)-1]-0.3), 
        #              color='red', 
        #              fontweight='bold', 
        #              horizontalalignment='center', 
        #              verticalalignment='top',
        #              transform=ccrs.PlateCarree(), 
        #              zorder=15)
        plt.title(f'Beginning Point:\n{lat[0]}, {lon[0]}', loc='left')
        plt.title(f'Ending Point:\n{lat[len(lat)-1]}, {lon[len(lon)-1]}', loc='right')
    else:
        print(lon, lat)
        plt.plot(float(lon), float(lat), marker='x', color='red', transform=ccrs.PlateCarree(), zorder=15)
        plt.annotate(f"Selected Point\n{lat}, {lon}", 
                     (float(lon), float(lat)-0.3), 
                     color='red', 
                     fontweight='bold', 
                     horizontalalignment='center', 
                     verticalalignment='top',
                     transform=ccrs.PlateCarree(), 
                     zorder=15)

    draw_logo(ax)
    plt.savefig(os.path.join(save_dir, "location_verification.png"), bbox_inches="tight", dpi=200)
    
    plt.close()

def remove_files(target_directory: str,
                 file_ext: str = 'idx',
                 verbose: bool = False,
                 verbose_callback: Callable[[str], Any] = print) -> int:

    """
    This function will take a target directory, and remove 
    all files with a given extension in that directory. Files
    removed can be logged/printed using a callback function.

    Inputs
        - target_directory, str:
                The directory to remove files from
        - file_ext, str: (optional: def 'idx'): 
                A string with the extension of files to
                target, defauling to index ('idx') files.
        - verbose, bool: (optional: def False)
                True/False to turn on verbose mode, which
                sends certain event strings to a callback
                function.
        - verbose_callback, Callable[[str], Any] (optional: def print)
                The callback function to use for verbose output, 
                which must take one string as an input. Default
                is printing to stdout.

    Returns
        Int: Number of files successfully removed, or
             -1 in the case of an error.
    """

    if target_directory[-1] != '/':
        target_directory += '/'

    if not os.path.exists(target_directory):
        if verbose:
            verbose_callback('Error: Directory to remove files from does not exist.')
        return -1

    files = glob(f'{target_directory}*.{file_ext}')
    num_files = len(files)
    removed_files = 0

    if not files:
        if verbose:
            verbose_callback(f'Error: No files found in directory, or wrong file extension (.{file_ext}).')
        return -1
    else:
        for file in files:
            os.remove(file)
            if verbose:
                verbose_callback(f"Sucessfully removed file '{file}'")
            removed_files += 1

        if verbose:
            verbose_callback(f"Removed {removed_files} files out of {num_files} selected.")

        return removed_files
