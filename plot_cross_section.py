import numpy as np
import pandas as pd
import glob
import io
import os
from copy import copy
from matplotlib import pyplot
from matplotlib.cm import get_cmap
from matplotlib.colors import from_levels_and_colors
import matplotlib.colors as mcolors
from cartopy import crs
from cartopy.feature import NaturalEarthFeature, COLORS
from netCDF4 import Dataset
from wrf import (getvar, to_np, get_cartopy, latlon_coords, vertcross, smooth2d,
                 cartopy_xlim, cartopy_ylim, interpline, CoordPair, ALL_TIMES, xy, ll_to_xy, interp2dxy)
import subprocess

def plot_wrf_cross_section()

if __name__ == "__main__":


fp = False

gif = []

gif = False

DEST_DIR = '/Users/rpurciel/Documents/Missoula MT Collapse/WRF/New/VV/'

files = sorted(glob.glob("/Users/rpurciel/Documents/Missoula MT Collapse/WRF/New/Data/*"))

#flight_path = pd.read_csv("/Users/mikhailk/Desktop/Learjet_MX_FP_X_Sec.csv", sep=',', header=0)

for i in files:
    # Open the NetCDF file
    filename = i
    ncfile = Dataset(filename)
    times = ncfile.variables['Times'][:]
    lat = ncfile.variables['XLAT'][:]
    lon = ncfile.variables['XLONG'][:]
    print(np.min(lat),np.max(lat))
    print(np.min(lon),np.max(lon))
    
    # Get the WRF variables
    ht = getvar(ncfile, "z", units="ft")
    ter = getvar(ncfile, "ter", units="ft")
    Theta = getvar(ncfile, "theta", timeidx=ALL_TIMES)
    wa = getvar(ncfile, "wa", units="ft/s", timeidx=ALL_TIMES)
    ua = getvar(ncfile, "ua", units="kt", timeidx=ALL_TIMES)
    va = getvar(ncfile, "va", units="kt", timeidx=ALL_TIMES)
    
    # Define the cross section start and end points
    
    # startlat =  20.093759 #Learjet Mexico
    # startlon = -100.245835
    
    # endlat = 20.587504
    # endlon = -100.829876
    
    
    # startlat = 35.719396 #NC Fires 2022 X_Sec2
    # startlon = -118.074383
    
    # endlat = 35.718029
    # endlon = -117.522345

    endlat = 46.533340 #2_3
    endlon = -115.42
    
    startlat = 46.533277 #2_3
    startlon = -116.42
    
    start_point = CoordPair(lat=startlat, lon=startlon) 
    end_point = CoordPair(lat=endlat, lon=endlon)
    
    # Get grid points for the Wind Cross Section    
    
    start = ll_to_xy(ncfile, startlat, startlon)
    end = ll_to_xy(ncfile, endlat, endlon)
        
    print(np.shape(Theta))
    print(np.shape(wa))
    
    # Compute the vertical cross-section interpolation.  Also, include the
    # lat/lon points along the cross-section in the metadata by setting latlon
    # to True.
    
    for it in range(0,len(times[:,0])):
        timestr=''
        for istrin in times[it,:]:
            timestr=timestr+istrin.decode('ascii')
        print(it,timestr)
    
        theta = Theta[:,:]
        wa = wa[:,:]
        ua = ua[:,:]
        va = va[:,:]
        
        #Create Vertical Velocity Cross Section 
        wa_cross = vertcross(wa, ht, wrfin=ncfile, start_point=start_point,
                            end_point=end_point, latlon=True, meta=True)
        
        # To remove the slight gap between the vertical velocity contours and terrain due to the
        # contouring of gridded data, a new vertical grid spacing, and model grid
        # staggering, fill in the lower grid cells with the first non-missing value
        # for each column.
        
        # Make a copy of the z cross data. Let's use regular numpy arrays for this.
        wa_cross_filled = np.ma.copy(to_np(wa_cross))
        
        # For each cross section column, find the first index with non-missing
        # values and copy these to the missing elements below.
        for i in range(wa_cross_filled.shape[-1]):
            column_vals = wa_cross_filled[:,i]
            # Let's find the lowest index that isn't filled. The nonzero function
            # finds all unmasked values greater than 0. Since 0 is a valid value
            # for dBZ, let's change that threshold to be -100 instead.
            first_idx = int(np.transpose((column_vals > -100).nonzero())[0])
            wa_cross_filled[0:first_idx, i] = wa_cross_filled[first_idx, i]
        
        ter_line = interpline(ter, wrfin=ncfile, start_point=start_point,
                          end_point=end_point)
        
        # Create Cross Section Lines
        xy_line = xy(theta, start_point=start, end_point=end)
        xy_linez = xy(ht, start_point=start, end_point=end)
        
        # Create Theta Contour via inter2dxy
        z_cross = interp2dxy(theta, xy_line)
        
        # Vertical Coordinate for Theta
        zz_cross = interp2dxy(ht, xy_linez)
        
        # z_cross = vertcross(theta, ht, wrfin=ncfile,
        #                     start_point=start_point,
        #                     end_point=end_point,
        #                     latlon=True, meta=True)
        
        # # Add back the attributes that xarray dropped from the operations above
        # z_cross.attrs.update(z_cross.attrs)
        # z_cross.attrs["description"] = "theta e cross section"
        # z_cross.attrs["units"] = "K"
        
        # # To remove the slight gap between the z cross contours and terrain due to the
        # # contouring of gridded data, a new vertical grid spacing, and model grid
        # # staggering, fill in the lower grid cells with the first non-missing value
        # # for each column.
        
        # # Make a copy of the z cross data. Let's use regular numpy arrays for this.
        # z_cross_filled = np.ma.copy(to_np(z_cross))
        
        # # For each cross section column, find the first index with non-missing
        # # values and copy these to the missing elements below.
        # for i in range(z_cross_filled.shape[-1]):
        #     column_vals = z_cross_filled[:,i]
        #     # Let's find the lowest index that isn't filled. The nonzero function
        #     # finds all unmasked values greater than 0. Since 0 is a valid value
        #     # for dBZ, let's change that threshold to be -200 dBZ instead.
        #     first_idx = int(np.transpose((column_vals > 0).nonzero())[0])
        #     z_cross_filled[0:first_idx, i] = z_cross_filled[first_idx, i]
        
    
        # Create U Winds Cross Section 
        ua_cross = vertcross(ua, ht, wrfin=ncfile, start_point=start_point,
                       end_point=end_point, latlon=True, meta=True)
    
        # To remove the slight gap between the vertical velocity contours and terrain due to the
        # contouring of gridded data, a new vertical grid spacing, and model grid
        # staggering, fill in the lower grid cells with the first non-missing value
        # for each column.
        
        # Make a copy of the z cross data. Let's use regular numpy arrays for this.
        ua_cross_filled = np.ma.copy(to_np(ua_cross))
        
        # For each cross section column, find the first index with non-missing
        # values and copy these to the missing elements below.
        for i in range(ua_cross_filled.shape[-1]):
            column_vals = ua_cross_filled[:,i]
            # Let's find the lowest index that isn't filled. The nonzero function
            # finds all unmasked values greater than 0. Since 0 is a valid value
            # for dBZ, let's change that threshold to be -100 instead.
            first_idx = int(np.transpose((column_vals > -100).nonzero())[0])
            ua_cross_filled[0:first_idx, i] = ua_cross_filled[first_idx, i]
        
        ter_line = interpline(ter, wrfin=ncfile, start_point=start_point,
                          end_point=end_point)
        
        # Create V Winds Cross Section 
        va_cross = vertcross(va, ht, wrfin=ncfile, start_point=start_point,
                           end_point=end_point, latlon=True, meta=True)
        
        # To remove the slight gap between the vertical velocity contours and terrain due to the
        # contouring of gridded data, a new vertical grid spacing, and model grid
        # staggering, fill in the lower grid cells with the first non-missing value
        # for each column.
        
        # Make a copy of the z cross data. Let's use regular numpy arrays for this.
        va_cross_filled = np.ma.copy(to_np(va_cross))
        
        # For each cross section column, find the first index with non-missing
        # values and copy these to the missing elements below.
        for i in range(va_cross_filled.shape[-1]):
            column_vals = va_cross_filled[:,i]
            # Let's find the lowest index that isn't filled. The nonzero function
            # finds all unmasked values greater than 0. Since 0 is a valid value
            # for dBZ, let's change that threshold to be -100 instead.
            first_idx = int(np.transpose((column_vals > -100).nonzero())[0])
            va_cross_filled[0:first_idx, i] = va_cross_filled[first_idx, i]
        
        ter_line = interpline(ter, wrfin=ncfile, start_point=start_point,
                          end_point=end_point)
    
        # Create the figure
        fig = pyplot.figure(figsize=(12,6))
        ax_cross = pyplot.axes()
        ax_cross.set_ylim([3000,11000])
        
        # Make the cross section plot for theta
        thetae_levels = np.arange(270.,760.,2.) #290 - 330 standard & 2 for the interval
        
        shape = z_cross.shape #Getting the shape of the Theta cross section
        xs = np.full((shape),np.arange(0,shape[1],1)) #Creates xs to match ys which matches z_cross' dimensions 
        
        # theta_contours = ax_cross.contour(xs,
        #                                   zz_cross, 
        #                                   to_np(z_cross),
        #                                   levels=thetae_levels,
        #                                   colors="black", alpha = 1.0, linewidths = .4, linestyles = 'dashed')
        # ax_cross.clabel(theta_contours,thetae_levels[::2], inline=1, fontsize=10, fmt="%i")
        
        # wa_cross_filled = wa_cross_filled * 100 * 1.96185 #Convert to ft/min
        print(wa_cross_filled.min(), wa_cross_filled.max())
        
        # Make the contour plot for vertical velocity
        levels = np.arange(-10, 11, 1)
        #levels = np.arange(-30,32,2)
        #levels[5] = .001
        vmin, vmax = levels.min(), levels.max()
        #vmin, vmax = -30, 30
        #norm = mcolors.Normalize(clip=False)
        palette = copy(pyplot.get_cmap("Spectral_r")) #Creating White for Zero
        palette.set_bad('white') #Creating White for Zero
        xs = np.arange(0, wa_cross.shape[-1], 1)
        ys = to_np(wa_cross.coords["vertical"])
        #np.putmask(wa_cross_filled, (wa_cross_filled > -5)&(wa_cross_filled < 5), 'nan')
        wa_contours = ax_cross.contourf(xs, ys, to_np(wa_cross_filled), levels = levels, 
                                        cmap='bwr', vmin = vmin, vmax = vmax, extend = 'both')
        
        # Make the contour line plot for vertical velocity
        levels = np.arange(-10, 11, 1) #np.arange(-500,600,100)
        #levels = np.arange(-30,35,5)
        #levels = [-1.5, -1, -.75, -.5, -.25, 0, .25, .5, .75, 1, 1.5] #Zephyr
        #levels = [-9, -8, -7, -6, -5 , -4, -3, -2, -1, 0, 1, 2, 3, 4] #Thomas Fire
        vmin, vmax = -500, 500
        #vmin, vmax = -30, 30
        xs = np.arange(0, wa_cross.shape[-1], 1)
        
        # ys = to_np(wa_cross.coords["vertical"])
        # np.putmask(wa_cross_filled, (wa_cross_filled > -5)&(wa_cross_filled < 5), 0)
        # wa_cross_filled = smooth2d(wa_cross_filled, 4)
        # wa_contour = ax_cross.contour(xs, ys, to_np(wa_cross_filled), levels = levels, colors = "saddlebrown")
        # ax_cross.clabel(wa_contour,levels[::1], inline=1, fontsize=10, fmt="%1.2i") #Labels every contour because each one is 10cm/s different
        
        # Make the quiver plot for winds
        #ax_cross.quiver(xs[::2], ys[::2], ua_cross_filled[::2,::2], va_cross_filled[::2,::2], pivot = 'middle', scale = 500) #Every other arrow with scaling 
        
        n = 2
        # Add the pressure level wind barbs, only plotting every 6th data point.
        # pyplot.barbs(xs[::n], ys[::n], to_np(ua_cross[::n, ::n]), to_np(va_cross[::n, ::n]), length=6)
        
        # Add the color bar
        pyplot.colorbar(wa_contours, ax=ax_cross)
        
        # Fill in the mountain area
        ht_fill = ax_cross.fill_between(xs, 0, to_np(ter_line),
                                        facecolor="saddlebrown")
        
        # Set the x-ticks to use latitude and longitude labels
        coord_pairs = to_np(wa_cross.coords["xy_loc"])
        x_ticks = np.arange(coord_pairs.shape[0])
        x_labels = [pair.latlon_str(fmt="{:.2f}, {:.2f}") for pair in to_np(coord_pairs)]
        
        # Set the y-ticks to be height.
        #vert_vals = to_np(wa_cross.coords["vertical"])
        #v_ticks = np.arange(vert_vals.shape[0])
        #ax_cross.set_yticks([0,2000,3000,4000,6000,8000,9000,10000,11000,12000,13000,14000])
        #ax_cross.set_yticklabels(vert_vals[::1], fontsize=8)
        
        # Set the desired number of x ticks below
        num_ticks = 7
        # thin = int((len(x_ticks) / num_ticks) + .5)
        thin = 12
        ax_cross.set_xticks(x_ticks[::thin])
        ax_cross.set_xticklabels(x_labels[::thin], rotation=45, fontsize=8)
        
        value = .3
        color = 'black'

        CRASH_lon, CRASH_lat = -114.285, 46.538100
        # Convert to xy grid for plotting on axes
        CRASH = ll_to_xy(ncfile, CRASH_lat, CRASH_lon)
        
        # Plotting the Point
        ax_cross.plot(CRASH[0],[3325], #Only using the first dimension of the CRASH array and second number is the elevation in FT
            color='red', marker='o')

        ax_cross.annotate('Building', (CRASH[0] - 13, 3175), horizontalalignment='center', color='red', fontsize=12, fontweight='bold', zorder=30)
        
        # # Lists for the x & y coordinates for points to be plotted on the cross section. X is along the x-axis and y is in feet. 
        # x_axis_points = [75,91]
        # y_axis_points = [[ter_line[round(x_axis_points[0])]],[ter_line[round(x_axis_points[1])]]]
        # labels = ['Stone', 'Conley']
        
        # for (i,j,k) in zip(x_axis_points, y_axis_points, labels):
        #     x_axis = i
        #     y_axis = j
        #     label = k
            
        #     # Plotting the Point
        #     ax_cross.plot([i],[j], #Only using the first dimension of the CRASH array and second number is the elevation in FT
        #           color='black', marker='o')
            
        #     # Adding Text Label for the Points
        #     ax_cross.annotate(k, xy = (i,j), xytext = (i+value,j), 
        #         horizontalalignment='left', color = color, fontsize = 14)
        
        if fp == True:
            
            lat_total = []
            lon_total = []
            
            for element in x_labels:
                
                point = element.split(",")
                lat_total.append(point[0])
                lon_total.append(point[1])
                
            lat_total_float = list(map(float, lat_total))
            lon_total_float = list(map(float, lon_total))
            
            # Lists for the x & y coordinates for points to be plotted on the cross section. X is along the x-axis and y is in feet. 
            x_axis_points = np.array(flight_path.Latitude)
            y_axis_points = np.array(flight_path.Longitude)
            labels = np.array(flight_path.Time)
            altitudes = np.array(flight_path.Altitude)
            
            for (i,j,k,l) in zip(x_axis_points, y_axis_points, labels, altitudes):
                x_axis = i
                y_axis = j
                label = k
                altitude = l
            
                lat_total_bool = np.isclose(lat_total_float, x_axis, rtol=1e-04, atol=1e-03) #Finds the array in x_y_arr that matches closest to CRASH array and returns True/False for each lat long
                lon_total_bool = np.isclose(lon_total_float, y_axis, rtol=1e-04, atol=1e-04) #Finds the array in x_y_arr that matches closest to CRASH array and returns True/False for each lat long
                
                lat_total_arr = np.where(lat_total_bool) #Finds the True Array
                lon_total_arr = np.where(lon_total_bool) #Finds the True Array
                
                lat_total_idx = lat_total_arr[0][0] #Getting the first index of the True array
                lon_total_idx = lon_total_arr[0][0] #Getting the first index of the True array
                
                if lat_total_idx != lon_total_idx:
                    idx = ((lat_total_idx + lon_total_idx) / 2)
                    
                else:
                    idx = lat_total_idx
                
                #Plotting the Point
                ax_cross.plot([idx],[l], #Only using the first dimension of the CRASH array and second number is the elevation in FT
                      color='black', marker='+', markersize = 5)
                
                #Adding Text Label for the Points
                ax_cross.annotate(k, xy = (idx,l), xytext = (idx+value,l), rotation = 45, 
                      horizontalalignment='left', color = color, fontsize = 8)
        
        ### Leer Jet Mexico X_Sec2
    
        #Plotting the Point
        # point = 27 #Camp was 74
        # ax_cross.plot([point],ter_line[round(point)], #Only using the first dimension of the CRASH array and second number is the elevation in FT
        #       color='black', marker='o')
        
        # #Adding Text Label for the Points
        # ax_cross.annotate("Crash", xy = (point,ter_line[round(point)]), xytext = (point+value,ter_line[round(point)]), 
        #       horizontalalignment='left', color = color, fontsize = 14)
        
        # #Plotting the Point
        # point = 91
        # ax_cross.plot([point],ter_line[round(point)], #Only using the first dimension of the CRASH array and second number is the elevation in FT
        #       color='black', marker='o')
        
        # #Adding Text Label for the Points
        # ax_cross.annotate("Conley", xy = (point,ter_line[round(point)]), xytext = (point+value,ter_line[round(point)]), 
        #       horizontalalignment='left', color = color, fontsize = 14)
        
        ### Leer Jet Mexico X_Sec2
        
        # ax_cross.vlines(12, 0, 20000, color = color, linestyles = 'dashed') #Fort Nelson for Domain 1 8.5
        # ax_cross.vlines(44, 0, 20000, color = color, linestyles = 'dashed') #Fort Nelson for Domain 2 43
        # ax_cross.vlines(155, 0, 1000000, color = color, linestyles = 'dashed') #Fort Nelson for Domain 2
        
        #ax_cross.grid()
        
        # Set the x-axis and  y-axis labels
        ax_cross.set_xlabel("Latitude, Longitude", fontsize=12)
        ax_cross.set_ylabel("Height (ft)", fontsize=12)
        
        #ax_cross.set_title("Vertical Cross Section of Vertical Velocity (cm/s); time="+timestr, fontweight = 'semibold', pad = 10)
        ax_cross.set_title("Vertical Cross Section of Vertical Velocity (ft/min) - "+timestr.replace("_", " ")+ " UTC", fontweight = 'semibold', pad = 10)
        
        buf = io.BytesIO()
        pyplot.savefig(os.path.join(DEST_DIR, "CrossSection.verticalvelocity."+timestr+".jpg") ,bbox_inches="tight",dpi=200)
        pyplot.close()
        buf.seek(0)
        buf.close()
        print("CrossSection.verticalvelocity."+timestr.split(':')[0]+".jpg") 
print("..I finished")

if gif == True:
    print("Making the GIF...")
    cmd = '/usr/local/bin/convert -delay 50 -resize 75% *.jpg /users/mikhailk/Desktop/VVCross.gif'
    subprocess.call(cmd, shell=True)
    print("...I Made the VVCross GIF")