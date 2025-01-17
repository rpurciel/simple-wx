## Library Info

<details>
<summary>Completeness</summary>

| Script | Completeness | Next Goal |
|---|:---:|---|---|
|`airmet_to_json_wizard.py`| 🟩🟩🟩🟩🟩🟩🟩🟩🟩⬜ | Clean code |
|`create_meteogram_data.py`| 🟩🟩🟩🟩🟩⬜⬜⬜⬜⬜ | Model support (ERA5/GFS/WRF/+), other variables, clean code |
|`create_movie.py`| 🟩🟩🟩🟩🟩🟩🟩⬜⬜⬜ | Clean code, options for more video formats |
|`create_raob_sounding.py`| 🟩🟩🟩🟩⬜⬜⬜⬜⬜⬜ | Model support (GFS/WRF/+), clean code |
|`create_sun_moon_pos_table.py`| 🟩🟩🟩🟩🟩🟩🟩🟩🟩⬜ | Elevation support |
|`create_twilight_times_listing.py`| 🟩🟩🟩🟩🟩🟩🟩🟩🟩⬜ | Elevation support, clean code |
|`download_aws.py`| 🟩🟩🟩🟩🟩🟩🟩⬜⬜⬜ | Clean code, model forecast data support |
|`download_eumetsat_from_cart.py`| 🟩🟩🟩🟩🟩🟩🟩🟩⬜⬜ | Clean code |
|`plot_airmet.py`| 🟩🟩🟩🟩⬜⬜⬜⬜⬜⬜ | Clean code, make more general purpose (CSIGs, etc) |
|`plot_cross_section.py`| 🟩⬜⬜⬜⬜⬜⬜⬜⬜⬜ | Build functionality |
|`plot_glm.py`| 🟩🟩🟩🟩🟩🟩🟩⬜⬜⬜ | Clean code, update |
|`plot_goes.py`| 🟩🟩🟩🟩🟩🟩🟩🟩⬜⬜ | Clean code, add in addl. composites |
|`plot_himawari.py`| 🟩🟩🟩🟩🟩🟩🟩⬜⬜⬜ | Clean code, add in addl. composites |
|`plot_meteosat.py`| 🟩🟩🟩🟩🟩🟩🟩⬜⬜⬜ | Clean code, add in addl. composites |
|`plot_plan_view.py`| 🟩🟩🟩🟩🟩⬜⬜⬜⬜⬜ | Model support (GFS/WRF/+), surface-level support, addl. variables |
|`plot_radar.py`| 🟩🟩🟩🟩🟩🟩⬜⬜⬜⬜ | Variable suppport |

</details>

## Downloading Scripts

<details>
<summary>download_aws.py</summary>

A script to download certain products from AWS S3 buckets.

> [!NOTE]
> Implemented products are: Multiband ABI and GLM data from GOES 16, 17, and 18,
> AHI data from Himawari 8 and 9, zero-hour analysis data from the GFS and HRRR models,
> Level 2 NEXRAD data by radar ID, and multi-radar multi-sensor data.

Output → Files downloaded from AWS of the selected product, filtered by time.

Options

| Option | Variable/Flag | Input | Description | 
|---|---|---|---|
| save_directory | _(required)_ | path | Directory to save downloaded files to |
| satellite | -s, --satellite | goes16, goes17, goes18, goes16-glm, goes17-glm, goes18-glm, hw8, hw9 | Specify a satellite product to download |
| model | -m, --model | gfs-anl, hrrr-anl | Specify a model product to download |
| radar | -r, --radar | radar_id (string) | Specify a radar to download data from ("mrms" for MRMS data) |
| start_time | --start-time | yyyy mm dd hh mm (int) | Download data starting at this time |
| end_time | --end-time | yyyy mm dd hh mm (int) | Download data until this time |
| around_time | --around-time | yyyy mm dd hh mm hr_del (int) | Download data around this time, with the length of the window specified by an hour delta |

</details>

<details>
<summary>download_eumetsat_from_cart.py</summary>

A script to download products from the EUMETSAT data store, via a saved cart file. 

Input → Saved "cart" file(s) from data.eumetsat.int, in XML format.

In order to get a saved cart file, the following steps must be taken:
1. Navigate to the desired product, and save the desired product times to your "cart".
2. Navigate to the cart page, and click "download cart".
3. Save the contents of the cart to an XML file.

Output → Files downloaded from EUMETSAT of the selected product, filtered by time.

Options

| Option | Variable/Flag | Input | Description | 
|---|---|---|---|
| input_file_directory | _(required)_ | path | Directory to read input files from |
| save_directory | _(required)_ | path | Directory to save downloaded files to |
| api_key | _(required)_ | api_key (string) | An API key to be used to download the files |

> [!IMPORTANT]
> An API key is required to download files using this script. It can be found at the following 
> page after making an account: https://api.eumetsat.int/api-key/. API keys are only valid for one hour.

</details>

## Plotting Scripts

<details>
<summary>plot_goes.py</summary>

A script to plot products from GOES satellites, including single-band
channels and certain composite products.

Input → Multi-band GOES data files ("ABI-MCMIPC")
	   (Default file type downloaded from AWS script)

Output → A series of images that plot the selected product, one image
    per input file.

Options

| Option | Variable/Flag | Input | Description | 
|---|---|---|---|
| input_file_directory | _(required)_ | path | Directory to read input files from |
| save_directory | _(required)_ | path | Directory to save output to |
| band | -b, --band | band no. (int) | Specify a band to plot |
| composite | -c, --composite | product (String) | Specify a composite product to plot |
| bbox | --bbox | N S E W (int) | Specify the bounding box for the image in decimal lat-lon |
| center | --center | lat lon rad (float) | Automatically center the bounding box of the image around a point, with a given radius |
| points_from_file | --points-from-file | path | Specify a CSV file to read in points from, to be annotated on the plot. CSV must be in format: lat,lon,marker,label |
| pallete | -p, --pallete | pallete (String) | Specify the color pallete to be used for the image (*single band only*)
| show_colorbar | -cb, --show-colorbar | _None_ | Flag to turn on display of colorbar, for either albedo or brightness temperature (*single band only*)
| pixel_value | -pv, --pixel-value | _None_ | Flag to turn on display of pixel values, for either albedo or brightness temperature (*single band only*)

Notes
- A satellite band **or** a composite product must be specified, not both.
</details>

<details>
<summary>plot_meteosat.py</summary>

A script to plot products from GOES satellites, including single-band
channels and certain composite products.

Input → Multi-band GOES data files ("ABI-MCMIPC")
	   (Default file type downloaded from AWS script)

Output → A series of images that plot the selected product, one image
    per input file.

Options

| Option | Variable/Flag | Input | Description | 
|---|---|---|---|
| input_file_directory | _(required)_ | path | Directory to read input files from |
| save_directory | _(required)_ | path | Directory to save output to |
| band | -b, --band | band no. (int) | Specify a band to plot |
| composite | -c, --composite | product (String) | Specify a composite product to plot |
| bbox | --bbox | N S E W (int) | Specify the bounding box for the image in decimal lat-lon |
| center | --center | lat lon rad (float) | Automatically center the bounding box of the image around a point, with a given radius |
| points_from_file | --points-from-file | path | Specify a CSV file to read in points from, to be annotated on the plot. CSV must be in format: lat,lon,marker,label |
| pallete | -p, --pallete | pallete (String) | Specify the color pallete to be used for the image (*single band only*)
| show_colorbar | -cb, --show-colorbar | _None_ | Flag to turn on display of colorbar, for either albedo or brightness temperature (*single band only*)
| pixel_value | -pv, --pixel-value | _None_ | Flag to turn on display of pixel values, for either albedo or brightness temperature (*single band only*)

Notes
- A satellite band **or** a composite product must be specified, not both.
</details>

<details>
<summary>plot_himawari.py</summary>

A script to plot products from GOES satellites, including single-band
channels and certain composite products.

Input → Multi-band GOES data files ("ABI-MCMIPC")
	   (Default file type downloaded from AWS script)

Output → A series of images that plot the selected product, one image
    per input file.

Options

| Option | Variable/Flag | Input | Description | 
|---|---|---|---|
| input_file_directory | _(required)_ | path | Directory to read input files from |
| save_directory | _(required)_ | path | Directory to save output to |
| band | -b, --band | band no. (int) | Specify a band to plot |
| composite | -c, --composite | product (String) | Specify a composite product to plot |
| bbox | --bbox | N S E W (int) | Specify the bounding box for the image in decimal lat-lon |
| center | --center | lat lon rad (float) | Automatically center the bounding box of the image around a point, with a given radius |
| points_from_file | --points-from-file | path | Specify a CSV file to read in points from, to be annotated on the plot. CSV must be in format: lat,lon,marker,label |
| pallete | -p, --pallete | pallete (String) | Specify the color pallete to be used for the image (*single band only*)
| show_colorbar | -cb, --show-colorbar | _None_ | Flag to turn on display of colorbar, for either albedo or brightness temperature (*single band only*)
| pixel_value | -pv, --pixel-value | _None_ | Flag to turn on display of pixel values, for either albedo or brightness temperature (*single band only*)

Notes
- A satellite band **or** a composite product must be specified, not both.
</details>

<details>
<summary>plot_plan_view.py</summary>

A script to plot products from GOES satellites, including single-band
channels and certain composite products.

Input → Multi-band GOES data files ("ABI-MCMIPC")
	   (Default file type downloaded from AWS script)

Output → A series of images that plot the selected product, one image
    per input file.

Options

| Option | Variable/Flag | Input | Description | 
|---|---|---|---|
| input_file_directory | _(required)_ | path | Directory to read input files from |
| save_directory | _(required)_ | path | Directory to save output to |
| band | -b, --band | band no. (int) | Specify a band to plot |
| composite | -c, --composite | product (String) | Specify a composite product to plot |
| bbox | --bbox | N S E W (int) | Specify the bounding box for the image in decimal lat-lon |
| center | --center | lat lon rad (float) | Automatically center the bounding box of the image around a point, with a given radius |
| points_from_file | --points-from-file | path | Specify a CSV file to read in points from, to be annotated on the plot. CSV must be in format: lat,lon,marker,label |
| pallete | -p, --pallete | pallete (String) | Specify the color pallete to be used for the image (*single band only*)
| show_colorbar | -cb, --show-colorbar | _None_ | Flag to turn on display of colorbar, for either albedo or brightness temperature (*single band only*)
| pixel_value | -pv, --pixel-value | _None_ | Flag to turn on display of pixel values, for either albedo or brightness temperature (*single band only*)

Notes
- A satellite band **or** a composite product must be specified, not both.
</details>

<details>
<summary>plot_radar.py</summary>

A script to plot products from GOES satellites, including single-band
channels and certain composite products.

Input → Multi-band GOES data files ("ABI-MCMIPC")
	   (Default file type downloaded from AWS script)

Output → A series of images that plot the selected product, one image
    per input file.

Options

| Option | Variable/Flag | Input | Description | 
|---|---|---|---|
| input_file_directory | _(required)_ | path | Directory to read input files from |
| save_directory | _(required)_ | path | Directory to save output to |
| band | -b, --band | band no. (int) | Specify a band to plot |
| composite | -c, --composite | product (String) | Specify a composite product to plot |
| bbox | --bbox | N S E W (int) | Specify the bounding box for the image in decimal lat-lon |
| center | --center | lat lon rad (float) | Automatically center the bounding box of the image around a point, with a given radius |
| points_from_file | --points-from-file | path | Specify a CSV file to read in points from, to be annotated on the plot. CSV must be in format: lat,lon,marker,label |
| pallete | -p, --pallete | pallete (String) | Specify the color pallete to be used for the image (*single band only*)
| show_colorbar | -cb, --show-colorbar | _None_ | Flag to turn on display of colorbar, for either albedo or brightness temperature (*single band only*)
| pixel_value | -pv, --pixel-value | _None_ | Flag to turn on display of pixel values, for either albedo or brightness temperature (*single band only*)

Notes
- A satellite band **or** a composite product must be specified, not both.
</details>

<details>
<summary>plot_glm.py</summary>

A script to plot data from the Geostationary Lightning Mapper (GLM) located on GOES satellites, 
for either single time-steps or cumulative time-steps.

Input → GLM data files ("GLM-L2-LCFA")
     (Default file type downloaded from AWS script)

Output → A series of images that plot the selected product, one image
    per input file.

Options

| Option | Variable/Flag | Input | Description | 
|---|---|---|---|
| input_file_directory | _(required)_ | path | Directory to read input files from |
| save_directory | _(required)_ | path | Directory to save output to |
| single | -s, --single | _None_ | Specify plotting of flashes at single time-steps only (i.e. every 20 sec.) |
| composite | -c, --composite | minutes (int) | Specify plotting of all flashes in a rolling time period (default is 20 min)|
| bbox | --bbox | N S E W (float) | Specify the bounding box for the image in decimal lat-lon |
| center | --center | lat lon rad (float) | Automatically center the bounding box of the image around a point, with a given radius |
| points_from_file | --points-from-file | path | Specify a CSV file to read in points from, to be annotated on the plot. CSV must be in format: lat,lon,marker,label |
| save_kmz | --save-kmz | _None_ | Flag to save all detected flashes in a KMZ |

</details>

<details>
<summary>plot_airmet.py</summary>

A script to create a movie from a series of images, in mp4 format.

Input → A series of images

Output → An mp4 video at a specified framerate

Options

| Option | Variable/Flag | Input | Description | 
|---|---|---|---|
| input_file_directory | _(required)_ | path | Directory to read input files from |
| save_directory | _(required)_ | path | Directory to save output to |
| title | -t, --title | title (string) | Specify the title of the movie (default is "Movie[.mp4]") |
| fps | --fps | fps (int) | Specify the frames per second of the video |

</details>


## Creating/Transforming Scripts

<details>
<summary>create_raob_sounding.py</summary>

A script to create a series of "soundings" (vertical profiles) of data from a model,
specifically for the program RAOB (in RAOB CSV format).

Input → Model data files

> [!Note]
> Currently only ERA5 and HRRR analysis files are supported.

Output → A sounding or soundings, depending on the processing mode selected.

Options

| Option | Variable/Flag | Input | Description | 
|---|---|---|---|
| input_file_directory | _(required)_ | path | Directory to read input files from |
| save_directory | _(required)_ | path | Directory to save output to |
| model | _(required)_ | _hrrr_, _gfs_, _era5_,  or _wrf_ | Specify which model is being used as input data |
| wrf_input_freq | --wrf-input-freq | freq_str (str) | If using WRF data, specify the temporal frequency of the data in Pandas "freq" strings. This should be the same as the output timing for that domain.|
| point_sounding | -pt, --point-sounding | lat lon year month day hour minute (float) | **Processing Mode 1:** Automatically select the file closest in time to the given time, and make a sounding at the given coordinates. Only creates a single sounding file. |
| time_height | -th, --time-height | lat lon (float) | **Processing Mode 2:** Iterate through all data files, and produce a sounding at a given point for each file. These soundings can then be directly imported into RAOB to make a time-height diagram. |
| cross_section | -cs, --cross-section | path | **Processing Mode 3:** Using a CSV file containing points and times, iterate through each entry and find the file closest to that time, and make a sounding at that given point. These soundings can then be directly imported into RAOB to make a cross-section diagram |
| skip_duplicates | -d, --skip-duplicates | _None_ | Flag to turn on skipping of creating multiple soundings at the same point. Only useful for cross-section processing mode because RAOB doesn't like multiple soundings at the same point |
| verification_plot | -vp, --verification-plot | _None_ | Flag to generate a plot that can be used to verify the location of the soundings generated |
| grid_points | -gp, --grid-points | _None_ | Flag to save all grid points that are used by the script in a KMZ |

</details>

<details>
<summary>create_meteogram_data.py</summary>

A script to extract surface-level data from a series of data files,
and place that data into a CSV for plotting as a meteogram. Each row
of the CSV has the variables for the point, and a time (which corresponds
to the valid time of the file.)

Input → Model data file(s), in either Grib/Grib2 or NetCDF format 
		(Currently only HRRR model supported)

Output → A CSV containing surface-level data from the specified point,
		for all of the valid times contained in the data files (not
		necissarily ordered). A location verification plot can also be
		produced if desired.

Options

| Option | Variable/Flag | Input | Description | 
|---|---|---|---|
| model | _(required)_ | hrrr, gfs, era5, wrf | Specify data from what model is to be used. Will error or produce undefined result if data format differs from specifed model |
| input_file_directory | _(required)_ | path | Directory to read input files from |
| save_directory | _(required)_ | path | Directory to save output to |
| location | -l, --loc | lat, lon | Specify a point to take the meteogram data from |
| temp_units | -tu, --temperature-units | unit (string) | Specify the unit to be used for temperature/dewpoint |
| wind_units | -wu, --wind-units | unit (string) | Specify the unit to be used for wind speed/gusts |
| pressure_units | -pu, --pressure-units | unit (String) | Specify the unit to be used for pressure |
| verification_plot | -vp, --verification-plot | _None_ | Generate a plot to verify location of meteogram |
| precip_type | -pt, --record-precip-type | _None_ | Record precip type (if available) as a string (*Not Yet Implemented*)
</details>

<details>
<summary>create_raob_sounding.py</summary>

A script to extract surface-level data from a series of data files,
and place that data into a CSV for plotting as a meteogram. Each row
of the CSV has the variables for the point, and a time (which corresponds
to the valid time of the file.)

Input → Model data file(s), in either Grib/Grib2 or NetCDF format 
		(Currently only HRRR model supported)

Output → A CSV containing surface-level data from the specified point,
		for all of the valid times contained in the data files (not
		necissarily ordered). A location verification plot can also be
		produced if desired.

Options

| Option | Variable/Flag | Input | Description | 
|---|---|---|---|
| model | _(required)_ | hrrr, gfs, era5, wrf | Specify data from what model is to be used. Will error or produce undefined result if data format differs from specifed model |
| input_file_directory | _(required)_ | path | Directory to read input files from |
| save_directory | _(required)_ | path | Directory to save output to |
| location | -l, --loc | lat, lon | Specify a point to take the meteogram data from |
| temp_units | -tu, --temperature-units | unit (string) | Specify the unit to be used for temperature/dewpoint |
| wind_units | -wu, --wind-units | unit (string) | Specify the unit to be used for wind speed/gusts |
| pressure_units | -pu, --pressure-units | unit (String) | Specify the unit to be used for pressure |
| verification_plot | -vp, --verification-plot | _None_ | Generate a plot to verify location of meteogram |
| precip_type | -pt, --record-precip-type | _None_ | Record precip type (if available) as a string (*Not Yet Implemented*)
</details>

<details>
<summary>create_sun_moon_pos_table.py</summary>

A script to extract surface-level data from a series of data files,
and place that data into a CSV for plotting as a meteogram. Each row
of the CSV has the variables for the point, and a time (which corresponds
to the valid time of the file.)

Input → Model data file(s), in either Grib/Grib2 or NetCDF format 
		(Currently only HRRR model supported)

Output → A CSV containing surface-level data from the specified point,
		for all of the valid times contained in the data files (not
		necissarily ordered). A location verification plot can also be
		produced if desired.

Options

| Option | Variable/Flag | Input | Description | 
|---|---|---|---|
| model | _(required)_ | hrrr, gfs, era5, wrf | Specify data from what model is to be used. Will error or produce undefined result if data format differs from specifed model |
| input_file_directory | _(required)_ | path | Directory to read input files from |
| save_directory | _(required)_ | path | Directory to save output to |
| location | -l, --loc | lat, lon | Specify a point to take the meteogram data from |
| temp_units | -tu, --temperature-units | unit (string) | Specify the unit to be used for temperature/dewpoint |
| wind_units | -wu, --wind-units | unit (string) | Specify the unit to be used for wind speed/gusts |
| pressure_units | -pu, --pressure-units | unit (String) | Specify the unit to be used for pressure |
| verification_plot | -vp, --verification-plot | _None_ | Generate a plot to verify location of meteogram |
| precip_type | -pt, --record-precip-type | _None_ | Record precip type (if available) as a string (*Not Yet Implemented*)
</details>

<details>
<summary>create_twilight_times_listing.py</summary>

A script to extract surface-level data from a series of data files,
and place that data into a CSV for plotting as a meteogram. Each row
of the CSV has the variables for the point, and a time (which corresponds
to the valid time of the file.)

Input → Model data file(s), in either Grib/Grib2 or NetCDF format 
		(Currently only HRRR model supported)

Output → A CSV containing surface-level data from the specified point,
		for all of the valid times contained in the data files (not
		necissarily ordered). A location verification plot can also be
		produced if desired.

Options

| Option | Variable/Flag | Input | Description | 
|---|---|---|---|
| model | _(required)_ | hrrr, gfs, era5, wrf | Specify data from what model is to be used. Will error or produce undefined result if data format differs from specifed model |
| input_file_directory | _(required)_ | path | Directory to read input files from |
| save_directory | _(required)_ | path | Directory to save output to |
| location | -l, --loc | lat, lon | Specify a point to take the meteogram data from |
| temp_units | -tu, --temperature-units | unit (string) | Specify the unit to be used for temperature/dewpoint |
| wind_units | -wu, --wind-units | unit (string) | Specify the unit to be used for wind speed/gusts |
| pressure_units | -pu, --pressure-units | unit (String) | Specify the unit to be used for pressure |
| verification_plot | -vp, --verification-plot | _None_ | Generate a plot to verify location of meteogram |
| precip_type | -pt, --record-precip-type | _None_ | Record precip type (if available) as a string (*Not Yet Implemented*)
</details>

<details>
<summary>create_movie.py</summary>

A script to create a movie from a series of images, in mp4 format.

Input → A series of images

Output → An mp4 video at a specified framerate

Options

| Option | Variable/Flag | Input | Description | 
|---|---|---|---|
| input_file_directory | _(required)_ | path | Directory to read input files from |
| save_directory | _(required)_ | path | Directory to save output to |
| title | -t, --title | title (string) | Specify the title of the movie (default is "Movie[.mp4]") |
| fps | --fps | fps (int) | Specify the frames per second of the video |

</details>
