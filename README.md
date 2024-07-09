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
|---|---|---|---|
| band | -b, --band | band no. (int) | Specify a band to plot |
| composite | -c, --composite | product (String) | Specify a composite product to plot |
|---|---|---|---|
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
<summary>meteogram_data.py</summary>

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
|---|---|---|---|
| location | -l, --loc | lat, lon | Specify a point to take the meteogram data from |
| temp_units | -tu, --temperature-units | unit (string) | Specify the unit to be used for temperature/dewpoint |
| wind_units | -wu, --wind-units | unit (string) | Specify the unit to be used for wind speed/gusts |
| pressure_units | -pu, --pressure-units | unit (String) | Specify the unit to be used for pressure |
| verification_plot | -vp, --verification-plot | _None_ | Generate a plot to verify location of meteogram |
| precip_type | -pt, --record-precip-type | _None_ | Record precip type (if available) as a string (*Not Yet Implemented*)
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
|---|---|---|---|
| single | -s, --single | _None_ | Specify plotting of flashes at single time-steps only (i.e. every 20 sec.) |
| composite | -c, --composite | minutes (int) | Specify plotting of all flashes in a rolling time period (default is 20 min)|
|---|---|---|---|
| bbox | --bbox | N S E W (float) | Specify the bounding box for the image in decimal lat-lon |
| center | --center | lat lon rad (float) | Automatically center the bounding box of the image around a point, with a given radius |
| points_from_file | --points-from-file | path | Specify a CSV file to read in points from, to be annotated on the plot. CSV must be in format: lat,lon,marker,label |
| save_kmz | --save-kmz | _None_ | Flag to save all detected flashes in a KMZ |

</details>

<details>
<summary>create_raob_sounding.py</summary>

A script to create a series of "soundings" (vertical profiles) of data from a model,
specifically for the program RAOB (in RAOB CSV format).

Input → Model data files

>![WARNING]
>Currently only ERA5 and HRRR analysis files are supported.

Output → A sounding or soundings, depending on the processing mode selected.

Options

| Option | Variable/Flag | Input | Description | 
|---|---|---|---|
| input_file_directory | _(required)_ | path | Directory to read input files from |
| save_directory | _(required)_ | path | Directory to save output to |
| model | _(required)_ | _hrrr_, _gfs_, _era5_,  or _wrf_ | Specify which model is being used as input data |
| wrf_input_freq | --wrf-input-freq | freq_str (str) | If using WRF data, specify the temporal frequency of the data in Pandas "freq" strings. This should be the same as the output timing for that domain.|
|---|---|---|---|
| point_sounding | -pt, --point-sounding | lat lon year month day hour minute (float) | **Processing Mode 1:** Automatically select the file closest in time to the given time, and make a sounding at the given coordinates. Only creates a single sounding file. |
| time_height | -th, --time-height | lat lon (float) | **Processing Mode 2:** Iterate through all data files, and produce a sounding at a given point for each file. These soundings can then be directly imported into RAOB to make a time-height diagram. |
| cross_section | -cs, --cross-section | path | **Processing Mode 3:** Using a CSV file containing points and times, iterate through each entry and find the file closest to that time, and make a sounding at that given point. These soundings can then be directly imported into RAOB to make a cross-section diagram |
|---|---|---|---|
| skip_duplicates | -d, --skip-duplicates | _None_ | Flag to turn on skipping of creating multiple soundings at the same point. Only useful for cross-section processing mode because RAOB doesn't like multiple soundings at the same point |
| verification_plot | -vp, --verification-plot | _None_ | Flag to generate a plot that can be used to verify the location of the soundings generated |
| grid_points | -gp, --grid-points | _None_ | Flag to save all grid points that are used by the script in a KMZ |

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
|---|---|---|---|
| title | -t, --title | title (string) | Specify the title of the movie (default is "Movie[.mp4]") |
| fps | --fps | fps (int) | Specify the frames per second of the video |

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
|---|---|---|---|
| title | -t, --title | title (string) | Specify the title of the movie (default is "Movie[.mp4]") |
| fps | --fps | fps (int) | Specify the frames per second of the video |

</details>

