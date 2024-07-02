##             plot_goes.py              ##

A script to plot GOES data.

Input: Multi-band GOES data files (ABI-MCMIPC)
	   (Default file type downloaded from AWS script)

Output: A series of 

##           meteogram_data.py           ##

A script to extract surface-level data from a series of data files,
and place that data into a CSV for plotting as a meteogram. Each row
of the CSV has the variables for the point, and a time (which corresponds
to the valid time of the file.)

Input 
: Model data file(s), in either Grib/Grib2 or NetCDF format 
		(Currently only HRRR model supported)

Output 
: A CSV containing surface-level data from the specified point,
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


positional arguments:
  {hrrr,gfs,era5,wrf}   specify data from what model is to be used. will error or produce undefined result if data format differs from
                        specifed model
  input_file_directory  directory to read input files from
  save_directory        directory to save meteogram to

options:
  -h, --help            show this help message and exit
  -l lat lon, --loc lat lon
                        specify a point to take the meteogram data from
  -tu TEMP_UNITS, --temperature-units TEMP_UNITS
                        specify the unit to be used for temperature/dewpoint (default: degC)
  -wu WIND_UNITS, --wind-units WIND_UNITS
                        specify the unit to be used for wind speed/gusts (default: kt)
  -pu PRESSURE_UNITS, --pressure-units PRESSURE_UNITS
                        specify the unit to be used for pressure (default: hectopascal)
  -vp, --verification-plot
                        generate a plot to verify location of meteogram
  -pt, --record-precip-type
                        record precip type as a string (default: off)
  -t name, --title name
                        specify a title for the output data csv
