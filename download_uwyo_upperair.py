import re
import requests
import os
import glob
import argparse
from time import sleep

import pandas as pd
from tqdm.auto import tqdm

SLEEP_TIME = 5
PLOT_FREQ = '12h'

REGION_READABLE = {
	'naconf': 'NorthAmerica',
	'us30': 'CONUS',
	'samer': 'SouthAmerica',
	'pac': 'Pacific',
	'npac': 'NorthPacific',
	'ant': 'Antarctica',
	'np': 'Arctic',
	'europe': 'Europe',
	'africa': 'Africa',
	'seasia': 'SEAsia',
	'mideast': 'MENA',
	'indocn': 'IndianOcean',
	'nh': 'NorthernHem',
	'sh': 'SouthernHem'
}

LEVEL_READABLE = {
	'10': '10hPa',
	'30': '30hPa',
	'50': '50hPa',
	'70': '70hPa',
	'100': '100hPa',
	'150': '150hPa',
	'200': '200hPa',
	'250': '250hPa',
	'300': '300hPa',
	'400': '400hPa',
	'500': '500hPa',
	'700': '700hPa',
	'850': '850hPa',
	'1000': '1000hPa',
	'0': 'Surface',
	'T': 'Thickness',
	'K': 'LiftedIDX',
	'W': 'PWat',
	'A': '500hPaVort'
}

def download_plot(timestamp: pd.Timestamp,
				  user_params: dict,
				  save_directory: str,
				  use_desc_name: bool) -> None:

	date_str = f'{str(timestamp.year)}-{str(timestamp.month).zfill(2)}-{str(timestamp.day).zfill(2)}'

	date_params = {'date': date_str,
				   'hour': str(timestamp.hour),}

	user_params.update(date_params)

	# print(user_params)

	gen_req = requests.get('http://www.weather.uwyo.edu/cgi-bin/uamap', 
						   params=user_params)

	# print(gen_req.url)
	# print(gen_req.content)

	img_loc = re.search(r"\/upperair\/maps\/[\s\S]+\.gif", 
						gen_req.text.replace('\n',''))


	img_req = requests.get(f'http://www.weather.uwyo.edu{img_loc[0]}')
	# print(img_req.url)


	if use_desc_name:
		file_name = f"UWYO_{REGION_READABLE[user_params['REGION']]}_{LEVEL_READABLE[user_params['LEVEL']]}_{str(timestamp.year)}{str(timestamp.month).zfill(2)}{str(timestamp.day).zfill(2)}{str(timestamp.hour).zfill(2)}UTC.gif"
	else:
		file_name = f"{str(timestamp.year)}{str(timestamp.month).zfill(2)}{str(timestamp.day).zfill(2)}{str(timestamp.hour).zfill(2)}.{user_params['LEVEL']}oa.{user_params['REGION']}.gif"

	if not os.path.exists(save_dir):
		os.makedirs(save_dir)

	with open(os.path.join(save_dir, file_name), 'wb') as img:
		img.write(img_req.content)
		img.close()


if __name__ == "__main__":

	parser = argparse.ArgumentParser(description=f'Download upper-air and surface charts from University of Wyoming')
	parser.add_argument('--region', 
						help='specify a region to plot', 
						type=str,
						choices=list(REGION_READABLE.keys()),
						dest='region',
						metavar=('region'),
						default=None)
	parser.add_argument('--level', 
						help='specify a plot level or other variable', 
						type=str,
						choices=list(LEVEL_READABLE.keys()),
						dest='level',
						metavar=('level_id'),
						default=None)
	parser.add_argument('--remove-obs', 
                        help='remove observations on plots', 
                        action='store_true',
                        dest='remove_obs',
                        default=False)
	parser.add_argument('--descriptive-name', 
                        help='use more descriptive name', 
                        action='store_true',
                        dest='desc_name',
                        default=False)
	parser.add_argument('--start-time',  
						help='download plots starting at this time',
						nargs=4,
						type=int,
						dest='start_time',
						metavar=('yyyy', 'mm', 'dd', 'hh'),
						default=None)
	parser.add_argument('--end-time',  
						help='download plots ending at this time',
						nargs=4,
						type=int,
						dest='end_time',
						metavar=('yyyy', 'mm', 'dd', 'hh'),
						default=None)
	parser.add_argument('save_directory',  
						help='directory to save plots to',
						type=str)

	args = parser.parse_args()
	print(args)

	## Quality Checks/Sanitizing of input arguments

	start_ts = pd.Timestamp(args.start_time[0], 
							args.start_time[1], 
							args.start_time[2], 
							args.start_time[3])

	end_ts = pd.Timestamp(args.end_time[0], 
						  args.end_time[1], 
						  args.end_time[2], 
						  args.end_time[3])

	dl_range = pd.period_range(start_ts.to_period(PLOT_FREQ), 
							   end_ts.to_period(PLOT_FREQ),
							   freq=PLOT_FREQ)
	# print(dl_range)

	if args.save_directory[-1:] != "/":
		save_dir = args.save_directory + "/"
	else:
		save_dir = args.save_directory

	user_settings = dict()

	if args.region:
		user_settings.update({'REGION': args.region})
	else:
		user_settings.update({'REGION': 'us30'}) #conus

	user_settings.update({'OUTPUT': 'gif'})

	if args.remove_obs:
		user_settings.update({'TYPE': 'an'})
	else:
		user_settings.update({'TYPE': ['obs', 'an']})
		
	if args.level:
		user_settings.update({'LEVEL': args.level})
	else:
		user_settings.update({'LEVEL': '0'}) #sfc

	print(f"Saving to '{save_dir}'")

	with tqdm(miniters=0, total=len(dl_range), desc=f"Generating and downloading {len(dl_range)} plots...", ascii=" 123456789#") as progress:
	   for time in dl_range:
	   		this_ts = time.to_timestamp()

	   		download_plot(this_ts,
	   					  user_settings,
	   					  save_dir,
	   					  True if args.desc_name else False)

	   		# sleep(SLEEP_TIME) #prevents possible timeouts?

	   		progress.update()







