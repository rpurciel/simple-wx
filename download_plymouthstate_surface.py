import re
import requests
import os
import glob
import argparse
from time import sleep

import pandas as pd
from tqdm.auto import tqdm

SLEEP_TIME = 5
PLOT_FREQ = '1h'

STATIC_PARAMS = {'density': '',
				 'sc': '1.0',
				 'ge': '912x650',
				 'pg': 'print',
				 'id': '',
				 'zoom': '.6',}

def download_plot(timestamp: pd.Timestamp,
				  user_params: dict,
				  save_directory: str) -> None:

	date_params = {'yy': str(timestamp.year)[2:],
				   'mm': str(timestamp.month).zfill(2),
				   'dd': str(timestamp.day).zfill(2),
				   'hh': str(timestamp.hour).zfill(2),}

	user_params.update(date_params)
	user_params.update(STATIC_PARAMS)

	# print(user_params)

	gen_req = requests.get('https://vortex.plymouth.edu/wxp/cgi-bin/sfc/gen-pltmap-a.cgi', 
						   params=user_params)

	# print(gen_req.url)
	# print(gen_req.content)

	img_loc = re.search(r"\/myowxp\/maps\/\d+.gif", 
						gen_req.text.replace('\n',''))


	img_req = requests.get(f'https://vortex.plymouth.edu{img_loc[0]}')
	# print(img_req.url)

	file_name = f'PlymouthState_Sfc_{user_params["va"]}_{user_params["re"].upper()}_{str(timestamp.year)}{str(timestamp.month).zfill(2)}{str(timestamp.day).zfill(2)}{str(timestamp.hour).zfill(2)}UTC.gif'

	if not os.path.exists(save_dir):
		os.makedirs(save_dir)

	with open(os.path.join(save_dir, file_name), 'wb') as img:
		img.write(img_req.content)
		img.close()


if __name__ == "__main__":

	parser = argparse.ArgumentParser(description=f'Download surface charts from plymouth state')
	parser.add_argument('--region', 
						help='specify a point to take the listing data from', 
						type=str,
						dest='reg',
						metavar=('reg'),
						default=None)
	parser.add_argument('--variable', 
						help='specify an elevation in meters at the observation point', 
						type=str,
						dest='var',
						metavar=('var_id'),
						default=None)
	parser.add_argument('--start-time',  
						help='download plots  starting at this time',
						nargs=4,
						type=int,
						dest='start_time',
						metavar=('yyyy', 'mm', 'dd', 'hh'),
						default=None)
	parser.add_argument('--end-time',  
						help='download plots starting at this time',
						nargs=4,
						type=int,
						dest='end_time',
						metavar=('yyyy', 'mm', 'dd', 'hh'),
						default=None)
	parser.add_argument('save_directory',  
						help='directory to save listing to',
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
							   end_ts.to_period(PLOT_FREQ))

	if args.save_directory[-1:] != "/":
		save_dir = args.save_directory + "/"
	else:
		save_dir = args.save_directory

	user_settings = dict()

	if args.reg:
		user_settings.update({'re': args.reg})
	else:
		user_settings.update({'re': 'us'})
		
	if args.var:
		user_settings.update({'va': args.var})
	else:
		user_settings.update({'va': 'temp'})

	print(f"Saving to '{save_dir}'")

	with tqdm(miniters=0, total=len(dl_range), desc=f"Generating and downloading {len(dl_range)} plots...", ascii=" 123456789#") as progress:
	   for time in dl_range:
	   		this_ts = time.to_timestamp()

	   		download_plot(this_ts,
	   					  user_settings,
	   					  save_dir)

	   		# sleep(SLEEP_TIME) #prevents possible timeouts?

	   		progress.update()







