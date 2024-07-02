from datetime import datetime
import requests
import json
import os
import glob
from threading import RLock as TRLock
from concurrent import futures
from concurrent.futures import ThreadPoolExecutor
from functools import partial
import argparse

import bs4 as soup
import lxml
from tqdm.auto import tqdm
from tqdm.contrib.concurrent import process_map, thread_map

def requests_download_multithread_worker(save_dir: str,
                                         api_token: str,
                                         file_url: str,
                                         file_name: str,
                                         progress_pos: int) -> None:

    try:
        dl_obj = requests.get(f'{file_url}?access_token={api_token}', stream=True)
    except Exception as err:
        print(err)
    file_size = int(dl_obj.headers.get("content-length", 0))

    dest_path = os.path.join(save_dir, file_name)

    with tqdm(unit='B', unit_scale=True, unit_divisor=1024, miniters=1, desc=f'{file_name}', total=file_size, leave=None) as progress:
        with open(dest_path, 'wb') as out_file:
            for chunk in dl_obj.iter_content(chunk_size=1024):
                progress.update(len(chunk))
                out_file.write(chunk)

def requests_download_multithread(save_dir: str,
                                  api_token: str,
                                  file_names: list[str, ...],
                                  file_urls: list[str, ...]) -> None:

    tqdm.set_lock(TRLock())
    try:
        with ThreadPoolExecutor(initializer=tqdm.set_lock, initargs=(tqdm.get_lock(),)) as executor:
            executor.map(partial(requests_download_multithread_worker, save_dir, api_token), file_urls, file_names, range(1, len(file_urls)+1, 1))
    except Exception as err:
        print(err)

def read_details_from_cart_file(cart_file: str) -> tuple[list[str, ...], list[str, ...]]:

    try:
        with open(cart_file, 'rb') as raw_file:
            file = soup.BeautifulSoup(raw_file, "xml")
    except Exception as err:
        print(f'Error opening cart file: {err}')

    file_urls = []
    file_names = []

    for item in file.find_all('file'):
        file_dest = item['name']
        file_url = item.url.text

        file_urls += [file_url]
        file_names += [file_dest]

    return file_urls, file_names

if __name__ == "__main__":

    parser = argparse.ArgumentParser(description=f'Download EUMETSAT/Meteosat Data via EUMETSAT Data Store Cart Files')

    parser.add_argument('api_token',  
                        help='a valid API token from the EUMETSAT data store (https://api.eumetsat.int/api-key/)',
                        type=str)
    parser.add_argument('input_file_directory',  
                        help='directory to read input cart files from',
                        type=str)
    parser.add_argument('save_directory',  
                        help='directory to download files to',
                        type=str)

    args = parser.parse_args()

    print(args)

    if not args.api_token:
        raise ValueError('A valid API token must be specified.')

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

    input_files = sorted(glob.glob(f'{input_dir}*.xml'))
    num_input_files = len(input_files)

    tot_files = 0

    for cart_file in input_files:

        file_urls, file_names = read_details_from_cart_file(cart_file)
        print(f"Downloading {len(file_urls)} files from EUMETSAT cart file ({cart_file})...")

        requests_download_multithread(save_dir, args.api_token, file_names, file_urls)

