###TODO: Smart logging for sub functions of download
###      Making everything into classes might make everything easier?

import requests
import sys
import os
import warnings
from pathlib import Path, PurePath
import math
import re
from datetime import datetime
import json
import argparse
import glob

import pandas as pd
import simplekml
from tqdm.auto import tqdm

warnings.filterwarnings("ignore")

DEF_STATE_FILTER = False

DEF_CARDINAL_DIR_TO_DEG_DICT = {
    "N" : 0,
    "NNE" : 22.5,
    "NE" : 45,
    "ENE" : 67.5,
    "E" : 90,
    "ESE" : 112.5,
    "SE" : 135,
    "SSE" : 157.5,
    "S" : 180,
    "SSW" : 202.5,
    "SW" : 225,
    "WSW" : 247.5,
    "W" : 270,
    "WNW" : 292.5,
    "NW" : 315,
    "NNW" : 337.5,
}

DEF_AIRMET_TYPE_TO_COND_DICT = {
	"SIERRA" : "IFR",
	"TANGO" : "TURB",
	"ZULU" : "ICE"
}

DEF_ANCILLARY_PATH_TO_VORS_RELATIVE_TO_SRC = "/Users/rpurciel/WeatherExtreme Ltd Dropbox/Ryan Purciel/Scripts/Single Use/ancillary/vors.csv"


def plot_kmz(save_dir, 
             subgroups, 
             airmet_type, 
             airmet_id, 
             airmet_raw_text, 
             valid_time, 
             iss_time, 
             **kwargs):

    if not os.path.exists(save_dir):
        os.makedirs(save_dir)

    verbose = False
    debug = False

    airmet_kml = simplekml.Kml()

    iss_time_str = iss_time.strftime("%Y-%m-%d %H:%M:%S UTC")
    valid_time_str = valid_time.strftime("%Y-%m-%d %H:%M:%S UTC")

    cond_str = DEF_AIRMET_TYPE_TO_COND_DICT.get(airmet_type)

    airmet_kml.document.name = f"{airmet_id} [{cond_str}]"
    airmet_kml.document.description = f"{iss_time_str} THRU {valid_time_str}\n{airmet_raw_text}"

    state_filter = []
    for arg, value in kwargs.items():
        if arg == 'filter_by_states':
            state_filter = value

    num_polygons = 0

    for group in subgroups:
        if verbose:
            print(f"PLOTTER: Iterating through following group:\n{group}")
        quals = group.get("qualifiers")
        airmet_flag = False
        llws_pot_flag = False

        if not quals:
            quals = []

        for qual in quals:
            #print(f"QUAL: {qual}")
            if (qual.find("AIRMET") != -1) or (qual.find("LLWS") != -1) or (qual.find("OTLK") != -1) or (qual.find("SIGMET") != -1):
                if verbose:
                    print(f"PLOTTER: Group is an AIRMET, should be plotted.")
                airmet_flag = True
                airmet_title = f"AIRMET {airmet_type}"

                if qual.find("LLWS") != -1:
                    if verbose:
                        print(f"PLOTTER: Group is an LLWS Potential Group, should be plotted.")
                    llws_pot_flag = True
                    airmet_title = f"LLWS POTENTIAL"

                if qual.find("OTLK") != -1:
                    if verbose:
                        print(f"PLOTTER: Group is an Outlook Group, should be plotted.")
                    otlk_flag = True
                    airmet_title = f"OUTLOOK"
                    otlk_valid_time = qual

                if qual.find("OTLK") != -1:
                    if verbose:
                        print(f"PLOTTER: Group is an Convective Sigmet, should be plotted.")
                    airmet_flag = True
                    airmet_title = f"Convective Sigmet"

            else:
                if verbose:
                    print(f"PLOTTER: Group is an freezing level or unknown. Skipping plotting...")

        if airmet_flag:
            if state_filter:
                if verbose:
                    print(f"!!! PLOTTER: State filtering turned ON\nPLOTTER: Plotting AIRMETS that only intersect the following states: {state_filter}")

                states = group.get("states")
                includes_state_flag = False

                if not states:
                    states = []

                for state in states:
                    if state in selected_states:
                        includes_state_flag = True

                if includes_state_flag:
                    if verbose:
                        print("PLOTTER: AIRMET includes a specified state. Plotting...")        
                    vors = group.get("vors")
                    desc = group.get("desc")
                    airmet_for = quals[0]

                    if llws_pot_flag:
                        desc_text = "FOR LLWS POTENTIAL\n" + desc
                    elif otlk_flag:
                        desc_text = otlk_valid_time + "\n" + desc
                    elif (llws_pot_flag != True) and (otlk_flag != True):
                        desc_text = "FOR " + airmet_for[7:] + "\n" + desc

                    airmet_kml, status = _add_poly_to_kml(airmet_kml, vors, airmet_type, airmet_title, desc_text)
                    num_polygons += status

                else:
                    if verbose:
                        print("PLOTTER: AIRMET does not include a specified state. Skipping...")
            else:
                if verbose:
                    print("PLOTTER: Plotting AIRMET...")
                states = group.get("states")
                vors = group.get("vors")
                desc = group.get("desc")
                airmet_for = quals[0]
                desc_text = "FOR " + airmet_for[7:] + "\n" + desc
                if llws_pot_flag:
                    desc_text = "FOR LLWS POTENTIAL\n" + desc

                airmet_kml, status = _add_poly_to_kml(airmet_kml, vors, airmet_type, airmet_title, desc_text)
                num_polygons += status

    if num_polygons == 0:
        error_str = "NoPolygons"

        return 0, 0, error_str

    else:
        iss_time_str = iss_time.strftime("%Y%m%d_%H%M%S")
        valid_time_str = valid_time.strftime("%Y%m%d_%H%M%S")
        airmet_id_for_file = airmet_id.replace(" ", "")

        file_name = f"{airmet_id_for_file}_{cond_str}_iss{iss_time_str}_valid{valid_time_str}"

        dest_path = os.path.join(save_dir, file_name + ".kmz")

        airmet_kml.savekmz(dest_path)

        return 1, 0, dest_path
    
def _add_poly_to_kml(kml, list_of_vors, airmet_type, airmet_title, desc, **kwargs):

    verbose = False
    debug = False

    if debug:
        print("DEBUG: Kwargs passed:", kwargs)

    airmet_points = []

    for vor in list_of_vors:
        if debug:
            print(f"DEBUG: Selected vor {vor}")
        dir_args = []
        if vor.find(" ") != -1:
            if debug:
                print(f"DEBUG: Combination VOR/DIR, splitting...")
            dir_vor = vor.split(" ")
            dist_carddir = dir_vor[0]
            vor = dir_vor[1]

            number_match = re.match(r"^\d*", dist_carddir)
            dist_nm = dist_carddir[:number_match.end()]
            card_dir = dist_carddir[number_match.end():]

            if debug:
                print(f"DEBUG: VOR distance: {dist_nm}\nDEBUG: VOR cardinal direction: {card_dir}")

            dir_args.append(dist_nm.strip())
            dir_args.append(card_dir.strip())

        point_lat, point_lon = _vor_dir_to_lat_lon(vor, dir_args)

        if debug:
        	print(f"DEBUG: Selected point: ({point_lat}, {point_lon})")

        airmet_points.append((point_lon, point_lat))

    #airmet_points = airmet_points[:-1]

    if debug:
    	print(f"DEBUG: AIRMET points:")
    	print(airmet_points)

    if airmet_type == "SIERRA":
        polygon_line_color = simplekml.Color.violet
        polygon_line_width = 5
        polygon_color = simplekml.Color.changealphaint(60, simplekml.Color.violet)
    elif airmet_type == "ZULU":
        polygon_line_color = simplekml.Color.cyan
        polygon_line_width = 5
        polygon_color = simplekml.Color.changealphaint(60, simplekml.Color.cyan)
    elif airmet_type == "TANGO":
        polygon_line_color = simplekml.Color.coral
        polygon_line_width = 5
        polygon_color = simplekml.Color.changealphaint(60, simplekml.Color.coral)
    else:
        polygon_line_color = simplekml.Color.red
        polygon_line_width = 5
        polygon_color = simplekml.Color.changealphaint(60, simplekml.Color.red)

    airmet_poly = kml.newpolygon(name=airmet_title,
                                 description=desc,
                                 outerboundaryis=airmet_points,)

    airmet_poly.style.linestyle.color = polygon_line_color
    airmet_poly.style.linestyle.width = polygon_line_width
    airmet_poly.style.polystyle.color = polygon_color

    return kml, 1

def _vor_dir_to_lat_lon(vor, *args):

    complex_vor_flag = False

    if args != ([],):
        args = args[0]
        distance_nm = int(args[0])
        cardinal = args[1]
        bearing_deg = DEF_CARDINAL_DIR_TO_DEG_DICT.get(cardinal)
        complex_vor_flag = True


    vor_path = DEF_ANCILLARY_PATH_TO_VORS_RELATIVE_TO_SRC

    vors_master_list = pd.read_csv("/Users/rpurciel/WeatherExtreme Ltd Dropbox/Ryan Purciel/Scripts/Single Use/ancillary/vors.csv", sep=",", na_values = ["0","M"], index_col=0)

    try:
        vor_data = vors_master_list.loc[vor]
    except:
        print(f"ERROR: VOR data not found for {vor}. Please add an issue on Github with details of this issue.")

    #print(vor_data)
    vor_lat = vor_data.lat
    vor_lon = vor_data.lon
    #print(f"LAT: {vor_lat}\nLON: {vor_lon}")

    if complex_vor_flag:

        eq_rad_km = 6378.137
        pol_rad_km = 6356.752

        bearing_rad = math.radians(bearing_deg)
        distance_km = 1.852 * distance_nm

        vor_lat_rad = math.radians(vor_lat)
        vor_lon_rad = math.radians(vor_lon)

        #Copied from earlier code
        #Sorry reader

        #Radius of earth at latitude
        vor_lat_earth_radius = (((((eq_rad_km**2)*math.cos(vor_lat_rad))**2)
                              +(((pol_rad_km**2)*math.sin(vor_lat_rad))**2))
                              /((eq_rad_km*math.cos(vor_lat_rad))**2
                              +(pol_rad_km*math.sin(vor_lat_rad))**2))**0.5

        #Latitude of airmet point in radians
        airmet_pt_lat_rad = math.asin(math.sin(vor_lat_rad)*math.cos(distance_km/vor_lat_earth_radius) +
                             math.cos(vor_lat_rad)*math.sin(distance_km/vor_lat_earth_radius)*math.cos(bearing_rad))

        #Longitude of airmet point in radians
        airmet_pt_lon_rad = vor_lon_rad + math.atan2(math.sin(bearing_rad)*math.sin(distance_km/vor_lat_earth_radius)*math.cos(vor_lat_rad),
                             math.cos(distance_km/vor_lat_earth_radius)-math.sin(vor_lat_rad)*math.sin(airmet_pt_lat_rad))

        airmet_pt_lat = math.degrees(airmet_pt_lat_rad)
        airmet_pt_lon = math.degrees(airmet_pt_lon_rad)

    else:

        airmet_pt_lat = vor_lat
        airmet_pt_lon = vor_lon

    return airmet_pt_lat, airmet_pt_lon

def _header_to_dict(header, **kwargs):
    
    main_block_match = re.search(r"\#([A-Z]|\d){4}\s.+", header)
    
    header = header[main_block_match.start():]
    amended = False
    if header.find("AMD") != -1:
        header = header.replace(" AMD", "")
        amended = True
    
    airmet_id_match = re.search(r"^\#([A-Z]|\d){4}\s", header)
    airmet_airport_match = re.search(r"\#&([A-Z]|\d){4}\s", header)
    airmet_iss_time_match = re.search(r"\d{6}\#", header)
    airmet_type_match = re.search(r"\#AIRMET\s\w+", header)
    airmet_conds_match = re.search(r"FOR\s(\w|\s)+VALID", header)
    airmet_valid_time_match = re.search(r"VALID\sUNTIL\s\d{6}", header)
    
    airmet_id = header[airmet_id_match.start():airmet_id_match.end()].replace("#", "").strip()
    try:
        airmet_airport = header[airmet_airport_match.start():airmet_airport_match.end()].replace("#&", "").strip()
    except:
        airmet_airport = ''
    airmet_iss_time = header[airmet_iss_time_match.start():airmet_iss_time_match.end()].replace("#", "").strip()
    airmet_type = header[airmet_type_match.start():airmet_type_match.end()].replace("#AIRMET ", "").strip()
    airmet_valid_time = header[airmet_valid_time_match.start():airmet_valid_time_match.end()].replace("VALID UNTIL ", "").strip()
    
    if amended:
        airmet_id = airmet_id + " AMD"
    
    airmet_conds_str = header[airmet_conds_match.start():airmet_conds_match.end()].replace("FOR ", "").replace(" VALID", "")
    airmet_conds_str = airmet_conds_str.replace("STG WNDS", "STG_WNDS").replace("MTN OBSCN", "MTN_OBSCN").replace("STG SFC WNDS", "STD_SFC_WNDS").replace("AND ", "")
    airmet_conds = airmet_conds_str.split(" ")
    
    iss_day = int(airmet_iss_time[:2])
    iss_hour = int(airmet_iss_time[2:4])
    iss_minute = int(airmet_iss_time[4:6])
    
    valid_day = int(airmet_valid_time[:2])
    valid_hour = int(airmet_valid_time[2:4])
    valid_minute = int(airmet_valid_time[4:6])
    
    header_dict = {
        "airmet_id" : airmet_id,
        "iss_airport" : airmet_airport,
        "iss_year" : year,
        "iss_month" : month,
        "iss_day" : iss_day,
        "iss_hour" : iss_hour,
        "iss_minute" : iss_minute,
        "iss_time_str" : airmet_iss_time,
        "valid_year" : year,
        "valid_month" : month,
        "valid_day" : valid_day,
        "valid_hour" : valid_hour,
        "valid_minute" : valid_minute,
        "valid_time_str" : airmet_valid_time,
        "airmet_type" : airmet_type,
        "conditions" : airmet_conds,
    }
    
    return header_dict
    
if __name__ == "__main__":

    parser = argparse.ArgumentParser(description=f'Plot AIRMETS in a JSON file')
    parser.add_argument('input_file_directory',  
                        help='directory to read input files from',
                        type=str)
    parser.add_argument('save_directory',  
                        help='directory to save output KMZs to',
                        type=str)

    args = parser.parse_args()
    print(args)

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

    input_files = sorted(glob.glob(input_dir + "*"))

    with tqdm(miniters=0, total=len(input_files), desc=f'Plotting AIRMETS from {len(input_files)} files...', ascii=" ░▒▓█") as progress:
        for input_file_path in input_files:
            with open(input_file_path) as file:
                data = file.read()

            try:
                main_list = json.loads(data)
            except Exception as e:
                print(f"Error in reading file {input_file_path}: {e}")

            for airmet in main_list:
                airmet_id = airmet.get("airmet_id")
                airmet_type = airmet.get("airmet_type")
                iss_time = datetime(airmet.get("iss_year"), airmet.get("iss_month"), airmet.get("iss_day"), airmet.get("iss_hour"), airmet.get("iss_minute"))
                valid_time = datetime(airmet.get("valid_year"), airmet.get("valid_month"), airmet.get("valid_day"), airmet.get("valid_hour"), airmet.get("valid_minute"))
                airmet_raw_text = airmet.get("raw_text")

                iss_str = "ISSUED " + iss_time.strftime("%Y-%m-%d %H:%M:%S")
                valid_str = "VALID UNTIL " + valid_time.strftime("%Y-%m-%d %H:%M:%S") 

                #print(airmet_id, airmet_type, cond_str)

                subgroups = airmet.get("subgroups")
                if not subgroups:
                    subgroups = []

                _, _, _, = plot_kmz(save_dir, subgroups, airmet_type, airmet_id, airmet_raw_text, valid_time, iss_time, debug=True)
                progress.update()


    






