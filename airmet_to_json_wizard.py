import json
from os import path
from typing import Any


def enter_list(desc: str = "items") -> list[str, ...]:
    list_entry = True
    list = []

    print(f"Enter list of {desc} (or X to exit)")
    while list_entry:
        val = input("Enter item: ")
        if (val == "X") or (val == "x"):
            return list
        else:
            list.append(val)

def process_string(string: str, 
                   type: str) -> list[str, ...]:
    string = string.replace("\n", "")
    if type == "states":
        string = string.replace("AND CSTL WTRS",  "CSTLWTRS")
        list = string.split(" ")
        if "CSTLWTRS" in list:
            list[list.index("CSTLWTRS")] = "CSTL WTRS"
    elif type == "vors":
        if "-" in string:
            string = string.replace("FROM ", "")
            list = string.split("-")
        else:
            string = string.replace("FROM ", "")
            list = string.split(" TO ")
    else:
        raise ValueError("Must specify type of input ('states' or 'vors')")

    return list

def enter_subgroups() -> list[dict[Any, ...], ...]:
    subgroup_entry = True
    subgroup_count = 1
    subgroup_list = []

    print(f"Enter information for subgroup #{subgroup_count}")
    while subgroup_entry:
        subgroup = {}
        subgroup['qualifiers'] = enter_list("qualifiers")
        subgroup['vors'] = process_string(input("VORs: "), "vors")
        subgroup['states'] = process_string(input("States: "), "states")
        subgroup['desc'] = input("Description Text: ").replace("\n", " ")
        subgroup_list.append(subgroup)
        subgroup_count += 1

        val = input("Add another subgroup? (y/n)")
        if (val == "N") or (val == "n"):
            return subgroup_list

def edit_header(header: dict[Any, ...]) ->:
    print(f"Available fields: {header.keys()}")
    while True:
        field = input("Enter Field ID to edit (or X to exit): ")
        if (field == "X") or (field == "x"):
            break
        elif field not in header.keys():
            print("Please enter a valid field.")
        else:
            print(f"Current: {field} = {header[field]}")
            header[field] = input("Enter new value: ")

    header = update_time_strs(header)
    return header

def update_time_strs(header):

    header['iss_time_str'] = f"{str(header['iss_day']).zfill(2)}{str(header['iss_hour']).zfill(2)}{str(header['iss_minute']).zfill(2)}"
    header['valid_time_str'] = f"{str(header['valid_day']).zfill(2)}{str(header['valid_hour']).zfill(2)}{str(header['valid_minute']).zfill(2)}"

    return header

#enter header

airmet_list = []
airmet_header = {}

print("Add information to airmet header:")
airmet_header["iss_airport"] = input("iss_airport: ")
airmet_header["iss_year"] = int(input("iss_year: "))
airmet_header["iss_month"] = int(input("iss_month: "))
airmet_header["iss_day"] = int(input("iss_day: "))
airmet_header["iss_hour"] = int(input("iss_hour: "))
airmet_header["iss_minute"] = int(input("iss_minute: "))
airmet_header["valid_year"] = int(input("valid_year: "))
airmet_header["valid_month"] = int(input("valid_month: "))
airmet_header["valid_day"] = int(input("valid_day: "))
airmet_header["valid_hour"] = int(input("valid_hour: "))
airmet_header["valid_minute"] = int(input("valid_minute: "))

airmet_header = update_time_strs(airmet_header)

airmet_count = 1

while True:
    airmet = {}
    print(f"Enter information for AIRMET #{airmet_count}")
    
    if airmet_count > 1:
        edit_head = input("Edit header? (y/n)")
        if (edit_head == "y") or (edit_head == "Y"):
            airmet_header = edit_header(airmet_header)
            
    airmet["airmet_id"] = input("airmet_id: ")
    airmet["airmet_type"] = input("airmet_type: ")
    
    airmet['conditions'] = enter_list("conditions")
    airmet['subgroups'] = enter_subgroups()

    airmet.update(airmet_header)
    airmet_count += 1
    airmet_list += [airmet]

    val = input("Add another AIRMET? (y/n)")
    if (val == "N") or (val == "n"):
            airmet_entry = False
            break
        
save_dir = input("JSON file save directory: ")
json_obj = json.dumps(airmet_list, indent=4)
 
with open(path.join(save_dir, 'airmet.json'), "w") as outfile:
    outfile.write(json_obj)