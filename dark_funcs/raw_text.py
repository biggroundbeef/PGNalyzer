import pdb
from numpy import argmin, array, datetime64, timedelta64, where, zeros
from pgnalyzer.mods.dates import update_str, update_raw

# Creating a function that uses the dataentry's Filterwheels as a check for dark measurements
def is_dark(dataentry:str, param_index) -> bool:
    return dataentry[param_index["Position of filterwheel #2, 0=filterwheel"]] == "3"

# Creating a function that takes full dark and separates it to calculation-critical information
def simplify_darks(self, dark_routine:list, param_index:dict) -> dict:
    simple_darks = {}
    for dark in dark_routine:
        date = dark[param_index["UT date and time for"]]
        vals = dark[param_index["Pixel 1"]:param_index["Unc 1"]]
        exp_time = dark[param_index["Integration time [ms]"]]
        scale = dark[param_index["Scale factor for data (to"]]
        simple_darks[date] = {
            "PIX"       : array(list(map(lambda val: float(val), vals))),
            "SCALE"     : float(scale),
            "EXP_TIME"  : float(exp_time)  
        }
    return simple_darks

# Creating a function that selected the dark closest to a sun measurment in time and exposure_time        
def find_best_dark(target_date: datetime64, exp_time: float, dark_entries: dict, param_index: dict) -> dict:
    max_date = (target_date + timedelta64(5, "m")).astype(datetime64)
    min_date = (target_date - timedelta64(5, "m")).astype(datetime64)

    all_dark_dates = list(map(lambda dark_entry: update_str(update_raw(dark_entry[param_index["UT date and time for"]])), dark_entries))

    dark_date_search = [
        date for date in all_dark_dates
        if (date - min_date).astype(float) > 0 and 
        (date - max_date).astype(float) < 0
    ]

    if len(dark_date_search) > 1:
        diffs = zeros(len(dark_date_search))
        for i, timestamp in enumerate(dark_date_search): 
            diffs[i] = timestamp - timedelta64(int(exp_time * 1e6), "ns")
        closest_date = dark_date_search[argmin(diffs)]
    else: closest_date = dark_date_search[0]

    best_dark_idx = where(all_dark_dates == closest_date)[0]

    return dark_entries[best_dark_idx].flatten()
