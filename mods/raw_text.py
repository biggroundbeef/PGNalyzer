import pdb
from tqdm import tqdm 
from numpy import arange, array, zeros
from pgnalyzer.dark_funcs.raw_text import is_dark, find_best_dark
from pgnalyzer.mods.dates import update_str, update_raw

# Creating a function that will locate the start of the Pandora Instrument's data
def find_data_start(in_path: str) -> str:
    with open(in_path, encoding = "latin1") as file:
        lines = file.readlines()
        intro_line_passed = False
        header_line_passed = False
        for i, line in enumerate(lines): 
            if "-------" in line and intro_line_passed:
                header_line_passed = True
            elif "-------" in line:
                intro_line_passed = True
            elif header_line_passed and intro_line_passed:
                intro_line_passed, header_line_passed = False, False
                return line[:2]
            
# Creating a function that uses findDataStart to located the final characters of the introduction("metadata") and header            
def find_intro_header(in_path: str, data_st_chars) -> list[str, str]:
    with open(in_path, encoding = "latin1") as file:
        intro_end = file.seek(file.read().find("-\nColumn"), 0) + len("-\n")
    with open(in_path, encoding = "latin1") as file:
        header_end = file.seek(file.read().find(f"-\n{data_st_chars}"), 0) + len("-\n")
    char_dif = header_end - intro_end
    return intro_end, char_dif

def find_lat_long(PGN_intro:str) -> list[str, str]:
    """
    This function uses a PGN introduction string to find Latitude and Longitude.
    """
    for line in PGN_intro.splitlines():
        if line.startswith("Location lat"):
            file_lat = line.split()[-1]
        elif line.startswith("Location long"):
            file_long = line.split()[-1]

    return file_lat, file_long

# Creating a function that uses Header string to find each unique column name for use as a header
def find_headers(PGN_header: str) -> list:
    headers = []
    # Splitting all the lines in the txt file into a list to create header names
    for i, line in enumerate(tqdm(PGN_header.split("\n"), desc = "Creating column headers: ")):
        
        if i == len(PGN_header.split("\n")) - 4:
            for j in arange(1, int(line.split()[1].split("-")[1][:-1]) - int(line.split()[1].split("-")[0]) + 2):
                headers.append(
                    f"Pixel {j}" 
                )

        elif i == len(PGN_header.split("\n")) - 3:
            for j in arange(1, int(line.split()[1].split("-")[1][:-1]) - int(line.split()[1].split("-")[0]) + 2):
                headers.append(
                    f"Unc {j}" 
                )

        elif line.lower().startswith("column"):
            headers.append(
                " ".join(line.split()[2:7]) 
            )

    return headers

def sort_routine(letters, info, routine_dict, param_index):
    first_letter = letters[:1]
    while True:
        try: 
            if first_letter in ["A", "C", "D", "E", "F", "H", "K", "L", "P", "R", "T", "W", "Z"]:
                routine_dict[first_letter][letters].append(info)
            
            elif len(info) < 30: routine_dict["CMT"][letters].append(info)
            
            elif letters.startswith("S"):
                if is_dark(info, param_index): routine_dict["DarkS"][letters].append(info)
                else: routine_dict[first_letter][letters].append(info)
            elif letters.startswith("M"):
                if is_dark(info, param_index): routine_dict["DarkM"][letters].append(info)
                else: routine_dict[first_letter][letters].append(info)

            else:
                routine_dict["UDF"][letters].append(info)
    
            return routine_dict
        except KeyError as ke:
            # print("KEY ERROR: ", ke) 
            if first_letter in ["A", "C", "D", "E", "F", "H", "K", "L", "P", "R", "T", "W", "Z"]: routine_dict[first_letter][letters] = []
            elif len(info) < 30: routine_dict["CMT"][letters] = []
            elif is_dark(info, param_index) and first_letter == "S": routine_dict["DarkS"][letters] = []
            elif is_dark(info, param_index) and first_letter == "M": routine_dict["DarkM"][letters] = []
            elif first_letter == "S": routine_dict[first_letter][letters] = []
            elif first_letter == "M": routine_dict[first_letter][letters] = []
            else: routine_dict["UDF"][letters] = []
            
# Creating a dictionary that contains all routines and measurements
def separate_routines(
    PGN_data: list, 
    param_index: dict, 
) -> dict:
    """
    Please refer to Page 25 of Blick Software User Manual for routine reference.
    """
    routines = {
        "A" : {}, "C" : {}, "D" : {},           
        "E" : {}, "F" : {}, "H" : {},
        "K" : {}, "L" : {}, "M" : {},
        "P" : {}, "R" : {}, "S" : {},
        "T" : {}, "W" : {}, "Z" : {},
        "DarkS" : {},   "DarkM" : {}, 
        "CMT"   : {},   "UDF"   : {}
    }
     
    for line_of_raw_data in tqdm(PGN_data, desc = "Iterating Through Lines of Data: "):
        all_line_info   = line_of_raw_data.split()
        routine_letters = all_line_info[0].upper()
        routines = sort_routine(routine_letters, all_line_info, routines, param_index)
    
    for i in routines.keys():
        for j in routines[i].keys():
            try: routines[i][j] = array(routines[i][j])
            except: pass

    return routines

 # Creating a function that opens a PGN L0 file and returns a python dictionary of file contents

# This function is a combination of all previous functions
def read_PGN(input_file_path: str) -> list[dict, dict]:
    """
    Using a PGN L0 as an input path, this function sorts through all the data
    abiding by the Blick Software User Manual. Outputs are a dictionary of 
    all the different routines and a dictionary of the important Pandora Spectrometer
    information (ex: Number/Location/etc).
    """
    
    data_start_chars = find_data_start(input_file_path)
    intro_end_char, header_end_char = find_intro_header(input_file_path, data_start_chars)
    
    with open(input_file_path, encoding = "latin1") as file:
        
        raw_info_intro_str   = file.read(intro_end_char)
        raw_info_header_str  = file.read(header_end_char)
        raw_info_data_list   = file.readlines()

        latitude, longitude  = find_lat_long(raw_info_intro_str)
        headers              = find_headers(raw_info_header_str)
        param_idx            = dict(zip(headers, arange(0, len(headers))))

        sorted_routines_out  = separate_routines(raw_info_data_list, param_idx)
        dora_data_out        = {"LATI" : latitude, "LONG" : longitude, "P_IDX" : param_idx}

        return sorted_routines_out, dora_data_out
    
# A function that brings all previous functions together to read the PGN file
def apply_dark_correction(in_data, in_dark_data, in_param_idx):

    bright_counts     = in_data[:, in_param_idx["Pixel 1"]:in_param_idx["Unc 1"]].astype(float)                    
    bright_Fs         = in_data[:, in_param_idx["Scale factor for data (to"]].astype(float)

    timestamps = array(list(map(
        lambda dataentry: update_str(update_raw(dataentry[in_param_idx["UT date and time for"]])), 
        in_data
    )))
    exp_times = array(
        list(map(
            lambda dataentry: float(dataentry[in_param_idx["Integration time [ms]"]]), 
            in_data
        )),
        dtype = float
    )

    ccs = zeros(len(bright_counts), dtype = object)
    used_darks = zeros(len(bright_counts), dtype = object)
    for i, timestamp in enumerate(tqdm(timestamps, desc = "Applying Dark Correction")):
        
        best_dark = find_best_dark(timestamp, exp_times[i], in_dark_data, in_param_idx)
        best_dark_pix = best_dark[in_param_idx["Pixel 1"]:in_param_idx["Unc 1"]].astype(float)
        best_dark_F = best_dark[in_param_idx["Scale factor for data (to"]].astype(float)
        
        ccs[i] = (
            ( (bright_counts[i] / bright_Fs[i]) - (best_dark_pix / best_dark_F) )
            / exp_times[i]
        )

        used_darks[i] = best_dark

    return {
        "CC"        : ccs,
        "DATES"     : timestamps,
        "BRIGHTS"   : in_data,
        "DARKS"     : used_darks,
        "P_IDX"     : in_param_idx
    }