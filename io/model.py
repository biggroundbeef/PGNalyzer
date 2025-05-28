
from pandas import read_table
from numpy import array

def get_model(input_file_path:str):
    data = read_table(input_file_path, names = ["wavelength", "irradiance"], delimiter = " ", encoding = "latin1")
    return array(list(map(lambda val: float(val), data["wavelength"][data["irradiance"].values != 0]))), array(list(map(lambda val: float(val), data["irradiance"][data["irradiance"].values != 0])))
