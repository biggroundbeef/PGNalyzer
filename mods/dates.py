from numpy import datetime64

def update_raw(input_date_raw:str) -> str:
    """
    Converts a date input from L0 data to a numpy datetime64 compatible string.

    Keyword Arguments:  
    input_date_raw 
        str
        This is a string containing a date from L0 data.  
        Ex: "20250525T143528.8Z"

    returns str
        Year    ->  input_date_raw[:4]
        Month   ->  input_date_raw[4:6]
        Day     ->  input_date_raw[6:8]
        Hour    ->  input_date_raw[9:11]
        Minute  ->  input_date_raw[11:13]
        Second  ->  input_date_raw[13:-1]
        Ex:"2025-05-25T14:35:28.8"
    """
    return f"{input_date_raw[:4]}-{input_date_raw[4:6]}-{input_date_raw[6:8]}T{input_date_raw[9:11]}:{input_date_raw[11:13]}:{input_date_raw[13:-1]}"

def update_str(input_date_str:str) -> datetime64:
    """
    Converts a date input from L0 data to a numpy datetime64 object
    
    Keyword Arguments:
    input_date_str 
        str
        This is a numpy datetime64 compatible string.
        Ex: "2025-05-25T14:35:28.8"
    
    returns numpy.datetime64
    """
    return datetime64(input_date_str)

# ChatGPT IDL function -> python    
def ct2lst(jd, lon):
    """
    Convert Julian Date and longitude to Local Sidereal Time (LST).
    LST is needed to compute the observer's right ascension.
    """
    T = (jd - 2451545.0) / 36525.0  # Julian centuries since J2000.0
    GMST = 280.46061837 + 360.98564736629 * (jd - 2451545.0) + T**2 * (0.000387933 - T / 38710000.0)
    GMST = GMST % 360.0  # Normalize to [0, 360)
    LST = (GMST + lon) % 360.0  # Add longitude to get local sidereal time
    return LST / 15.0  # Convert degrees to hours

def quick_test():
    reference_date = "20250505T050448.8Z"

    fixed_str = update_raw(reference_date)
    numpy_datetime64 = update_str(fixed_str)

    print("Success!")

if __name__ == "__main__":
    quick_test()