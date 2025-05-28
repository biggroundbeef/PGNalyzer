
from numpy import pi, deg2rad, sin, cos, rad2deg, arccos
from astropy.time import Time
from astropy.coordinates import get_sun
from pgnalyzer.mods.dates import ct2lst

def get_solar_zenith_angle(date, lati, long, inDegrees:bool = True):
    hour2radian = pi / 12 # there are 12 hours in pi radians 
    jd = Time(date, format = "datetime64").jd
    time = Time(jd, format = "jd")
    lst = ct2lst(jd, long) # This uses an AI generated conversion of IDL's ct2lst -> python

    ra_zen = lst * hour2radian
    dec_zen = deg2rad(lati)

    sun_pos = get_sun(time)
    ra_sun = sun_pos.ra.to_value(unit = "rad")
    dec_sun = sun_pos.dec.to_value(unit = "rad")

    th1 = pi/2 - dec_zen
    th2 = pi/2 - dec_sun
    ph1 = ra_zen
    ph2 = ra_sun

    sth1 = sin(th1)
    cth1 = cos(th1)
    sph1 = sin(ph1)
    cph1 = cos(ph1)

    sth2 = sin(th2)
    cth2 = cos(th2)
    sph2 = sin(ph2)
    cph2 = cos(ph2)

    x1 = sth1*cph1
    y1 = sth1*sph1
    z1 = cth1

    x2 = sth2*cph2
    y2 = sth2*sph2
    z2 = cth2

    sza = arccos(x1*x2 + y1*y2 + z1*z2)
    if inDegrees:
        return rad2deg(sza)
    else: return sza

