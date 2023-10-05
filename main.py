import matplotlib.pyplot as plt
from skyfield.api import Star, load
from skyfield.data import hipparcos
from astroquery.vizier import Vizier
from astropy import units as u
import numpy as np
from astropy.coordinates import SkyCoord, Angle, angular_separation
from astropy.units.format.utils import warnings
import astroquery as aq
import pandas as pd

warnings.simplefilter('ignore', category = u.UnitsWarning)

table = pd.DataFrame(columns = [ "Star Name", "Year", "RA (hms)", "RA (degree)", "Dec (hms)", "Dec (degree)", "Angular Separation (degree)", "Proper Motions (arcsec / yr)", "Proper Motion AVG (arcsec / yr)" ])

# Haversine formula to find the angular separation between two points
def haversine(ra1, dec1, ra2, dec2):
    """
    Calculate the great circle distance in kilometers between two points 
    on the earth (specified in decimal degrees)
    """
    # convert decimal degrees to radians 

    ra1, dec1, ra2, dec2 = map(np.radians, [ra1, dec1, ra2, dec2])

    # haversine formula 
    dra = ra2 - ra1 
    ddec = dec2 - dec1 

    a = np.sin(ddec/2)**2 + np.cos(dec1) * np.cos(dec2) * np.sin(dra/2)**2
    c = 2 * np.arcsin(np.sqrt(a))
    return c

ts = load.timescale()

# Load the planets data
planets = load('de421.bsp')
earth = planets["earth"]

# Helper function to find the proper motion
def dd(y0, yf, n, star_name):
    star = Vizier.query_object(star_name)
    HIP = star["I/239/hip_main"]["HIP"]
    ys = np.linspace(y0, yf, n)
    years = yf - y0

    with load.open(hipparcos.URL) as f:
        df = hipparcos.load_dataframe(f)

    
    ra_list = []
    dec_list = []
    ra_deg_list = []
    dec_deg_list = []
    as_list = []
    pm_list = []

    for y in ys:
        # years, month, day, hour, minute, seconds
        t = ts.tt(y, 1, 1, 0, 0, 0)
        star = Star.from_dataframe(df.loc[HIP])
        astro = earth.at(t).observe(star)
        ra, dec, _ = astro.radec()
        ra_list.append(ra)
        dec_list.append(dec)
        ra_deg = ra.hours * 15
        dec_deg = dec.degrees
        ra_deg_list.append(ra_deg[0])
        dec_deg_list.append(dec_deg[0])

    pm = 0

    for i in range(len(ra_deg_list) - 1):
        ra1, ra2 = ra_deg_list[i], ra_deg_list[i + 1]
        dec1, dec2 = dec_deg_list[i], dec_deg_list[i + 1]

        AS = Angle(haversine(ra1, dec1, ra2, dec2), unit = u.rad)
        as_list.append(AS.to(u.deg).value)
        pm = pm + AS.to(u.arcsec) / (years * u.year)
        pm_list.append(pm)
 
    ys = [str(y) for y in ys]
    # table = pd.DataFrame(columns = [ "Star Name", "Years", "RA", "Dec", "Angular Separation", "Proper Motion" ])
    # table.loc[len(table)] = [star_name, ys, sum(ra_list)/len(ra_list), sum(dec_list)/len(dec_list), sum(as_list)/len(as_list), pm.value] 
    table.loc[len(table)] = [star_name, ys, ra_list, ra_deg_list, dec_list, dec_deg_list, as_list, pm_list, pm.value] 
    return pm.value

# Function that returns the proper motion of star(s) in arcsec / year
def proper_motion(y0, yf, n, star_name):
    
    if type(star_name) == list:
        pm_list = []
        for name in star_name:
            pm_list.append(dd(y0, yf, n, name))
        return pm_list
    else:
        return dd(y0, yf, n, star_name)

pm = proper_motion(2000, 2001, 2, [ "Proxima Centauri", "Barnard's Star", "Vega"])
print(table.to_markdown())
