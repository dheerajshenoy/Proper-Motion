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

table = pd.DataFrame(columns = [ "Star Name", "Years", "RA (degree)", "Dec (degree)", "Angular Separation (rad)", "Proper Motions (arcsec / yr)", "Proper Motion AVG (arcsec / yr)" ])

# Haversine formula to find the angular separation between two points
def haversine(lon1, lat1, lon2, lat2):
    """
    Calculate the great circle distance in kilometers between two points 
    on the earth (specified in decimal degrees)
    """
    # convert decimal degrees to radians 
    lon1, lat1, lon2, lat2 = map(np.radians, [lon1, lat1, lon2, lat2])

    # haversine formula 
    dlon = lon2 - lon1 
    dlat = lat2 - lat1 
    a = np.sin(dlat/2)**2 + np.cos(lat1) * np.cos(lat2) * np.sin(dlon/2)**2
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
    as_list = []
    pm_list = []

    for y in ys:
        t = ts.tt(y)
        star = Star.from_dataframe(df.loc[HIP])

        astro = earth.at(t).observe(star)

        ra, dec, _ = astro.radec()

        ra_deg = ra.hours * 15
        dec_deg = dec.degrees

        ra_list.append(ra_deg[0])
        dec_list.append(dec_deg[0])

    
    # ra1, ra2 = ra_list[0], ra_list[1]
    # dec1, dec2 = dec_list[0], dec_list[1]
    #
    # AS = Angle(haversine(ra1, dec1, ra2, dec2), unit = u.rad)
    # pm = AS.to(u.arcsec)[0] / (years * u.year)

    # Calculate average of angular separation
    pm = 0

    for i in range(len(ra_list) - 1):
        ra1, ra2 = ra_list[i], ra_list[i + 1]
        dec1, dec2 = dec_list[i], dec_list[i + 1]

        AS = Angle(haversine(ra1, dec1, ra2, dec2), unit = u.rad)
        as_list.append(AS.to(u.deg).value)
        pm = pm + AS.to(u.arcsec) / (years * u.year)
        pm_list.append(pm)
        
    ys = [str(y) for y in ys]
    # table = pd.DataFrame(columns = [ "Star Name", "Years", "RA", "Dec", "Angular Separation", "Proper Motion" ])
    # table.loc[len(table)] = [star_name, ys, sum(ra_list)/len(ra_list), sum(dec_list)/len(dec_list), sum(as_list)/len(as_list), pm.value] 
    table.loc[len(table)] = [star_name, ys, ra_list, dec_list, as_list, pm_list, pm.value] 
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

# pm = proper_motion(2000, 2010, 5, [ "Proxima Centauri", "Barnard's Star", "61 Cygni A"])
pm = proper_motion(2000, 2010, 2, "Proxima Centauri")
print(table.to_markdown())

