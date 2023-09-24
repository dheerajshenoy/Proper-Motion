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

table = pd.DataFrame(columns = [ "Star Name", "Year", "RA (hms)", "RA (degree)", "Angular Separation (deg)", "Proper Motions (arcsec / yr)", "Proper Motion AVG (arcsec / yr)" ])

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
    ra_deg_list = []
    as_list = []
    pm_list = []

    for y in ys:
        # years, month, day, hour, minute, seconds
        t = ts.tt(y, 1, 1, 0, 0, 0)
        star = Star.from_dataframe(df.loc[HIP])
        astro = earth.at(t).observe(star)
        ra, _, _ = astro.radec()
        ra_list.append(ra)
        ra_deg = ra.hours * 15
        ra_deg_list.append(ra_deg[0])

    pm = 0

    for i in range(len(ra_deg_list) - 1):
        ra1, ra2 = ra_deg_list[i], ra_deg_list[i + 1]

        AS = Angle(ra2 - ra1, unit = u.deg)
        as_list.append(AS.to(u.deg).value)
        pm = pm + AS.to(u.arcsec) / (years * u.year)
        pm_list.append(pm)
    
    PM = sum(pm_list)/len(pm_list)
    ys = [str(y) for y in ys]
    # table = pd.DataFrame(columns = [ "Star Name", "Years", "RA", "Dec", "Angular Separation", "Proper Motion" ])
    # table.loc[len(table)] = [star_name, ys, sum(ra_list)/len(ra_list), sum(dec_list)/len(dec_list), sum(as_list)/len(as_list), pm.value] 
    table.loc[len(table)] = [star_name, ys, ra_list, ra_deg_list, as_list, pm_list, PM] 
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
pm = proper_motion(2000, 2010, 2, ["Proxima Centauri", "Vega", "Fang"])
print(table.to_markdown())
