import matplotlib.pyplot as plt
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
    
    print("RA1 = ", ra1)
    print("DEC1 = ", dec1)
    print("RA2 = ", ra2)
    print("DEC2 = ", dec2)
    print("")

    ra1, dec1, ra2, dec2 = map(np.radians, [ra1, dec1, ra2, dec2])
    
    ra1 = np.round(ra1, 6)
    dec1 = np.round(dec1, 6)
    ra2 = np.round(ra2, 6)
    dec2 = np.round(dec2, 6)


    print("RA1 = ", ra1)
    print("DEC1 = ", dec1)
    print("RA2 = ", ra2)
    print("DEC2 = ", dec2)
    print("")

    # haversine formula 
    dra = np.round(ra2 - ra1, 6)
    ddec = np.round(dec2 - dec1, 6)


    print("DRA = ", dra)
    print("DDEC = ", ddec)
    
    print("")

    s1 = np.round(np.sin(ddec/2)**2, 10)

    print("S_DDEC = ", s1)

    s2 = np.round(np.sin(dra/2)**2, 10)

    print("S_DRA = ", s2)
    print("")

    a = np.round(s1 + np.cos(dec1) * np.cos(dec2) * s2, 10)
    c = 2 * np.arcsin(np.sqrt(a))
    return c

def proper_motion(ra1, dec1, ra2, dec2, dt):

    print("RA1 = ", ra1)
    print("DEC1 = ", dec1)
    print("RA2 = ", ra2)
    print("DEC2 = ", dec2)
    print("")
    
    ra1 = np.round(Angle(ra1, unit = u.deg).hour * 15, 5)
    dec1 = np.round(Angle(dec1).degree, 5)
    ra2 = np.round(Angle(ra2, unit = u.degree).hour * 15, 5)
    dec2 = np.round(Angle(dec2).degree, 5)

    AS = np.round(haversine(ra1, dec1, ra2, dec2), 7)
    print("AS (Rad) = ", AS)
    AS = np.degrees(AS)
    print("AS (Deg) = ", AS)
    AS = AS * 3600
    print("AS (arcseconds) = ", AS)

    pm = AS / dt

    print("Proper motion (arcsecond / yr) = ", pm)

print("Barnard's Star")
proper_motion("17h57m47.13s", "4d41m34.1s", "17h57m46.64s", "4d43m17.1s", 10) # Bernard's Star
print("\nProxima Centauri")
proper_motion("14h29m46.17", "-62d40m36.9s", "14h29m40.96s", "-62d40m22.9", 10) # Proxima Centauri
print("\nVega")
proper_motion("18h36m54.56s", "38d47m0.1s", "18h36m54.64s", "38d47m3.4s", 10) # Vega
print("\nSirius")
proper_motion("6h45m10.36s", "-16d42m59.2s", "6h45m9.93s", "-16d43m12.1s", 10) # Sirius
