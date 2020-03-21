#  Copyright (c) 2020  Marc van der Sluys - marc.vandersluys.nl
#  
#  This file is part of the SolarEnergy Python package, containing a Python module do simple modelling in the field of
#  solar energy.  See: https://github.com/MarcvdSluys/SolarEnergy
#  
#  This is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License
#  as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
#  
#  This software is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied
#  warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.
#  
#  You should have received a copy of the GNU General Public License along with this code.  If not, see
#  <http://www.gnu.org/licenses/>.


"""Functions for solar energy dealing with (solar) radiation.

References:
  - M. van der Sluys, Celestial mechanics in a nutshell, https://cmians.sourceforge.io (2014).
      
"""


import sys
import datetime as dt
import pytz as tz
import numpy as np

import soltrack as st


def computeSunPos(lon,lat, year,month,day, hour,minute=0,second=0, timezone='UTC', debug=False):
    """Compute the Sun local position (azimuth, altitude and distance) for the given geographical location and
    date and time using SolTrack.
    
    Args:
        lon (float):  Geographic longitude to compute the Sun position for (rad).
        lat (float):  Geographic latitude to compute the Sun position for (rad).
    
        year (int):    Year (CE) to compute the Sun position for.
        month (int):   Month to compute the Sun position for.
        day (int):     Day of month to compute the Sun position for.
    
        hour (int):    Hour of day to compute the Sun position for (local time!).
        minute (int):  Minute to compute the Sun position for (optional; default = 0).
        second (int):  Second to compute the Sun position for (optional; default = 0).
    
        timezone (timezone):  Time zone for which date and time are provided (optional; default = 'UTC').
    
        debug (bool):  Switch to write detailed output to stdout (optional; default = False).
    
    Returns:
        tuple (float,float,float):  Tuple containing (azimuth, altitude, distance):
    
            - azimuth (float):   Azimuth of the Sun (rad; south = 0 on the northern hemisphere).
            - altitude (float):  Altitude of the Sun (rad).
            - distance (float):  Distance Sun-Earth (AU).
    
    """
    
    myTZ = tz.timezone(timezone)   # My timezone
    myTime = dt.datetime(int(year),int(month),int(day), int(hour),int(minute),int(second))  # Time w/o timezone
    lt = myTZ.localize(myTime, is_dst=None)  # Mark as local time
    utc = lt.astimezone(tz.utc)              # Convert to UTC
    
    # Set up geographical location (in degrees, since useDegrees=True) in a SolTrack Location dataclass object:
    loc = st.Location(lon,lat)  # longitude (>0: east of Greenwich),  latitude (>0: northern hemisphere), in radians
    
    # Set (UT!) date and time in a SolTrack Time dataclass object:
    time = st.Time.datetime2st(utc)
    
    # Compute positions - returns a st.Position object:
    pos = st.computeSunPosition(loc, time, computeDistance=True)
    
    # Write results to standard output:
    if(debug):
        r2d = 180/np.pi  # Convert radians to degrees
        print("Location:  %0.3lf E, %0.3lf N"  % (loc.longitude*r2d, loc.latitude*r2d))
        print("Date:      %4d %2d %2d"         % (time.year, time.month, time.day))
        print("Time:      %2d %2d %9.6lf"      % (time.hour, time.minute, time.second))
        print("JD:        %0.11lf"             % (pos.julianDay))
        print()
        
        print("Corrected azimuth, altitude:  %10.6lf° %10.6lf°" % (pos.azimuthRefract*r2d, pos.altitudeRefract*r2d))
        print("Distance:                     %10.6lf AU"        % (pos.distance))
        print()
    
    return pos.azimuthRefract, pos.altitudeRefract, pos.distance



def angleSunPanels(spAz,spIncl, sunAz,sunAlt):
    """Compute the cosine of the angle between the orientation vector of the solar panels and the position vector
       of the Sun.
    
    This is the cosine of the angle under which the direct sunlight hits the solar panels.  Multiply it with
    the DNI to obtain the direct insolation on the panels.
    
    Args:
        spAz (float):    Azimuth in which the solar panels are facing (rad; south = 0 on the northern hemisphere).
        spIncl (float):  Inclination ('zenith angle') of the solar panels w.r.t. the horizontal (rad).
        sunAz (float):   Azimuth of the Sun (rad; south = 0 on the northern hemisphere).
        sunAlt (float):  Altitude of the Sun (rad).
    
    Returns:
        float:  The cosine between the normal vector of the solar panels and the position vector of the Sun (rad).

    """
    
    cosTheta = np.sin(sunAlt) * np.cos(spIncl)  +  np.cos(sunAlt) * np.sin(spIncl) * np.cos(sunAz - spAz)
    
    return cosTheta



def airmass(sunAlt):
    """Compute airmass as a function of Sun altitude.
    
    Args:
        sunAlt (float):  Altitude of the Sun (rad).
    
    Returns:
        float:  Airmass at sea level (AM~1 if the Sun is in the zenith, AM~38 near the horizon).
    
    """
    
    if(sunAlt < -0.00989):
        airmass = 1000 * (0.15 + abs(sunAlt))  # Very bad, but still getting worse for lower Sun, for solvers
    else:
        airmass = (1.002432*np.sin(sunAlt)**2 + 0.148386*np.sin(sunAlt) + 0.0096467) / \
                  (np.sin(sunAlt)**2*np.sin(sunAlt) + 0.149864*np.sin(sunAlt)**2 + 0.0102963*np.sin(sunAlt) + 0.000303978)
        airmass = max( airmass, 1 )   # Air mass cannot be lower than 1
    
    return airmass



def extintionFactor(airmass):
    """Compute the atmospheric extinction factor for sunlight from the air mass.
    
    Args:
        airmass (float):  Airmass at sea level (AM~1 if the Sun is in the zenith, AM~38 near the horizon).
    
    Returns:
        float:  The extinciton factor for sunlight in the atmosphere.  Divide the solar constant by this number 
                to obtain the DNI.
    
    """
    
    if(airmass > 38.2):
        extFac = np.sqrt(sys.float_info.max) * (0.15 + airmass)  # Very bad, but still getting worse for higher airmass, for solvers
    else:
        coefs = [ 9.1619283e-2, 2.6098406e-1,-3.6487512e-2, 6.4036283e-3,-8.1993861e-4, 6.9994043e-5,-3.8980993e-6,
                  1.3929599e-7, -3.0685834e-9, 3.7844273e-11,-1.9955057e-13]
        
        AMpow = 1.0                      # AM^0
        ext = coefs[0]                   # c_1 * AM^0
        for iCoef in range(1,len(coefs)):
            AMpow *= airmass             # AM^(i-1)
            ext += coefs[iCoef] * AMpow  # + c_i * AM^(i-1)
            
        extFac = np.exp(ext)
    
    return extFac


