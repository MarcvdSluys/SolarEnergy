# SolarEnergy #

![PyPI](https://img.shields.io/pypi/v/solarenergy?color=%230A0)
![PyPI - Downloads](https://img.shields.io/pypi/dm/solarenergy)
[![Code check](https://github.com/MarcvdSluys/SolarEnergy/actions/workflows/code-check.yml/badge.svg)](https://github.com/MarcvdSluys/SolarEnergy/actions/workflows/code-check.yml)
[![Documentation
Status](https://readthedocs.org/projects/solarenergy/badge/?version=latest)](https://solarenergy.readthedocs.io/en/latest/?badge=latest)
![PyPI - Licence](https://img.shields.io/pypi/l/solarenergy?color=%230A0)

A Python module to do simple modelling in the field of solar energy.  The code is being developed by [Marc van
der Sluys](http://marc.vandersluys.nl) of the department of astrophysics of the Radboud University Nijmegen,
the Netherlands and the department of Sustainable energy of the HAN University of Applied Sciences in Arnhem,
the Netherlands.  The SolarEnergy package can be used under the conditions of the EUPL 1.2 licence.


## Installation ##

This package can be installed using `pip install solarenergy`.  This should automatically install the
dependency packages `pytz`, `numpy` `pandas` and `soltrack` if they haven't been installed already.  If you
are installing by hand, ensure that these packages are installed as well.


## Modes and speed ##

SolarEnergy start with the computation of the position of the Sun for a given instant (scalar) or for a series
of instances using an array or vector.  The latter is faster than the former (for ~2 instances or more) and
may be easier to use, depending on the problem given.

The table below shows the speed with which the position of the Sun is computed (in number of positions per
second) as a function of the size of the dataset.  The runs were timed on a single hyperthreaded core capped
at 3.4GHz and the minimum time of 10 runs is displayed.  The scalar runs scale linearly, the speed peaks
around 1e5-1e6 elements for vectors.  Timings below are for timezone-naive datetimes (and hence UTC), if
timezone-aware datetimes are used, the calculations take about 4.4 times longer(!)

| Mode   | Number of instances | Time (s) | Speed (/s) |
|--------|---------------------|----------|------------|
| scalar | 1e3                 | 0.680    | 1471       |
|        |                     |          |            |
| vector | 10                  | 8.86e-4  | 11,287     |
| vector | 100                 | 1.15e-3  | 86,957     |
| vector | 1000                | 2.19e-3  | 456,621    |
| vector | 1e4                 | 1.29e-2  | 775,194    |
| vector | 1e5                 | 6.85e-2  | 1,459,854  |
| vector | 1e6                 | 0.656    | 1,524,390  |
| vector | 1e7                 | 8.56     | 1,168,224  |



## Example for a single calculation ##

```python
import solarenergy as se
import numpy as np

# Location of my solar panels:
geoLon =  5*se.d2r  # Geographic longitude (>0 for eastern hemisphere; ° -> rad)
geoLat = 52*se.d2r  # Geographic latitude  (>0 for northern hemisphere; ° -> rad)

# Orientation of my solar panels:
spAz   = -2*se.d2r  # Azimuth ('wind direction') of my panels are facing.  Note: South=0, W=90° (pi/2 rad) in the northern hemisphere!  (rad)
spIncl = 28*se.d2r  # Inclination of my panels w.r.t. the horizontal  (rad)

# An hour past noon local time on 1 March 2020:
myTZ  = 'Europe/Amsterdam'
year  = 2020
month = 3
day   = 1
hour  = 13

# Compute Sun position (uses SolTrack behind the scenes):
sunAz,sunAlt,sunDist = se.sun_position_from_date_and_time(geoLon,geoLat, year,month,day, hour, timezone=myTZ)

AM        = se.airmass(sunAlt)                                  # Air mass for this Sun altitude
extFac    = se.extinction_factor(AM)                            # Extinction factor at sea level for this airmass
cosTheta  = se.cos_angle_sun_panels(spAz,spIncl, sunAz,sunAlt)  # cos of the angle with which Sun hits my panels
theta     = np.arccos(cosTheta)                                 # Angle with which Sun hits my panels

Iext      = se.sol_const / sunDist**2                           # Extraterrestrial radiation = Solar constant, scaled with distance
DNIcs     = Iext / extFac                                       # DNI for a clear sky
dirRad    = DNIcs * cosTheta                                    # Insolation of direct sunlight on my panels


# Print input and output:
print("Location:                    %0.3lf E, %0.3lf N"  % (geoLon*se.r2d, geoLat*se.r2d))
print("Date:                        %4d-%2.2d-%2.2d"     % (year, month, day))
print("Time:                        %2d:00"              % (hour))
print()

print("Sun azimuth:                 %7.3lf°"             % (sunAz*se.r2d))
print("Sun altitude:                %7.3lf°"             % (sunAlt*se.r2d))
print("Sun distance:                %7.4lf AU"           % (sunDist))
print()

print("Air mass:                    %7.3lf"              % (AM))
print("Extinction factor:           %7.3lf"              % (extFac))
print("Sun-panels angle:            %7.1lf°"             % (theta*se.r2d))
print()

print("Solar constant:              %7.1lf W/m²"         % (se.sol_const))
print("Extraterrestrial radiation:  %7.1lf W/m²"         % (Iext))
print("DNI (clear sky):             %7.1lf W/m²"         % (DNIcs))
print("Direct insolation:           %7.1lf W/m²"         % (dirRad))
print()
```

## Example for a range of instances ##

```python
#!/bin/env python

"""Example Python script using the SolarEnergy module for a range of instances."""

import solarenergy as se
from solarenergy import d2r,r2d

# Location of my solar panels:
geoLon =  5*d2r  # Geographic longitude (>0 for eastern hemisphere; ° -> rad)
geoLat = 52*d2r  # Geographic latitude  (>0 for northern hemisphere; ° -> rad)

# Orientation of my solar panels:
spAz   = -2*d2r  # Azimuth ('wind direction') of my panels are facing.  Note: South=0, W=90° (pi/2 rad) in the northern hemisphere!  (rad)
spIncl = 28*d2r  # Inclination of my panels w.r.t. the horizontal  (rad)


import pandas as pd
dates = pd.date_range(pd.to_datetime('2022-03-21'), pd.to_datetime('2022-03-22'), freq='1h')  # DatetimeIndex 0-24h
dates = dates.tz_localize('Europe/Amsterdam')
df = pd.DataFrame(data=dates, columns=['dtm'])  # Create a Pandas DataFrame with the datetimes as first column

# Compute Sun position (uses SolTrack behind the scenes) and add it as three columns to the df:
df['sunAz'],df['sunAlt'],df['sunDist'] = se.sun_position_from_datetime(geoLon,geoLat, df['dtm'])

df['I_ext']     = 1361.5 / df.sunDist**2                                 # Extraterrestrial radiation (at the top of the atmosphere; AM0)

df['AM']        = se.airmass(df.sunAlt)                                  # Air mass for this Sun altitude
df['extFac']    = se.extinction_factor(df.AM)                            # Extinction factor at sea level for this airmass
df['DNI']       = df.I_ext / df.extFac                                   # DNI for a clear sky

df['cosTheta']  = se.cos_angle_sun_panels(spAz,spIncl, df.sunAz,df.sunAlt)  # cos of the angle with which Sun hits my panels
df['dirRad']    = df.DNI * df.cosTheta                                   # Insolation of direct sunlight on my panels

df.sunAz  *= r2d  # Convert azimuth and ...
df.sunAlt *= r2d  # ... altitude to degrees

print(df[df.sunAlt > 0])  # Print the results for the hours where the Sun is above the horizon
```


## SolarEnergy pages ##

* [Pypi](https://pypi.org/project/solarenergy/): SolarEnergy Python package
* [GitHub](https://github.com/MarcvdSluys/SolarEnergy): SolarEnergy source code
* [ReadTheDocs](https://solarenergy.readthedocs.io): SolarEnergy documentation


## Author and licence ##

* Author: Marc van der Sluys
* Contact: http://marc.vandersluys.nl
* Licence: [EUPL 1.2](https://www.eupl.eu/1.2/en/)


## References ##

* This Python code is adapted from the Fortran implementation in
  [libTheSky](http://libthesky.sourceforge.net/), which contains many references.
* [Celestial mechanics in a nutshell (CMiaNS)](https://cmians.sourceforge.io/)
