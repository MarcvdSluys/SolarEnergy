# SolarEnergy

![PyPI](https://img.shields.io/pypi/v/solarenergy?color=%230A0)
![PyPI - Downloads](https://img.shields.io/pypi/dm/solarenergy)
[![Code check](https://github.com/MarcvdSluys/SolarEnergy/actions/workflows/code-check.yml/badge.svg)](https://github.com/MarcvdSluys/SolarEnergy/actions/workflows/code-check.yml)
[![Documentation
Status](https://readthedocs.org/projects/solarenergy/badge/?version=latest)](https://solarenergy.readthedocs.io/en/latest/?badge=latest)
![PyPI - Licence](https://img.shields.io/pypi/l/solarenergy?color=%230A0)

A Python module to do simple modelling in the field of solar energy.  The code is being developed by [Marc van
der Sluys](http://marc.vandersluys.nl) of the department of astrophysics of the Radboud University Nijmegen,
the Netherlands and the department of Sustainable energy of the HAN University of Applied Sciences in Arnhem,
the Netherlands, now at the Netherlands Institute for Nuclear and High-Energy Physics (Nikhef) and the
Institute for Gravitational and Subatomic Physics (GRASP) at Utrecht University in the Netherlands. The
SolarEnergy package can be used under the conditions of the EUPL 1.2 licence.


## Installation

This package can be installed using `pip install solarenergy`.  This should automatically install the
dependency packages `numpy`, `pytz`, `astroconst` and `soltrack` (>=0.2.0) if they haven't been installed
already.  If you are installing by hand, ensure that these packages are installed as well (if you're not using
a Python version older than 3.7, you will need to install `dataclasses` in addition).



## Code examples

### Code example for a numer or range (array, vector) of instances

In this mode, we prepare a list of datetimes (Pandas Series, DatetimeIndex, ndarrays of datetime64, ...) and
feed that to SolarEnergy at once, for better performance and easier use.

```python
"""Example Python script using the SolarEnergy module for a range of instances."""

# Location of my solar panels:
from solarenergy import d2r,r2d
geoLon =  5*d2r  # Geographic longitude (>0 for eastern hemisphere; ° -> rad)
geoLat = 52*d2r  # Geographic latitude  (>0 for northern hemisphere; ° -> rad)

# Orientation of my solar panels:
spAz   = -2*d2r  # Azimuth ('wind direction') of my panels are facing.  Note: South=0, W=90° (pi/2 rad) in the northern hemisphere!  (rad)
spIncl = 28*d2r  # Inclination of my panels w.r.t. the horizontal  (rad)


import pandas as pd
dates = pd.date_range('2022-03-21', pd.to_datetime('2022-03-22'), freq='1h', tz='Europe/Amsterdam')  # DatetimeIndex 0-24h
df = pd.DataFrame(data=dates, columns=['dtm'])  # Create a Pandas DataFrame with the datetimes as first column

# Compute Sun positions (using SolTrack behind the scenes) and add them as three columns to the df:
import solarenergy as se
df['sunAz'],df['sunAlt'],df['sunDist'] = se.sun_position_from_datetime(geoLon,geoLat, df['dtm'])

df['I_ext']     = 1361.5 / df.sunDist**2                                 # Extraterrestrial radiation (at the top of the atmosphere; AM0)

df['AM']        = se.airmass(df.sunAlt)                                  # Air mass for this Sun altitude
df['extFac']    = se.extinction_factor(df.AM)                            # Extinction factor at sea level for this airmass
df['DNI']       = df.I_ext / df.extFac                                   # DNI for a clear sky

df['cosTheta']  = se.cos_angle_sun_panels(spAz,spIncl, df.sunAz,df.sunAlt)  # cos of the angle with which Sun hits my panels
df['dirRad']    = df.DNI * df.cosTheta                                   # Insolation of direct sunlight on my panels

df.sunAz  *= r2d  # Convert azimuth and ...
df.sunAlt *= r2d  # ... altitude to degrees for printing

print(df[df.sunAlt > 0])  # Print the results for the hours where the Sun is above the horizon
```


### Code example for a single calculation

Note that in most cases, the vector option is preferred (see the [code example](#code-example) above, and see
[Performance](#performance) for details).  The code example below is given for completeness.

```python
"""Example Python script using the SolarEnergy module for a single instance."""

# Location of my solar panels:
from solarenergy import d2r,r2d  # Convert between radians and degrees
geoLon =  5*d2r  # Geographic longitude (>0 for eastern hemisphere; ° -> rad)
geoLat = 52*d2r  # Geographic latitude  (>0 for northern hemisphere; ° -> rad)

# Orientation of my solar panels:
spAz   = -2*d2r  # Azimuth ('wind direction') of my panels are facing.  Note: South=0, W=90° (pi/2 rad) in the northern hemisphere!  (rad)
spIncl = 28*d2r  # Inclination of my panels w.r.t. the horizontal  (rad)

# An hour past noon local time on 1 March 2020:
myTZ  = 'Europe/Amsterdam'
year  = 2020
month = 3
day   = 1
hour  = 13


# Compute Sun position (uses SolTrack behind the scenes):
import solarenergy as se
sunAz,sunAlt,sunDist = se.sun_position_from_date_and_time(geoLon,geoLat, year,month,day, hour, timezone=myTZ)

I_ext     = 1361.5 / sunDist**2                                 # Extraterrestrial radiation (at the top of the atmosphere; AM0)

AM        = se.airmass(sunAlt)                                  # Air mass for this Sun altitude
extFac    = se.extinction_factor(AM)                            # Extinction factor at sea level for this airmass
cosTheta  = se.cos_angle_sun_panels(spAz,spIncl, sunAz,sunAlt)  # cos of the angle with which Sun hits my panels

DNI       = I_ext / extFac                                      # DNI for a clear sky
dirRad    = DNI * cosTheta                                      # Insolation of direct sunlight on my panels


# Print input and output:
import numpy as np
print("Location:           %0.3lf E, %0.3lf N"  % (geoLon*r2d, geoLat*r2d))
print("Date:               %4d-%2.2d-%2.2d"     % (year, month, day))
print("Time:               %2d:00"              % (hour))
print()

print("Sun azimuth:        %7.3lf°"   % (sunAz*r2d))
print("Sun altitude:       %7.3lf°"   % (sunAlt*r2d))
print("Sun distance:       %7.4lf AU" % (sunDist))
print()

print("I_ext:              %7.1lf W/m²"    % (I_ext))
print()

print("Air mass:           %7.3lf"         % (AM))
print("Extinction factor:  %7.3lf"         % (extFac))
print("DNI:                %7.1lf W/m²"    % (DNI))
print()

print("Sun-panels angle:   %7.1lf°"        % (np.arccos(cosTheta)*r2d))
print("Direct insolation:  %7.1lf W/m²"    % (dirRad))
```

## Performance

SolarEnergy starts with the computation of the position of the Sun for a given instant (scalar) or for a
series of instances using an array or vector.  The latter is faster than the former (for ~2 instances or more)
and may be easier to use in the majority of applications, depending on the problem given.

The table below shows the speed with which the position of the Sun is computed (in number of positions per
second) as a function of the size of the dataset.  The runs were timed on a single hyperthreaded core capped
at 3.4GHz and the minimum time of 10 runs is displayed.  The scalar runs scale linearly, the speed peaks
around 1<sup>5</sup>-1<sup>6</sup> elements for vectors.  Timings below are for timezone-naive datetimes (and
hence UTC), if timezone-aware datetimes are used, the calculations take about 4.4 times longer(!)

| Mode   | N<sub>calc</sub> | Time (s) | Speed (/s) |
|--------|------------------|----------|------------|
| scalar | 1e3              | 0.616    | 1623       |
|        |                  |          |            |
| vector | 1                | 7.95e-4  | 1258       |
| vector | 1e1              | 8.79e-4  | 11,377     |
| vector | 1e2              | 0.001037 | 96,432     |
| vector | 1e3              | 0.00257  | 389,105    |
| vector | 1e4              | 0.0134   | 746,269    |
| vector | 1e5              | 0.0687   | 1,455,604  |
| vector | 1e6              | 0.667    | 1,499,250  |
| vector | 1e7              | 8.56     | 1,168,224  |


## SolarEnergy pages

* [Pypi](https://pypi.org/project/solarenergy/): SolarEnergy Python package
* [GitHub](https://github.com/MarcvdSluys/SolarEnergy): SolarEnergy source code
* [ReadTheDocs](https://solarenergy.readthedocs.io): SolarEnergy documentation


## Author and licence

* Author: Marc van der Sluys
* Contact: http://marc.vandersluys.nl
* Licence: [EUPL 1.2](https://www.eupl.eu/1.2/en/)


## References

* This Python code is adapted from the Fortran implementation in
  [libTheSky](http://libthesky.sourceforge.net/), which contains many references.
* [Celestial mechanics in a nutshell (CMiaNS)](https://cmians.sourceforge.io/)
* [SolTrack](https://pypi.org/project/soltrack/): a free, fast and simple Python package to compute the position of the Sun, as well as its rise and set times.
