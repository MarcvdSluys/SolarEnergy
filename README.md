# SolarEnergy #

![PyPI](https://img.shields.io/pypi/v/solarenergy?color=%230A0)
![PyPI - Downloads](https://img.shields.io/pypi/dm/solarenergy)
[![Code check](https://github.com/MarcvdSluys/SolarEnergy/actions/workflows/code-check.yml/badge.svg)](https://github.com/MarcvdSluys/SolarEnergy/actions/workflows/code-check.yml)
[![Documentation
Status](https://readthedocs.org/projects/solarenergy/badge/?version=latest)](https://solarenergy.readthedocs.io/en/latest/?badge=latest)
![PyPI - License](https://img.shields.io/pypi/l/solarenergy?color=%230A0)

A Python module to do simple modelling in the field of solar energy.  The code is being developed by [Marc
van der Sluys](http://han.vandersluys.nl/en/) of the department of astrophysics of the Radboud University
Nijmegen, the Netherlands and the department of Sustainable energy of the HAN University of Applied Sciences
in Arnhem, the Netherlands.  The SolarEnergy package can be used under the conditions of the GPLv3
licence.


## Installation ##

This package can be installed using `pip install solarenergy`.  This should automatically install the
dependency packages `pytz`, `numpy` and `soltrack` if they haven't been installed already.  If you are
installing by hand, ensure that these packages are installed as well.


## Example use ##

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

## SolarEnergy pages ##

* [Pypi](https://pypi.org/project/solarenergy/): SolarEnergy Python package
* [GitHub](https://github.com/MarcvdSluys/SolarEnergy): SolarEnergy source code
* [ReadTheDocs](https://solarenergy.readthedocs.io): SolarEnergy documentation


## Author and licence ##

* Author: Marc van der Sluys
* Contact: http://han.vandersluys.nl/en/
* Licence: [GPLv3+](https://www.gnu.org/licenses/gpl.html)


## References ##

* This Python code is adapted from the Fortran implementation in
  [libTheSky](http://libthesky.sourceforge.net/), which contains many references.
* [Celestial mechanics in a nutshell (CMiaNS)](https://cmians.sourceforge.io/)
