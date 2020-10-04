# SolarEnergy #

A Python module to do simple modelling in the field of solar energy.  The code is being developped by [Marc
van der Sluys](http://han.vandersluys.nl/en/) of the department of astrophysics of the Radboud University
Nijmegen, the Netherlands and the department of Sustainable energy of the HAN University of Applied Sciences
in Arnhem, the Netherlands.


## Installation ##

This package can be installed using `pip install solarenergy`.  This should automatically install the
dependency packages `pytz`, `numpy` and `soltrack` if they haven't been installed already.  If you are
installing by hand, ensure that these packages are installed as well.


## Example use ##

```python
import solarenergy as se
import numpy as np

r2d = 180/np.pi  # Multiplication factor to convert radians to degrees
d2r = 1/r2d      # Multiplication factor to convert degrees to radians

# Location of my solar panels:
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
sunAz,sunAlt,sunDist = se.computeSunPos(geoLon,geoLat, year,month,day, hour, timezone=myTZ)

AM        = se.airmass(sunAlt)                               # Air mass for this Sun altitude
extFac    = se.extinctionFactor(AM)                          # Extinction factor at sea level for this airmass
cosTheta  = se.cosAngleSunPanels(spAz,spIncl, sunAz,sunAlt)  # cos of the angle with which Sun hits my panels

solConst  = 1361.5 / sunDist**2                              # Solar constant, scaled with solar distance
DNI       = solConst / extFac                                # DNI for a clear sky
dirRad    = DNI * cosTheta                                   # Insolation of direct sunlight on my panels


# Print input and output:
print("Location:           %0.3lf E, %0.3lf N"  % (geoLon*r2d, geoLat*r2d))
print("Date:               %4d-%2.2d-%2.2d"     % (year, month, day))
print("Time:               %2d:00"              % (hour))
print()

print("Sun azimuth:        %7.3lf°"   % (sunAz*r2d))
print("Sun altitude:       %7.3lf°"   % (sunAlt*r2d))
print("Sun distance:       %7.4lf AU" % (sunDist))
print()

print("Air mass:           %7.3lf"         % (AM))
print("Extinction factor:  %7.3lf"         % (extFac))
print("Sun-panels angle:   %7.1lf°"        % (np.arccos(cosTheta)*r2d))
print()

print("Solar constant:     %7.1lf W/m²"    % (solConst))
print("DNI:                %7.1lf W/m²"    % (DNI))
print("Direct insolation:  %7.1lf W/m²"    % (dirRad))
print()
```

## SolarEnergy pages ##

* [Pypi](https://pypi.org/project/solarenergy/): SolarEnergy Python package
* [GitHub](https://github.com/MarcvdSluys/SolarEnergy): SolarEnergy source code


## Author and licence ##

* Author: Marc van der Sluys
* Contact: http://han.vandersluys.nl/en/
* Licence: [GPLv3+](https://www.gnu.org/licenses/gpl.html)


## References ##

* This Python code is adapted from the Fortran implementation in
  [libTheSky](http://libthesky.sourceforge.net/), which contains many references.
* [Celestial mechanics in a nutshell (CMiaNS)](https://cmians.sourceforge.io/)
