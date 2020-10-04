#  Copyright (c) 2020  Marc van der Sluys - marc.vandersluys.nl
#  
#  This file is part of the SolarEnergy Python package, containing a Python module to do simple modelling in
#  the field of solar energy.  See: https://github.com/MarcvdSluys/SolarEnergy
#  
#  This is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License
#  as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
#  
#  This software is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied
#  warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.
#  
#  You should have received a copy of the GNU General Public License along with this code.  If not, see
#  <http://www.gnu.org/licenses/>.


"""Constants for solar energy.

"""


# Modules:
import math as m

pi   = m.pi;           """pi"""
pi2  = 2*pi;           """2 pi"""
pio2 = pi/2;           """pi/2"""

r2d  = m.degrees(1);   """Radians to degrees"""
d2r  = m.radians(1);   """Degrees to radians"""

solConst = 1361.5;     """Solar constant: ~1361-1362 W/m² - https://en.wikipedia.org/wiki/Solar_constant"""



