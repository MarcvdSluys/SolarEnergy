#!/bin/env python
# -*- coding: utf-8 -*-

#  Copyright (c) 2020-2021  Marc van der Sluys - marc.vandersluys.nl
#  
#  This file is part of the SolarEnergy Python package, containing a Python module to do simple modelling in
#  the field of solar energy.  See: https://github.com/MarcvdSluys/SolarEnergy
#  
#  SolarEnergy has been developed by Marc van der Sluys of the Department of Astrophysics at the Radboud
#  University Nijmegen, the Netherlands and the department of Sustainable energy of the HAN University of
#  applied sciences in Arnhem, the Netherlands.
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
import numpy as np

pi   = np.pi;           """pi"""
pi2  = 2*pi;            """2 pi"""
pio2 = pi/2;            """pi/2"""

r2d  = np.degrees(1);   """Radians to degrees"""
d2r  = np.radians(1);   """Degrees to radians"""

solConst = 1361.5;      """Obsolescent; use sol_const instead!"""
sol_const = 1361.5;     """Solar constant: ~1361-1362 W/m² - Wikipedia: https://en.wikipedia.org/wiki/Solar_constant"""


# Test code:
if __name__ == '__main__':
    print(pi, r2d, sol_const)
