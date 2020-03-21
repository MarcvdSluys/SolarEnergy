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


"""SolarEnergy module

SolarEnergy contains a Python module do simple modelling in the field of solar energy.  SolarEnergy can be used
under the conditions of the GPLv3 licence.  These pages contain the API documentation.  For more information
on the Python package, licence, source code and data files, see the 
[SolarEnergy GitHub page](https://github.com/MarcvdSluys/SolarEnergy).

The SolarEnergy code is based on the [libTheSky](http://libthesky.sourceforge.net) Fortran library.
Information on the theory behind this code can be found in the document
[Celestial mechanics in a nutshell](https://cmians.sourceforge.io).

"""


name = "solarenergy"

from .radiation import *



