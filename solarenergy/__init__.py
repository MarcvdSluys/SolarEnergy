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


"""SolarEnergy module

SolarEnergy contains a Python module to do simple modelling in the field of solar energy.  The code
is being developed by `Marc van der Sluys <http://han.vandersluys.nl/en/>`_ of the department of Astrophysics at
the Radboud University Nijmegen, the Netherlands and the department of Sustainable energy of the HAN University
of Applied Sciences in Arnhem, the Netherlands.  SolarEnergy can be used under the conditions of the GPLv3 licence.
These pages contain the API documentation.  For more information on the Python package, licence, source code
and data files, see the `SolarEnergy GitHub page <https://github.com/MarcvdSluys/SolarEnergy>`_.

The SolarEnergy code is based on the `libTheSky <http://libthesky.sourceforge.net>`_ Fortran library.
Information on the theory behind this code can be found in the document
`Celestial mechanics in a nutshell <https://cmians.sourceforge.io>`_.

"""


name = 'solarenergy'

from .radiation import *
