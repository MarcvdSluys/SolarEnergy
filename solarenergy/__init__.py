# -*- coding: utf-8 -*-
# SPDX-License-Identifier: EUPL-1.2
#  
#  Copyright (c) 2020-2022  Marc van der Sluys - marc.vandersluys.nl
#   
#  This file is part of the SolarEnergy Python package, containing a Python module to do simple modelling in
#  the field of solar energy.  See: https://github.com/MarcvdSluys/SolarEnergy
#   
#  SolarEnergy has been developed by Marc van der Sluys of the Department of Astrophysics at the Radboud
#  University Nijmegen, the Netherlands and the department of Sustainable energy of the HAN University of
#  applied sciences in Arnhem, the Netherlands.
#  
#  This is free software: you can redistribute it and/or modify it under the terms of the
#  European Union Public Licence 1.2 (EUPL 1.2).
#  
#  This software is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
#  without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
#  See the EU Public Licence for more details.
#  
#  You should have received a copy of the European Union Public Licence along with this code.
#  If not, see <https://www.eupl.eu/1.2/en/>.


"""SolarEnergy module

SolarEnergy contains a Python module to do simple modelling in the field of solar energy.  The code is being
developed by `Marc van der Sluys <http://marc.vandersluys.nl>`_ of the department of Astrophysics at the
Radboud University Nijmegen, the Netherlands and the department of Sustainable energy of the HAN University of
Applied Sciences in Arnhem, the Netherlands.  SolarEnergy can be used under the conditions of the EUPL 1.2
licence.  These pages contain the API documentation.  For more information on the Python package, licence,
source code and data files, see the `SolarEnergy GitHub page <https://github.com/MarcvdSluys/SolarEnergy>`_.

The SolarEnergy code is based on the `libTheSky <http://libthesky.sourceforge.net>`_ Fortran library.
Information on the theory behind this code can be found in the document
`Celestial mechanics in a nutshell <https://cmians.sourceforge.io>`_.

"""


name = 'solarenergy'

from .radiation import *
from .solar_panels import *
