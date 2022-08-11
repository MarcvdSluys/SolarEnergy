#!/bin/env python
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


"""Functions for solar energy dealing with solar (PV) panels/modules.
"""


# Allow relative imports from __main__() when running this file (PEP 366):
if __name__ == '__main__' and __package__ is None:
    __package__ = 'solarenergy'


from dataclasses import dataclass


@dataclass
class SolarPanels:
    """Dataclass containing solar-panel parameters."""
    
    # Geographic location of the solar panels:
    geo_lon:  float = 0.0;  """Geographic longitude of the panels (rad; >0 for northern hemisphere)"""
    geo_lat:  float = 0.0;  """Geographic latitude of the panels (rad; >0 for east of Greenwich)"""
    
    # Orientation of the solar panels:
    az:       float = 0.0;  """'Azimuth' of the panel normal vector  (rad; 0=S, π/2=W)"""
    incl:     float = 0.0;  """'Zenith angle' of the panel normal vector  (rad; 0=horizontal, π/2=vertical)"""
    
    # Size and capacity of the solar panels:
    area:     float = 0.0;  """Surface area of solar panels (m2)"""
    eff:      float = 0.0;  """Efficiency of solar panels (0-1)"""
    t_coef:   float = 0.0;  """PV temperature coefficient (/K; typically -0.005)"""
    p_max:    float = 0.0;  """Maximum electrical power of solar panels or inverter (W)"""
    

# Test code:
if __name__ == '__main__':
    # Solar-panel data:
    _sp = SolarPanels()
