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
    

def pv_cell_temperature(temp_a, glob_insol, v_wind,  temp_c_noct=45, eta_c_noct=0.15, t_coef=-0.0045, glob_insol_noct=800, temp_a_noct=20, v_wind_noct=1, tau_alp=0.9):
    """Estimate the PV-cell temperature as a function of ambient temperature, insolation and wind velocity.

    Parameters:
      temp_a (float):           Ambient temperature (°C or K - ensure all temperatures use the same unit).
      glob_insol (float):       Global projected insolation on PV cell (W/m2).
      v_wind (float):           Wind velocity (m/s).
    
      temp_c_noct (float):      Cell temperature for nominal operating cell temperature (NOCT; °C or K).
      eta_c_noct (float):       Cell/module efficiency for NOCT (-).
      t_coef (float):           Temperature coefficient for Pmpp and η_cell (/K).
      glob_insol_noct (float):  Projected global insolation on the PV module for NOCT (W/m2).
      temp_a_noct (float):      Ambient temperature for NOCT (°C or K).
      v_wind_noct (float):      Wind velocity for NOCT (m/s).
    
      tau_alp (float):          Tau alpha (τα): the optical transmission-absorption coefficient for the PV module (0-1).
    
    Returns:
      (float):  PV cell temperature (°C or K).
    
    Note:
      - This follows Duffie & Beckman (2013).
      - NOCT stands for normal-operation cell temperature and concerns the specifications of the PV module for 
        realistic conditions.  If not available, use STC specs everywhere instead.
      - Note that all temperatures should be either in °C or K; they should not be mixed.
    """
    
    # Heat loss to the environment, current vs. NOCT, affected by wind:
    u_l_fac = (5.7+3.8*v_wind_noct)/(5.7+3.8*v_wind)  # U_L,N / U_L
    
    temp_c = (
        temp_a + ( temp_c_noct - temp_a_noct ) * u_l_fac * glob_insol/glob_insol_noct
        * ( 1 - eta_c_noct/tau_alp * ( 1 - t_coef * temp_c_noct ) )
    ) / ( 
        1 + ( temp_c_noct - temp_a_noct ) *  u_l_fac * glob_insol/glob_insol_noct * (t_coef * eta_c_noct)/tau_alp
    )
    
    return temp_c


def pv_efficiency(temp_c, eta_c_stc=0.15, t_coef=-0.0045, temp_c_stc=25):
    """Compute the instantaneous PV effiency for a given cell temperature.
    
    We assume that the efficiency varies linearly with temperature.
    
    Parameters:
      temp_c (float):      Temperature of the PV cell (°C or K).
      eta_c_stc (float):   Efficiency of PV (+inverter if desired) for standard test conditions (STC) (0-1).
      t_coef (float):      Temperature coefficient for P_mpp and η (/K; <0).
      temp_c_stc (float):  Temperature of the cells under Standard Test Conditions (°C or K, same as temp_c).
     
    Returns:
      float:           Efficiency of the solar cell.
    """
    
    eta_c = eta_c_stc * (1 + t_coef * (temp_c - temp_c_stc))
    
    return eta_c


# Test code:
if __name__ == '__main__':
    # Solar-panel data:
    _sp = SolarPanels()
