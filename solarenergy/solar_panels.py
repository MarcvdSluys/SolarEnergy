#!/bin/env python
# -*- coding: utf-8 -*-
# SPDX-License-Identifier: EUPL-1.2
#  
#  Copyright (c) 2020-2023  Marc van der Sluys - marc.vandersluys.nl
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
    geo_lon:  float =    0.0;  """Geographic longitude of the panels (rad; >0 for northern hemisphere)"""
    geo_lat:  float =    0.0;  """Geographic latitude of the panels (rad; >0 for east of Greenwich)"""
    tz:       str   =    None;  """Time zone where the solar panels sit (e.g. Europe/Berlin)"""
    
    # Orientation of the solar panels:
    az:       float =    0.0;  """'Azimuth' of the panel normal vector  (rad; 0=S, π/2=W)"""
    incl:     float =    0.0;  """'Zenith angle' of the panel normal vector  (rad; 0=horizontal, π/2=vertical)"""
    
    # Size and capacity of the solar panels:
    area:     float =    0.0;  """Surface area of solar PV panels (m2; typically 1.6m2 per panel)"""
    p_max:    float =    0.0;  """Maximum electrical power of solar PV panels or inverter (kW)"""
    
    # Time dependence of efficiency:
    eff0:     float =    0.0;  """Original efficiency of solar panels + inverter, at installation (0-1; e.g. 0.15 for 15%)"""
    deff_dt:  float =    0.0;  """Linear degradation of efficiency factor over time (yr^-1; -5e-3: degrades to 90% after 20 years)"""
    year:     float =   2015;  """Installation year (e.g. 2015.25 for 2015-04-01)"""
    
    # Other parameters: temperature and angle dependence:
    t_coef:   float = -0.005;  """PV temperature coefficient (/K; typically -0.005)"""
    n_refr:   float =   1.43;  """Refractive index of PV cover (typically 1.43; air: 1.000293)"""
    
    # Inverter model and serial number:
    inv_model:  str =   None;  """Model or type of the inverter"""
    inv_sn:     str =   None;  """Serial number of the inverter"""
    name:       str =   None;  """PV plant name"""
    

def read_solar_panel_specs(cfg_file='.solar_panels.cfg', rel_to_home=True, to_rad=True):
    """Read solar-panel specifications from a configuration file.
    
    Parameters:
      cfg_file (str):      Configuration file to read specs from (by default relative to home directory).
      rel_to_home (bool):  File path/name is relative to home directory.
      to_rad (bool):       Convert angles from degrees to radians (location, orientation).
    
    Returns:
      (SolarPanels):  Dataclass of type SolarPanels containing the specifications.
    """
    
    sp = SolarPanels()
    
    
    # Read configuration file:
    import configparser
    config = configparser.ConfigParser(inline_comment_prefixes=('#'))  # Allow inline comments
    if rel_to_home:
        from pathlib import Path as _Path
        config.read(str(_Path.home())+'/'+cfg_file)
    else:
        config.read(str(cfg_file))
    
    # Use fallback to allow missing sections and keys, and use the default values instead.
    
    # Section Geographic location of the solar panels:
    sp.geo_lon  = config.getfloat('Location', 'geo_lon', fallback=sp.geo_lon)         # Geographic longitude of the panels (rad; >0 for northern hemisphere)
    sp.geo_lat  = config.getfloat('Location', 'geo_lat', fallback=sp.geo_lat)         # Geographic latitude of the panels (rad; >0 for east of Greenwich)
    sp.tz       = config.get('Location',     'timezone', fallback=sp.tz)              # Timezone where the solar panels sit (e.g. Europe/Paris)
    
    # Section Orientation of the solar panels:
    sp.az       = config.getfloat('Orientation', 'az',   fallback=sp.az)              # 'Azimuth' of the panel normal vector  (rad; 0=S, π/2=W)
    sp.incl     = config.getfloat('Orientation', 'incl', fallback=sp.incl)            # 'Zenith angle' of the panel normal vector  (rad; 0=horizontal, π/2=vertical)
    
    # Section Size and capacity of the solar panels:
    sp.area     =  config.getfloat('Capacity', 'area',   fallback=sp.area)            # Surface area of solar PV panels (m2; typically 1.6m2 per panel)
    sp.p_max    =  config.getfloat('Capacity', 'p_max',  fallback=sp.p_max)           # Maximum electrical power of solar PV panels or inverter (kW)
    
    # Section Time dependence of efficiency:
    sp.eff0     =  config.getfloat('TimeDependence', 'eff0',    fallback=sp.eff0)     # Original efficiency of solar panels + inverter, at installation (0-1; e.g. 0.15 for 15%)
    sp.deff_dt  =  config.getfloat('TimeDependence', 'deff_dt', fallback=sp.deff_dt)  # Linear degradation of efficiency factor over time (yr^-1; -5e-3: degrades to 90% after 20 years)
    sp.year     =  config.getfloat('TimeDependence', 'year',    fallback=sp.year)     # Installation year (e.g. 2015.25 for 2015-04-01)
    
    # Section Other parameters: temperature and angle dependence:
    sp.t_coef   =  config.getfloat('Other', 't_coef',  fallback=sp.t_coef)            # PV temperature coefficient (/K; typically -0.005)
    sp.n_refr   =  config.getfloat('Other', 'n_refr',  fallback=sp.n_refr)            # Refractive index of PV cover (typically 1.43; air: 1.000293)
    
    # Section Inverter: model and serial number:
    sp.inv_model =  config.get('Inverter', 'model',    fallback=sp.inv_model)         # Inverter model or type
    sp.inv_sn    =  config.get('Inverter', 'sn',       fallback=sp.inv_sn)            # Inverter serial number
    sp.name      =  config.get('Inverter', 'name',     fallback=sp.name)              # PV plant name
    
    
    if to_rad:
        from astroconst import d2r
        sp.geo_lon *= d2r
        sp.geo_lat *= d2r
        
        sp.az   *= d2r
        sp.incl *= d2r
        
    return sp


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
