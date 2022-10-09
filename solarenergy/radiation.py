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


"""Functions for solar energy dealing with (solar) radiation.

References:
  - Marc van der Sluys, Celestial mechanics in a nutshell, https://cmians.sourceforge.io (2014-2022).
"""


# Allow relative imports from __main__() when running this file (PEP 366):
if __name__ == '__main__' and __package__ is None:
    __package__ = 'solarenergy'


import numpy as np
from astroconst import pi2,pio2, d2r,r2d, sol_const


def sun_position_from_date_and_time(geo_lon,geo_lat, year,month,day, hour,minute=0,second=0, timezone='UTC', debug=False):
    """Compute the Sun local position (azimuth, altitude and distance) for the given geographical location and
    date and time (ymd, hms) using SolTrack.
    
    Parameters:
        geo_lon  (float):  Geographic longitude to compute the Sun position for (rad).
        geo_lat  (float):  Geographic latitude to compute the Sun position for (rad).
    
        year     (int):    Year (CE) to compute the Sun position for.
        month    (int):    Month to compute the Sun position for.
        day      (int):    Day of month to compute the Sun position for.
    
        hour     (int):    Hour of day to compute the Sun position for (local time!).
        minute   (int):    Minute to compute the Sun position for (optional; default = 0).
        second   (int):    Second to compute the Sun position for (optional; default = 0).
    
        timezone (str):    Timezone for which date and time are provided (optional; default = 'UTC').
    
        debug    (bool):   Switch to write detailed output to stdout (optional; default = False).
    
    Returns:
        tuple (float,float,float):  Tuple containing (azimuth, altitude, distance):
    
            - azimuth  (float):  Azimuth of the Sun (rad; south = 0 on the northern hemisphere).
            - altitude (float):  Altitude of the Sun (rad).
            - distance (float):  Distance Sun-Earth (AU).
    """
    
    # Create a timezone-aware datetime object:
    import datetime as dt
    import pytz as tz
    myTZ   = tz.timezone(timezone)   # My timezone
    myTime = dt.datetime(int(year),int(month),int(day), int(hour),int(minute),int(second))  # Timezone-naive time
    lt     = myTZ.localize(myTime, is_dst=None)  # Mark as local time
    utc    = lt.astimezone(tz.utc)               # Convert to UTC
    
    azimuth, altitude, distance = sun_position_from_datetime(geo_lon,geo_lat, utc, debug)
    
    return azimuth, altitude, distance



def sun_position_from_datetime(geo_lon,geo_lat, date_time, utc=False, debug=False):
    """Compute the Sun local position (azimuth, altitude and distance) for the given geographical location and
    a Python datetime using SolTrack.
    
    Parameters:
        geo_lon (float):       Geographic longitude to compute the Sun position for (rad).
        geo_lat (float):       Geographic latitude to compute the Sun position for (rad).
    
        date_time (datetime):  Date and time to compute the Sun position for (timezone aware and/or UTC, i.e. must be UTC if timezone naive.  CHECK: must be UTC if a list or (numpy) array?).
        
        utc (bool):            Specify that the input is in UTC.  This can speed up calls with a single datetime (as opposed to arrays).  Defaults to False.
        debug (bool):          Switch to write detailed output to stdout (optional; default = False).
    
    Returns:
        tuple (float,float,float):  Tuple containing (arrays of) (azimuth, altitude, distance):
    
            - azimuth  (float):  Azimuth of the Sun (rad; south = 0 on the northern hemisphere).
            - altitude (float):  Altitude of the Sun (rad).
            - distance (float):  Distance Sun-Earth (AU).
    """
    
    # Create a SolTrack instance for the desired location and specify preferences:
    from soltrack import SolTrack
    st = SolTrack(geo_lon,geo_lat, compute_refr_equatorial=False)  # No need for equatorial coordinates
    
    # Set date and time:
    st.set_date_time(date_time, utc=utc)
    
    # Compute the Sun's position:
    st.compute_position()
    
    
    if debug:
        print('Location:  ', st.geoLongitude*r2d, st.geoLatitude*r2d)
        print('Date/time: ', date_time)
        print('JD:        ', st.julianDay)
        print()
        
        print('Azimuth:   ', st.azimuth*r2d)
        print('Altitude:  ', st.altitude*r2d)
        print('Distance:  ', st.distance)
        print()
        
    return st.azimuth, st.altitude, st.distance



def cos_angle_sun_panels(sp_az,sp_incl, sun_az,sun_alt):
    """Compute the cosine of the angle between the orientation vector of the solar panels and the position vector
       of the Sun.
    
    This is the cosine of the angle under which the direct sunlight hits the solar panels.  Multiply it with
    the DNI to obtain the direct insolation on the panels.  See Celestial mechanics in a nutshell, Sect. 4.3:
    Insolation on an inclined surface (http://CMiaNS.sf.net).
    
    Parameters:
        sp_az   (float):  Azimuth in which the solar panels are facing (rad; e.g. north or south = 0, same as sun_az).
        sp_incl (float):  Inclination ('zenith angle') of the solar panels w.r.t. the horizontal (rad).
        sun_az  (float):  Azimuth of the Sun (rad; e.g. north or south = 0, same as sp_az).
        sun_alt (float):  Altitude of the Sun (rad).
    
    Returns:
        float:  The cosine between the normal vector of the solar panels and the position vector of the Sun (rad).
                Note that this value is zero (indicating radiation from behind the panels) or positive.
                
    """
    
    cos_theta = np.sin(sun_alt) * np.cos(sp_incl)  +  np.cos(sun_alt) * np.sin(sp_incl) * np.cos(sun_az - sp_az)
    cos_theta = np.maximum(cos_theta, 0)  # Return 0 rather than negative values (indicating radiation from behind).
    
    return cos_theta


def airmass(sun_alt, return_value_below_horizon=False):
    """Compute airmass as a function of Sun altitude.
    
    Parameters:
        sun_alt                    (float):  True altitude of the Sun (uncorrected for atmospheric refraction; rad), can be an array.
        return_value_below_horizon (bool):   Return a very large value when the Sun is below the horizon, larger
                                             when the Sun is lower.  This can be useful for solvers.  Default: False.
    
    Returns:
        float:  Airmass at sea level (AM~1 if the Sun is in the zenith, AM~38 near the horizon).
    
    References:
        A.T. Young, "AIR-MASS AND REFRACTION," Applied Optics, vol. 33, pp. 1108-1110, Feb 1994.
    """
    
    sun_alt = np.asarray(sun_alt)
    scalar_input = False
    if sun_alt.ndim == 0:
        sun_alt = sun_alt[np.newaxis]  # Convert scalar to 1D array
        scalar_input = True
    
    airmass = np.empty(sun_alt.shape)
    
    # Sun below the horizon:
    sel = (sun_alt < -0.00989)
    if return_value_below_horizon:
        airmass[sel] = 1000 * (0.15 + abs(sun_alt[sel]))  # Very bad, but still getting worse for even lower Sun, for solvers
    else:
        airmass[sel] = float('inf')
    
    # Sun above the horizon:
    sel = (sun_alt >= -0.00989)
    airmass[sel] = (1.002432 * np.sin(sun_alt[sel]) ** 2 + 0.148386 * np.sin(sun_alt[sel]) + 0.0096467) / \
        (np.sin(sun_alt[sel]) ** 2 * np.sin(sun_alt[sel]) + 0.149864 * np.sin(sun_alt[sel]) ** 2 + 0.0102963 * np.sin(sun_alt[sel]) + 0.000303978)
    airmass[sel] = np.maximum(airmass[sel], 1)  # Air mass cannot be lower than 1
    
    if scalar_input:
        return np.ndarray.item(airmass)
    
    return airmass



def extinction_factor(airmass, return_value_below_horizon=False):
    """Compute the atmospheric extinction factor for sunlight from the air mass.
    
    Parameters:
        airmass                    (float):  Airmass at sea level (AM~1 if the Sun is in the zenith, AM~38 near the horizon), can be an array.
        return_value_below_horizon (bool):   Return a very large value when the Sun is below the horizon, larger
                                             when the Sun is lower.  This can be useful for solvers.  Default: False.
    
    Returns: 
        float: The extinciton factor for sunlight in the atmosphere.  Divide the extraterrestrial (AM0)
               radiation (or, if unknown, the solar constant) by this number to obtain the DNI.
    
    """
    
    coefs = [ 9.1619283e-2,  2.6098406e-1,-3.6487512e-2,  6.4036283e-3,-8.1993861e-4, 6.9994043e-5,-3.8980993e-6,
              1.3929599e-7, -3.0685834e-9, 3.7844273e-11,-1.9955057e-13]  # Fit coefficients
    
    airmass = np.asarray(airmass)
    scalar_input = False
    if airmass.ndim == 0:
        airmass = airmass[np.newaxis]  # Convert scalar to 1D array
        scalar_input = True
    
    ext_fac = np.empty(airmass.shape)
    
    # Sun below the horizon:
    sel = (airmass > 38.2)
    if return_value_below_horizon:
        import sys
        ext_fac[sel] = np.sqrt(sys.float_info.max) * (0.15 + airmass[sel])  # Very bad, but still getting worse for even higher airmass, for solvers
    else:
        ext_fac[sel] = float('inf')
    
    # Sun above the horizon:
    AMpow = np.ones(airmass.shape)             # AM^0 = 1
    ext = np.ones(airmass.shape) * coefs[0]    # c_1 * AM^0
    sel = (airmass <= 38.2)
    for iCoef in range(1, len(coefs)):
        AMpow[sel] *= airmass[sel]             # AM^(i-1)
        ext[sel] += coefs[iCoef] * AMpow[sel]  # + c_i * AM^(i-1)
    
    ext_fac[sel] = np.exp(ext[sel])
    
    if scalar_input:
        return np.ndarray.item(ext_fac)
    
    return ext_fac


def diffuse_radiation_projection_perez87(doy, sun_alt, surf_incl, theta, beam_norm, dif_horiz, return_components=False):
    """Compute diffuse radiation on an inclined surface using the 1987 Perez model
    
    This function is adapted from the libTheSky Fortran implementation (libthesky.sf.net).
    
    See Perez et al. Solar Energy Vol. 39, Nr. 3, p. 221 (1987) - references to equations and tables are to
    this paper.  Most equations can be found in the Nomenclature section at the end of the paper (p.230).  I
    use a and c here, not b and d.
    
    
    Parameters:
      doy       (int):     Day of year (Nday)
      sun_alt   (float):   Altitude of the Sun (radians, may be an array)
      
      surf_incl (float):   Surface inclination wrt horizontal (radians) - 0 = horizontal, pi/2 = vertical
      theta     (float):   Angle between surface normal vector and Sun position vector (radians, may be an array)
      
      beam_norm (float):   Beam (direct) normal radiation = DNI (W/m2; in the direction of the Sun, may be an array)
      dif_horiz (float):   Diffuse radiation on a horizontal surface = DHI (W/m2, may be an array)

      return_components (bool):  Return circumsolar and horizon parts separately (i.e., three return values).
      
    Returns:
      float:    Diffuse irradiation on the inclined surface (W/m2) (may be an array) - if return_components=False
    
      tuple:    Diffuse irradiation + components as (float1, float2, float3), may be arrays - if return_components=True:

                - float1:    Total diffuse irradiation on the inclined surface (W/m2)
                - float2:    Circumsolar diffuse irradiation on the inclined surface (W/m2)
                - float3:    Horizon-band diffuse irradiation on the inclined surface (W/m2)
      
    """
    
    # *** Compute the brightness coefficients for the circumsolar (F1) and horizon (F2) regions ***
    
    arrSize = np.size(sun_alt)  # Size (length) of the 1D numpy arrays (1 if no arrays)
    
    # 'External' (AM0) radiation:
    AM0rad = 1370 * (1 + 0.00333 * np.cos(pi2/365 * doy))
    
    # Air mass:
    if arrSize == 1:  # Scalar
        if sun_alt < -3.885*d2r:
            Mair = 99
        elif sun_alt < 10*d2r:
            Mair = 1 / ( np.sin(sun_alt) + 0.15 * (sun_alt*r2d + 3.885)**(-1.253) )
        else:
            Mair = 1 / np.sin(sun_alt)
    else:  # Array
        Mair = np.ones(arrSize) * 36.51  # Air mass is 36.51 (value for sun_alt=0), unless...
        Mair[sun_alt >= 10*d2r]    = 1 / np.sin(sun_alt[sun_alt >= 10*d2r])
        Mair[sun_alt <  10*d2r]    = 1 / ( np.sin(sun_alt[sun_alt < 10*d2r]) + 0.15 * (sun_alt[sun_alt < 10*d2r]*r2d + 3.885)**(-1.253) )
        Mair[sun_alt < -3.885*d2r] = 99
    
    Delta = dif_horiz * Mair / AM0rad  # Brightness of overcast sky - par. 2.2.4 (a)
    
    
    # Cloudliness: epsilon;  epsilon ~ 1: overcast, epsilon -> infinity: clear  (epsilon ~ 1/fraction of covered sky)
    #   Needed for correct row in Table 1
    if arrSize == 1:  # Scalar
        if dif_horiz <= 0:  # Division by zero
            if beam_norm <= 0:  # No direct light: 0/0
                epsilon = 0     # -> completely overcast - first row of Table 1
            else:               # Some direct daylight: x/0 = large
                epsilon = 99    # -> completely clear, should be >11 for last row of Table 1
        else:
            epsilon = (dif_horiz + beam_norm) / dif_horiz  # Overcast: epsilon ~ 1,  clear: epsilon -> infinity
    else:  # Array
        epsilon = np.zeros(arrSize)                        # Set epsilon = 0 by default = completely overcast - first row of Table 1
        epsilon[(dif_horiz <= 0) & (beam_norm > 0)] = 99   # Some direct daylight: x/0 = large -> completely clear, should be >11 for last row of Table 1
        epsilon[dif_horiz > 0] = (dif_horiz[dif_horiz > 0] + beam_norm[dif_horiz > 0]) / dif_horiz[dif_horiz > 0]  # Overcast: epsilon ~ 1,  clear: epsilon -> infinity
    
    
    # Table 1:
    f11=0;  f12=1;  f13=2;  f21=3; f22=4; f23=5
    
    if arrSize == 1:  # Scalar
        if epsilon <= 1.056:
            F = [-0.011,  0.748, -0.080, -0.048,  0.073, -0.024]
        elif epsilon <= 1.253:
            F = [-0.038,  1.115, -0.109, -0.023,  0.106, -0.037]
        elif epsilon <= 1.586:
            F = [ 0.166,  0.909, -0.179,  0.062, -0.021, -0.050]
        elif epsilon <= 2.134:
            F = [ 0.419,  0.646, -0.262,  0.140, -0.167, -0.042]
        elif epsilon <= 3.230:
            F = [ 0.710,  0.025, -0.290,  0.243, -0.511, -0.004]
        elif epsilon <= 5.980:
            F = [ 0.857, -0.370, -0.279,  0.267, -0.792,  0.076]
        elif epsilon <= 10.080:
            F = [ 0.734, -0.073, -0.228,  0.231, -1.180,  0.199]
        else:
            F = [ 0.421, -0.661,  0.097,  0.119, -2.125,  0.446]
            
        zeta = pio2 - sun_alt                             # Zenith angle = pi/2 - sun_alt
        F1 = F[f11]  +  F[f12] * Delta  +  F[f13] * zeta  # Circumsolar brightness coefficient
        F2 = F[f21]  +  F[f22] * Delta  +  F[f23] * zeta  # Horizon brightness coefficient
        
    else:  # Array
        F = np.empty((arrSize, 6))
        F[                      (epsilon <=  1.056), :]  =  [-0.011,  0.748, -0.080, -0.048,  0.073, -0.024]
        F[ (epsilon >  1.056) & (epsilon <=  1.253), :]  =  [-0.038,  1.115, -0.109, -0.023,  0.106, -0.037]
        F[ (epsilon >  1.253) & (epsilon <=  1.586), :]  =  [ 0.166,  0.909, -0.179,  0.062, -0.021, -0.050]
        F[ (epsilon >  1.586) & (epsilon <=  2.134), :]  =  [ 0.419,  0.646, -0.262,  0.140, -0.167, -0.042]
        F[ (epsilon >  2.134) & (epsilon <=  3.230), :]  =  [ 0.710,  0.025, -0.290,  0.243, -0.511, -0.004]
        F[ (epsilon >  3.230) & (epsilon <=  5.980), :]  =  [ 0.857, -0.370, -0.279,  0.267, -0.792,  0.076]
        F[ (epsilon >  5.980) & (epsilon <= 10.080), :]  =  [ 0.734, -0.073, -0.228,  0.231, -1.180,  0.199]
        F[ (epsilon > 10.080), :                      ]  =  [ 0.421, -0.661,  0.097,  0.119, -2.125,  0.446]
        
        zeta = pio2 - sun_alt                                      # Zenith angle = pi/2 - sun_alt
        F1 = F[:, f11]  +  F[:, f12] * Delta  +  F[:, f13] * zeta  # Circumsolar brightness coefficient
        F2 = F[:, f21]  +  F[:, f22] * Delta  +  F[:, f23] * zeta  # Horizon brightness coefficient
    
    
    
    
    # *** Compute the mean solid angles occupied by the circumsolar region (C and A, needed for Eq.8) ***
    
    alpha = 25*d2r  # Half angle of the circumsolar region (degrees -> radians; below Eq.9)
    
    
    # Solid angle of the circumsolar region weighted by incidence on the HORIZONTAL (variable C, subscript H;
    #   see Nomenclature, under c):
    # psiH:
    if arrSize == 1:  # Scalar
        if zeta > pio2 - alpha:
            psiH = 0.5 * (pio2 - zeta + alpha) / alpha  # Dimensionless ratio
        else:
            psiH = 1
    else:  # Array
        psiH = np.ones(arrSize)
        psiH[zeta > pio2 - alpha] = 0.5 * (pio2 - zeta[zeta > pio2 - alpha] + alpha) / alpha  # Dimensionless ratio
    
    # chiH:
    if arrSize == 1:  # Scalar
        if zeta < pio2 - alpha:
            chiH = np.cos(zeta)  # = np.sin(sun_alt)
        else:
            chiH = psiH * np.sin(psiH*alpha)
    else:  # Array
        chiH = np.cos(zeta)      # = np.sin(sun_alt)
        chiH[zeta >= pio2 - alpha] = psiH[zeta >= pio2 - alpha] * np.sin(psiH[zeta >= pio2 - alpha] * alpha)
    
    C = 2 * (1 - np.cos(alpha)) * chiH  # Solid angle of the circumsolar region, weighted by HORIZONTAL incidence
    
    
    # Solid angle of the circumsolar region weighted by incidence on the SLOPE (variable A, subscript C;
    #   see Nomenclature, under c):
    # psiC:
    psiC = 0.5 * (pio2 - theta + alpha) / alpha
    
    # chiC:
    if arrSize == 1:  # Scalar
        if theta < pio2 - alpha:
            chiC = psiH * np.cos(theta)
        elif theta < pio2 + alpha:
            chiC = psiH * psiC * np.sin(psiC*alpha)
        else:
            chiC = 0
    
    else:  # Array
        chiC = np.zeros(arrSize)
        chiC[theta < pio2 + alpha] = psiH[theta < pio2 + alpha] * psiC[theta < pio2 + alpha] * np.sin(psiC[theta < pio2 + alpha] * alpha)
        chiC[theta < pio2 - alpha] = psiH[theta < pio2 - alpha] * np.cos(theta[theta < pio2 - alpha])
    
    A = 2 * (1 - np.cos(alpha)) * chiC  # Solid angle of the circumsolar region, weighted by SLOPE incidence
    
    
    
    # Diffuse radiation from circumsolar (F1) and horizon (F2) regions on the inclined surface (Eq.8):
    diff_incl_csl = dif_horiz * ( 0.5 * (1 + np.cos(surf_incl)) * (1 - F1)  +  F1 * A/C )  # Circumsolar
    diff_incl_hzl = dif_horiz * ( F2 * np.sin(surf_incl) )                                 # Horizon band
    
    diff_incl_csl = np.maximum(diff_incl_csl, 0)  # Components are sometimes negative
    diff_incl_hzl = np.maximum(diff_incl_hzl, 0)
    diff_incl     = diff_incl_csl + diff_incl_hzl
    
    # Assign optional return values:
    if return_components:
        return diff_incl, diff_incl_csl, diff_incl_hzl
    else:
        return diff_incl


def clearsky_bird(sun_alt, i_ext=1353,sun_dist=1, press=1013,  uo=0.34,uw=1.42, ta5=0.2661,ta3=0.3538,ba=0.84,k1=0.1, rg=0.2):
    """A simplified clear-sky model for direct and diffuse insolation on horizontal surfaces.
    
    A.k.a. as "the Bird model".
    
    This function is adapted from the libTheSky Fortran implementation (libthesky.sf.net).
    
    See Bird & Hulstrom, A simplified clear-sky model for direct and diffuse insolation on horizontal
    surfaces, SERI/TR-642-761 (1981).
    
    Note that the value of Taa does not agree with tabulated values from the paper, and hence neither do
    dependent values (except for AM~1).  When I substitute their values for Taa, everything matches perfectly.
    Error in their formula, or (hopefully!) in their table?
    
    Parameters:
      sun_alt   (float):  Sun altitude above the horizon (rad)
    
      i_ext     (float):  Extraterrestrial radiation (at the top of the atmosphere; AM0; W/m^2 - optional, default: 1353 (1361.5))
      sun_dist  (float):  Sun distance (AU - optional, default: 1)
      press     (float):  Air pressure at the observer's site, corrected for altitude (hPa - optional, default: 1013)
    
      uo        (float):  Ozone abundance in a vertical column (cm - optional, default: 0.34)
      uw        (float):  Percipitable water-vapor abundance in a vertical column (cm - optional, default: 1.42)
    
      ta5       (float):  Aerosol optical depth from surface in vertical path at 500 nm (optional, default: 0.2661)
      ta3       (float):  Aerosol optical depth from surface in vertical path at 380 nm (optional, default: 0.3538)
      ba        (float):  Aerosol forward-scattering ratio  (optional, 0.82-0.86, default: 0.84)
      k1        (float):  Aerosol-absorptance constant (optional, rural: 0.0933, urban: 0.385, default: 0.1)
    
      rg        (float):  Ground albedo (optional, fraction - default: 0.2)
    
    
    Returns:
      tuple (float,float,float,float):  Tuple containing (i_tot, i_dir, i_dif, i_gr):
        
      - i_tot  (float):  Total (global) insolation on a horizontal surface (GHI; W/m^2)
      - i_dir  (float):  Direct (beam) insolation on a horizontal surface (BHI; W/m^2)
      - i_dif  (float):  Diffuse insolation on a horizontal surface (DHI; W/m^2)
      - i_gr   (float):  Ground-reflection insolation from a horizontal surface (W/m^2)
    """
    
    Z = pio2 - sun_alt  # Solar zenith angle
    cosZ = np.cos(Z)    # Save a few CPU cycles
    
    
    # Relative air mass for the solar vector:
    AM  = 1/(cosZ + 0.15 * (93.885-Z*r2d)**(-1.25))  # Air mass
    AMp = AM * press / 1013                          # Pressure-corrected air mass
    
    
    # TRANSMISSION EQUATIONS:
    # Rayleigh scattering:
    Tr = np.exp( -0.0903 * AMp**0.84 * (1 + AMp - AMp**1.01) )
    
    # Ozone:
    Xo = uo*AM  # Amount of ozone in the direction of the Sun
    To = 1  -  0.1611 * Xo * (1+139.48*Xo)**(-0.3035)  -  0.002715 * Xo / (1 + 0.044*Xo + 0.0003*Xo**2)  # Transmittance of ozone absorptance
    
    # Uniformly mixed gases (CO2, O2):
    Tum = np.exp(-0.0127 * AMp**0.26)  # Transmittance of mixed-gas absorptance
    
    # Water vapor:
    Xw = AM*uw  # Amount of water vapor in the direction of the Sun
    Tw = 1 - 2.4959 * Xw / ((1 + 79.034*Xw)**0.6828 + 6.385*Xw)             # Transmittance of water-vapor absorptance - Tw = 1-Aw
    
    # Daily turbidity:
    Tau = 0.2758*ta3 + 0.35*ta5                                             # Broadband turbidity: aerosol optical depth from surface in a vertical column
    Ta  = np.exp( -Tau**0.873  *  (1 + Tau - Tau**0.7088)  *  AM**0.9108 )  # Transmittance of aerosol absorptance and scattering
    Taa = 1 - k1 * (1 - AM + AM**1.06) * (1-Ta)                             # Transmittance of aerosol absorptance - this does not agree with tabulated values from the paper (except for AM~1).  When I substitute their values for Taa, everything matches perfectly.  Error in their formula, or in their table?
    Tas = Ta/Taa                                                            # Transmittance of aerosol scattering
    Rs  = 0.0685 + (1-ba) * (1-Tas)                                         # Sky albedo
    
    
    # IRRADIANCE EQUATIONS:
    # Direct radiation on a horizontal surface:
    tmpVar = i_ext * cosZ  *  To * Tum * Tw  # Save a few CPU cycles
    i_dir = 0.9662 * tmpVar  *  Tr * Ta  /  np.square(sun_dist)
    
    # Diffuse (atmosphere-scattered) radiation on a horizontal surface:
    i_dif  = 0.79 *  tmpVar        * Taa *  (0.5*(1-Tr) + ba*(1-Tas)) / (1 - AM + AM**1.02)
    
    # Total (direct + diffuse + ground->sky scattered) radiation on a horizontal surface:
    i_tot = (i_dir+i_dif) / (1 - rg*Rs)
    
    # MvdS: add ground-reflected radiation from a horizontal surface:
    i_gr  = (i_dir+i_dif)*rg
    
    return i_tot, i_dir, i_dif, i_gr


def diffuse_radiation_from_global_radiation_and_sunshine(glob_horiz, sun_frac, sun_alt, i_ext=sol_const):
    """Compute the diffuse horizontal radiation from the global horizontal radiation, the fraction of sunshine
    and the Sun altitude.
    
    Parameters:
      glob_horiz (float):  Global horizontal radiation (W/m2).
      sun_frac   (float):  Fraction of sunshine (e.g. fraction of cloud cover) (-; 0-1).
      sun_alt    (float):  Sun altitude above the horizon (rad).
      i_ext      (float):  Extraterrestrial radiation (W/m2).  Defaults to the solar constant.
    
    Returns:
      tuple (float,float,float):  Tuple containing (dif_horiz, beam_horiz, beam_norm):
      
        - dif_horiz  (float):  Diffuse horizontal radiation = DHI (W/m2).
        - beam_horiz (float):  Beam (direct) horizontal radiation = BHI (W/m2).
        - beam_norm  (float):  Beam (direct) normal radiation = BNI = DNI (W/m2).
    """
    
    beam_norm  = i_ext / extinction_factor(airmass(sun_alt)) * sun_frac  # (Mean) DNI
    beam_horiz = beam_norm*np.sin(sun_alt)                               # Beam horizontal radiation
    dif_horiz  = glob_horiz - beam_horiz                                 # Diffuse horizontal radiation
    
    return dif_horiz, beam_horiz, beam_norm


def reflectance_transmittance(ang_i, n_ref1,n_ref2, comp_transmittance=False,comp_polarised=False):
    """Compute the reflectance and transmittance as function of incidence angle and refactive indices.
    
    Compute the reflectance for the transition from a medium with refractive index n_ref1 to one with n_ref2,
    under an incidence angle ang_i.  Optionally, compute the transmittance, and the polarised components.  The
    media are assumed to be non-magnetic.
    
    Parameters:
      ang_i  (float):  Angle of incidence (rad).
      n_ref1 (float):  Refractive index of initial medium (-).
      n_ref2 (float):  Refractive index of second medium (-).
    
      comp_transmittance (bool):  Compute and return the transmittance and its angle.
      comp_polarised (bool):      Compute and return the polarised reflectances and transmittances.
    
    Returns:
      (tuple):  Tuple consisting of one or more values, depending on the input parameters:
                
                - comp_transmittance=False, comp_polarised=False:  (r_unp);
                - comp_transmittance=False, comp_polarised=True:   (r_unp, r_prp,r_par);
                - comp_transmittance=True, comp_polarised=False:   (r_unp, t_unp, ang_t);
                - comp_transmittance=True, comp_polarised=True:    (r_unp, t_unp, ang_t,  r_prp,r_par,
                                                                   t_prp,t_par);
                
                With the following variables:
                  r_unp (float):   Unpolarised reflectance (-);
                  
                  t_unp (float):   Unpolarised transmittance (-);
                  ang_t (float):   Angle of transmittance (rad);
                  
                  r_prp (float):   Perpendicular polarised reflectance (-);
                  r_par (float):   Parallel polarised reflectance (-);
                  t_prp (float):   Perpendicular polarised transmittance (-);
                  t_par (float):   Parallel polarised transmittance (-).
    
    See:
      - libSUFR, optics.f90: libsufr.sf.net
      - Hecht, Optics, 3rd Ed. (1998), p.113ff
      - https://en.wikipedia.org/wiki/Fresnel_equations#Power_or_intensity_equations
    """
    
    var = n_ref1/n_ref2 * np.sin(ang_i)  # Argument for Snell's law
    
    # Default values (for total internal reflection):
    r_prp = np.ones_like(ang_i)
    r_par = np.ones_like(ang_i)
    ang_t = np.zeros_like(ang_i)
    
    SEL = (var <= 1) & (np.abs(ang_i) <= pio2)  # Selection of no total internal reflection and a valid input value
    ang_t[SEL] = np.arcsin(var[SEL])            # Angle of transmittance - Snell's law
    
    # Reuse variables:
    cos_ang_i = np.cos(ang_i)
    cos_ang_t = np.cos(ang_t)
    
    r_prp[SEL] = np.square( (n_ref1 * cos_ang_i[SEL] - n_ref2 * cos_ang_t[SEL]) /
                            (n_ref1 * cos_ang_i[SEL] + n_ref2 * cos_ang_t[SEL]) )
    r_par[SEL] = np.square( (n_ref1 * cos_ang_t[SEL] - n_ref2 * cos_ang_i[SEL]) /
                            (n_ref1 * cos_ang_t[SEL] + n_ref2 * cos_ang_i[SEL]) )
    
    
    # Unpolarised reflectance:
    r_unp = 0.5 * (r_prp + r_par)
    
    # Assign optional return values:
    t_unp = 1 - r_unp
    t_prp = 1 - r_prp
    t_par = 1 - r_par
    
    # print(ang_i, var, SEL, r_prp,r_par, r_unp)
    
    # Number of return values depends on input parameters:
    if comp_transmittance:
        if comp_polarised:
            return r_unp, t_unp, ang_t,  r_prp,r_par, t_prp,t_par
        else:
            return r_unp, t_unp, ang_t
        
    else:  # No transmittance
        if comp_polarised:
            return r_unp, r_prp,r_par
        else:
            return r_unp
            

# Obsolescent function aliases; wrapers around new functions for now, remove later.
def computeSunPos(geo_lon,geo_lat, year,month,day, hour,minute=0,second=0, timezone='UTC', debug=False):
    """Obsolete version of sun_position_from_date_and_time().  Use that function instead!"""
    print('ERROR: computeSunPos() is an obsolete version of sun_position_from_date_and_time().  Use that function instead!')
    exit(1)


def cosAngleSunPanels(sp_az,sp_incl, sun_az,sun_alt):
    """Obsolete version of cos_angle_sun_panels().  Use that function instead!"""
    print('ERROR: cosAngleSunPanels() is an obsolete version of cos_angle_sun_panels().  Use that function instead!')
    exit(1)


def extinctionFactor(airmass, return_value_below_horizon=False):
    """Obsolete version of extinction_factor().  Use that function instead!"""
    print('ERROR: extinctionFactor() is an obsolete version of extinction_factor().  Use that function instead!')
    exit(1)


def diffuse_radiation_projection_Perez87(doy, alt, surf_incl, theta, beam_norm,dif_horiz):
    """Obsolete version of diffuse_radiation_projection_perez87().  Use that function instead!"""
    print('ERROR: diffuse_radiation_projection_Perez87() is an obsolete version of diffuse_radiation_projection_perez87().  Use that function instead!')
    exit(1)


def diffuseRad_from_globalRad_sunshine(glob_horiz, sun_frac, sun_alt, i_ext=sol_const):
    """Obsolete version of diffuse_radiation_from_global_radiation_and_sunshine().  Use that function instead!"""
    print('ERROR: diffuseRad_from_globalRad_sunshine() is an obsolete version of diffuse_radiation_from_global_radiation_and_sunshine().  Use that function instead!')
    exit(1)
    

# Test code:
if __name__ == '__main__':
    print(cos_angle_sun_panels(0.0,40*d2r, 0.0,30*d2r))
