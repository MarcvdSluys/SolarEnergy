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


"""Functions for solar energy dealing with (solar) radiation.

References:
  - M. van der Sluys, Celestial mechanics in a nutshell, https://cmians.sourceforge.io (2014).
      
"""


import sys
import datetime as dt
import pytz as tz
import numpy as np

import soltrack as st


def computeSunPos(lon,lat, year,month,day, hour,minute=0,second=0, timezone='UTC', debug=False):
    """Compute the Sun local position (azimuth, altitude and distance) for the given geographical location and
    date and time using SolTrack.
    
    Args:
        lon (float):  Geographic longitude to compute the Sun position for (rad).
        lat (float):  Geographic latitude to compute the Sun position for (rad).
    
        year (int):    Year (CE) to compute the Sun position for.
        month (int):   Month to compute the Sun position for.
        day (int):     Day of month to compute the Sun position for.
    
        hour (int):    Hour of day to compute the Sun position for (local time!).
        minute (int):  Minute to compute the Sun position for (optional; default = 0).
        second (int):  Second to compute the Sun position for (optional; default = 0).
    
        timezone (timezone):  Time zone for which date and time are provided (optional; default = 'UTC').
    
        debug (bool):  Switch to write detailed output to stdout (optional; default = False).
    
    Returns:
        tuple (float,float,float):  Tuple containing (azimuth, altitude, distance):
    
            - azimuth (float):   Azimuth of the Sun (rad; south = 0 on the northern hemisphere).
            - altitude (float):  Altitude of the Sun (rad).
            - distance (float):  Distance Sun-Earth (AU).
    
    """
    
    myTZ = tz.timezone(timezone)   # My timezone
    myTime = dt.datetime(int(year),int(month),int(day), int(hour),int(minute),int(second))  # Time w/o timezone
    lt = myTZ.localize(myTime, is_dst=None)  # Mark as local time
    utc = lt.astimezone(tz.utc)              # Convert to UTC
    
    # Set up geographical location (in degrees, since useDegrees=True) in a SolTrack Location dataclass object:
    loc = st.Location(lon,lat)  # longitude (>0: east of Greenwich),  latitude (>0: northern hemisphere), in radians
    
    # Set (UT!) date and time in a SolTrack Time dataclass object:
    time = st.Time.datetime2st(utc)
    
    # Compute positions - returns a st.Position object:
    pos = st.computeSunPosition(loc, time, computeDistance=True)
    
    # Write results to standard output:
    if(debug):
        r2d = 180/np.pi  # Convert radians to degrees
        print("Location:  %0.3lf E, %0.3lf N"  % (loc.longitude*r2d, loc.latitude*r2d))
        print("Date:      %4d %2d %2d"         % (time.year, time.month, time.day))
        print("Time:      %2d %2d %9.6lf"      % (time.hour, time.minute, time.second))
        print("JD:        %0.11lf"             % (pos.julianDay))
        print()
        
        print("Corrected azimuth, altitude:  %10.6lf° %10.6lf°" % (pos.azimuthRefract*r2d, pos.altitudeRefract*r2d))
        print("Distance:                     %10.6lf AU"        % (pos.distance))
        print()
    
    return pos.azimuthRefract, pos.altitudeRefract, pos.distance



def cosAngleSunPanels(spAz,spIncl, sunAz,sunAlt):
    """Compute the cosine of the angle between the orientation vector of the solar panels and the position vector
       of the Sun.
    
    This is the cosine of the angle under which the direct sunlight hits the solar panels.  Multiply it with
    the DNI to obtain the direct insolation on the panels.  See Celestial mechanics in a nutshell, Sect. 4.3:
    Insolation on an inclined surface (http://CMiaNS.sf.net).
    
    Args:
        spAz (float):    Azimuth in which the solar panels are facing (rad; e.g. north or south = 0, same as sunAz).
        spIncl (float):  Inclination ('zenith angle') of the solar panels w.r.t. the horizontal (rad).
        sunAz (float):   Azimuth of the Sun (rad; e.g. north or south = 0, same as spAz).
        sunAlt (float):  Altitude of the Sun (rad).
    
    Returns:
        float:  The cosine between the normal vector of the solar panels and the position vector of the Sun (rad).
                Note that this value is zero (indicating radiation from behind the panels) or positive.

    """
    
    cosTheta = np.sin(sunAlt) * np.cos(spIncl)  +  np.cos(sunAlt) * np.sin(spIncl) * np.cos(sunAz - spAz)
    cosTheta = max(cosTheta, 0)  # Return 0 rather than negative values (indicating radiation from behind).
    
    return cosTheta



def airmass(sunAlt, returnValueBelowHorizon=False):
    """Compute airmass as a function of Sun altitude.
    
    Args:
        sunAlt (float):  Altitude of the Sun (rad).
        returnValueBelowHorizon (bool): Return a very large value when the Sun is below the horizon, larger
                                        when the Sun is lower.  This can be useful for solvers.  Default: False.
    
    Returns:
        float:  Airmass at sea level (AM~1 if the Sun is in the zenith, AM~38 near the horizon).
    
    """
    
    if(sunAlt < -0.00989):
        if(returnValueBelowHorizon):
            airmass = 1000 * (0.15 + abs(sunAlt))  # Very bad, but still getting worse for lower Sun, for solvers
        else:
            airmass = float('inf')
    else:
        airmass = (1.002432*np.sin(sunAlt)**2 + 0.148386*np.sin(sunAlt) + 0.0096467) / \
                  (np.sin(sunAlt)**2*np.sin(sunAlt) + 0.149864*np.sin(sunAlt)**2 + 0.0102963*np.sin(sunAlt) + 0.000303978)
        airmass = max( airmass, 1 )   # Air mass cannot be lower than 1
    
    return airmass



def extinctionFactor(airmass, returnValueBelowHorizon=False):
    """Compute the atmospheric extinction factor for sunlight from the air mass.
    
    Args:
        airmass (float):  Airmass at sea level (AM~1 if the Sun is in the zenith, AM~38 near the horizon).
        returnValueBelowHorizon (bool): Return a very large value when the Sun is below the horizon, larger
                                        when the Sun is lower.  This can be useful for solvers.  Default: False.
    
    Returns: 
        float: The extinciton factor for sunlight in the atmosphere.  Divide the extraterrestrial (AM)
               radiation (or, if unknown, solar constant) by this number to obtain the DNI.

    """
    
    if(airmass > 38.2):
        if(returnValueBelowHorizon):
            extFac = np.sqrt(sys.float_info.max) * (0.15 + airmass)  # Very bad, but still getting worse for higher airmass, for solvers
        else:
            extFac = float('inf')
    else:
        coefs = [ 9.1619283e-2, 2.6098406e-1,-3.6487512e-2, 6.4036283e-3,-8.1993861e-4, 6.9994043e-5,-3.8980993e-6,
                  1.3929599e-7, -3.0685834e-9, 3.7844273e-11,-1.9955057e-13]
        
        AMpow = 1.0                      # AM^0
        ext = coefs[0]                   # c_1 * AM^0
        for iCoef in range(1,len(coefs)):
            AMpow *= airmass             # AM^(i-1)
            ext += coefs[iCoef] * AMpow  # + c_i * AM^(i-1)
            
        extFac = np.exp(ext)
        
    return extFac


def diffuse_radiation_projection_Perez87(DoY, alt, surfIncl, theta, Gbeam_n,Gdif_hor):
    
    """Compute diffuse radiation on an inclined surface using the 1987 Perez model
    
    This function is adapted from the libTheSky Fortran implementation (libthesky.sf.net).
    
    See Perez et al. Solar Energy Vol. 39, Nr. 3, p. 221 (1987) - references to equations and tables are to
    this paper.  Most equations can be found in the Nomenclature section at the end of the paper (p.230).  I
    use a and c here, not b and d.
    
    
    Parameters:
      DoY       (int):     Day of year (Nday)
      alt       (float):   Altitude of the Sun (radians)
      
      surfIncl  (float):   Surface inclination wrt horizontal (radians) - 0 = horizontal, pi/2 = vertical
      theta     (float):   Angle between surface normal vector and Sun position vector (radians)
      
      Gbeam_n   (float):   Beam (direct) normal radiation (W/m2; in the direction of the Sun)
      Gdif_hor  (float):   Diffuse radiation on a horizontal surface (W/m2)
      
    Returns:
      float:    Diffuse irradiation on the inclined surface (W/m2) (output)

    """
    
    from solarenergy.constants import pi2,pio2, d2r,r2d
    
    # *** Compute the brightness coefficients for the circumsolar (F1) and horizon (F2) regions ***
    
    # 'External' (AM0) radiation:
    AM0rad = 1370 * (1 + 0.00333 * np.cos(pi2/365 * DoY))
    
    # Air mass:
    if(alt < -3.885*d2r):
        Mair = 99
    elif(alt < 10*d2r):
        Mair = 1 / ( np.sin(alt) + 0.15 * (alt*r2d + 3.885)**(-1.253) )
    else:
        Mair = 1 / np.sin(alt)
        
    Delta = Gdif_hor * Mair / AM0rad  # Brightness of overcast sky - par. 2.2.4 (a)
    
    
    # Cloudliness: epsilon;  epsilon ~ 1: overcast, epsilon -> infinity: clear  (epsilon ~ 1/fraction of covered sky)
    #   Needed for correct row in Table 1
    if(Gdif_hor <= 0):  # Division by zero
        if(Gbeam_n <= 0):  # No direct light: 0/0
            epsilon = 0    # -> completely overcast - first row of Table 1
        else:              # Some direct daylight: x/0 = large
            epsilon = 99   # -> completely clear, should be >11 for last row of Table 1
    else:
        epsilon = (Gdif_hor + Gbeam_n) / Gdif_hor  # Overcast: epsilon ~ 1,  clear: epsilon -> infinity
    
    
    
    # Table 1
    f11=0;  f12=1;  f13=2;  f21=3; f22=4; f23=5
    if(epsilon <= 1.056):
        F = [-0.011,  0.748, -0.080, -0.048,  0.073, -0.024]
    elif(epsilon <= 1.253):
        F = [-0.038,  1.115, -0.109, -0.023,  0.106, -0.037]
    elif(epsilon <= 1.586):
        F = [ 0.166,  0.909, -0.179,  0.062, -0.021, -0.050]
    elif(epsilon <= 2.134):
        F = [ 0.419,  0.646, -0.262,  0.140, -0.167, -0.042]
    elif(epsilon <= 3.230):
        F = [ 0.710,  0.025, -0.290,  0.243, -0.511, -0.004]
    elif(epsilon <= 5.980):
        F = [ 0.857, -0.370, -0.279,  0.267, -0.792,  0.076]
    elif(epsilon <= 10.080):
        F = [ 0.734, -0.073, -0.228,  0.231, -1.180,  0.199]
    else:
        F = [ 0.421, -0.661,  0.097,  0.119, -2.125,  0.446]
    
    
    zeta = pio2 - alt  # Zenith angle = pi/2 - alt
    F1 = F[f11]  +  F[f12] * Delta  +  F[f13] * zeta  # Circumsolar brightness coefficient
    F2 = F[f21]  +  F[f22] * Delta  +  F[f23] * zeta  # Horizon brightness coefficient
    
    
    
    
    # *** Compute the mean solid angles occupied by the circumsolar region (C and A, needed for Eq.8) ***
    
    alpha = 25*d2r  # Half angle of the circumsolar region (degrees -> radians; below Eq.9)
    
    
    # Solid angle of the circumsolar region weighted by incidence on the HORIZONTAL (variable C, subscript H;
    #   see Nomenclature, under c):
    # psiH:
    if(zeta > pio2 - alpha):
        psiH = 0.5 * (pio2 - zeta + alpha) / alpha  # Dimensionless ratio
    else:
        psiH = 1
    
    
    # chiH:
    if(zeta < pio2 - alpha):
        chiH = np.cos(zeta)  # = np.sin(alt)
    else:
        chiH = psiH * np.sin(psiH*alpha)
    
    
    C = 2 * (1 - np.cos(alpha)) * chiH  # Solid angle of the circumsolar region, weighted by HORIZONTAL incidence
    
    
    # Solid angle of the circumsolar region weighted by incidence on the SLOPE (variable A, subscript C;
    #   see Nomenclature, under c):
    # psiC:
    psiC = 0.5 * (pio2 - theta + alpha) / alpha
    
    # chiC:
    if(theta < pio2 - alpha):
        chiC = psiH * np.cos(theta)
    elif(theta < pio2 + alpha):
        chiC = psiH * psiC * np.sin(psiC*alpha)
    else:
        chiC = 0
    
    
    A = 2 * (1 - np.cos(alpha)) * chiC  # Solid angle of the circumsolar region, weighted by SLOPE incidence
    
    
    
    # Diffuse radiation from circumsolar (F1) and horizon (F2) regions on the inclined surface (Eq.8):
    Gdif_inc_csl = Gdif_hor * ( 0.5 * (1 + np.cos(surfIncl)) * (1 - F1)  +  F1 * A/C )  # Circumsolar
    Gdif_inc_hzl = Gdif_hor * ( F2 * np.sin(surfIncl) )                                         # Horizon band
    
    Gdif_inc_csl = max(Gdif_inc_csl, 0)  # Components are sometimes negative
    Gdif_inc_hzl = max(Gdif_inc_hzl, 0)
    Gdif_inc = Gdif_inc_csl + Gdif_inc_hzl
    
    # Assign optional return values:
    # if(present(Gdif_inc_cs)) Gdif_inc_cs = Gdif_inc_csl
    # if(present(Gdif_inc_hz)) Gdif_inc_hz = Gdif_inc_hzl
    
    return Gdif_inc



def clearsky_bird(alt, Iext=1353,Rsun=1, Press=1013,  Uo=0.34,Uw=1.42, Ta5=0.2661,Ta3=0.3538,Ba=0.84,K1=0.1, Rg=0.2):
    
    """A simplified clear-sky model for direct and diffuse insolation on horizontal surfaces.
    
    A.k.a. as "the Bird model".
    
    This function is adapted from the libTheSky Fortran implementation (libthesky.sf.net).
    
    See Bird & Hulstrom, A simplified clear-sky model for direct and diffuse insolation on horizontal
    surfaces, SERI/TR-642-761 (1981).
    
    Note that the value of Taa does not agree with tabulated values from the paper, and hence neither do
    dependent values (except for AM~1).  When I substitute their values for Taa, everything matches perfectly.
    Error in their formula, or (hopefully!) in their table?

    Parameters:
      alt Sun altitude above the horizon (rad)
    
      Iext   Extraterrestrial radiation (at the top of the atmosphere; AM0; W/m^2 - optional, default: 1353 (1361.5))
      Rsun   Sun distance (AU - optional, default: 1)
      Press  Air pressure at the observer's site, corrected for altitude (hPa - optional, default: 1013)
    
      Uo     Ozone abundance in a vertical column (cm - optional, default: 0.34)
      Uw     Percipitable water-vapor abundance in a vertical column (cm - optional, default: 1.42)
    
      Ta5    Aerosol optical depth from surface in vertical path at 500 nm (optional, default: 0.2661)
      Ta3    Aerosol optical depth from surface in vertical path at 380 nm (optional, default: 0.3538)
      Ba     Aerosol forward-scattering ratio  (optional, 0.82-0.86, default: 0.84)
      K1     Aerosol-absorptance constant (optional, rural: 0.0933, urban: 0.385, default: 0.1)
    
      Rg     Ground albedo (optional, fraction - default: 0.2)
    
    
    Returns:
      tuple (float,float,float,float):  Tuple containing (rv1, rv2):
        
      - Itot (float):  Total insolation on a horizontal surface (W/m^2)
      - Idir (float):  Direct (beam) insolation on a horizontal surface (W/m^2)
      - Idif (float):  Diffuse insolation on a horizontal surface (W/m^2)
      - Igr  (float):  Ground-reflection insolation from a horizontal surface (W/m^2)

    """
    
    from solarenergy.constants import pio2, r2d
    
    Z = pio2 - alt  # Solar zenith angle
    cosZ = np.cos(Z)   # Save a few CPU cycles
    
    
    # Relative air mass for the solar vector:
    AM  = 1/(cosZ + 0.15 * (93.885-Z*r2d)**(-1.25))  # Air mass
    AMp = AM * Press / 1013                          # Pressure-corrected air mass
    
    
    # TRANSMISSION EQUATIONS:
    # Rayleigh scattering:
    Tr = np.exp( -0.0903 * AMp**0.84 * (1 + AMp - AMp**1.01) )
    
    # Ozone:
    Xo = Uo*AM  # Amount of ozone in the direction of the Sun
    To = 1  -  0.1611 * Xo * (1+139.48*Xo)**(-0.3035)  -  0.002715 * Xo / (1 + 0.044*Xo + 0.0003*Xo**2)  # Transmittance of ozone absorptance
    
    # Uniformly mixed gases (CO2, O2):
    Tum = np.exp(-0.0127 * AMp**0.26)  # Transmittance of mixed-gas absorptance
    
    # Water vapor:
    Xw = AM*Uw  # Amount of water vapor in the direction of the Sun
    Tw = 1 - 2.4959 * Xw / ((1 + 79.034*Xw)**0.6828 + 6.385*Xw)             # Transmittance of water-vapor absorptance - Tw = 1-Aw
    
    # Daily turbidity:
    Tau = 0.2758*Ta3 + 0.35*Ta5                                             # Broadband turbidity: aerosol optical depth from surface in a vertical column
    Ta  = np.exp( -Tau**0.873  *  (1 + Tau - Tau**0.7088)  *  AM**0.9108 )  # Transmittance of aerosol absorptance and scattering
    Taa = 1 - K1 * (1 - AM + AM**1.06) * (1-Ta)                             # Transmittance of aerosol absorptance - this does not agree with tabulated values from the paper (except for AM~1).  When I substitute their values for Taa, everything matches perfectly.  Error in their formula, or in their table?
    Tas = Ta/Taa                                                            # Transmittance of aerosol scattering
    Rs  = 0.0685 + (1-Ba) * (1-Tas)                                         # Sky albedo
    
    
    # IRRADIANCE EQUATIONS:
    # Direct radiation on a horizontal surface:
    tmpVar = Iext * cosZ  *  To * Tum * Tw  # Save a few CPU cycles
    Idir = 0.9662 * tmpVar  *  Tr * Ta  /  Rsun**2
    
    # Diffuse (scattered) radiation on a horizontal surface:
    Idif  = 0.79 *  tmpVar        * Taa *  (0.5*(1-Tr) + Ba*(1-Tas)) / (1 - AM + AM**1.02)
    
    # Total (direct+diffuse) radiation on a horizontal surface:
    Itot = (Idir+Idif) / (1 - Rg*Rs)
    
    # Ground-reflected radiation from a horizontal surface:
    Igr  = Itot - (Idir+Idif)
    
    return Itot, Idir, Idif, Igr


def diffuseRad_from_globalRad_sunshine(Gglob_hor, sunFrac, sunAlt, Iext):
    """Compute the diffuse horizontal radiation from the global horizontal radiation and the Sun altitude.
    
    Parameters:
      Gglob_hor (float):  Global horizontal radiation (W/m2).
      sunFrac   (float):  Fraction of sunshine (e.g. fraction of cloud cover) (-; 0-1).
      sunAlt    (float):  Sun altitude above the horizon (rad).
    
    Returns:
      tuple (float,float,float):  Tuple containing (Gdif_hor, Gbeam_hor, DNI):
      
        - Gdif_hor  (float):  Diffuse horizontal radiation (W/m2).
        - Gbeam_hor (float):  Beam (direct) horizontal radiation (W/m2).
        - DNI       (float):  DNI = direct (beam) normal irradiation (W/m2).
    """
    
    DNI = Iext / extinctionFactor(airmass(sunAlt)) * sunFrac  # (Mean) DNI
    Gbeam_hor = DNI*np.sin(sunAlt)
    Gdif_hor  = Gglob_hor - Gbeam_hor
    
    return Gdif_hor, Gbeam_hor, DNI
