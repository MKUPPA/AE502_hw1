import math
import numpy as np
import pandas as pd
from datetime import date, timedelta

np.set_printoptions(precision=5, suppress=True, floatmode='fixed')

# compute all dates between start and end date

def get_dates(sdate, edate):
    
    a = pd.date_range(sdate,edate-timedelta(days=1),freq='3d').strftime('%Y,%m,%d').tolist()

    dates = []
    for num1, iter1 in enumerate(a):
        dates.append(date(int(iter1.split(',')[0]), int(iter1.split(',')[1]), int(iter1.split(',')[2])))
    
    return np.array(dates)


# compute dt gives date 1 and date 2

def numOfsec(date1, date2):
    
    """
    Input
    ----------
    date format - date(year, month, day)
    
    Output
    ----------
    dt in seconds
    """
    
    ndays =  (date2-date1).days
    nsec = ndays*(24*60*60)
    
    return nsec

# au to km

def au_to_km(x):
    
    """
    Inputs
    ------
    x in au
    
    """
    
    return x*(1.496e+8)

# au/day to km/s

def au_per_day_to_km_per_s(x):
    
    """
    Inputs
    ------
    x in au/day
    
    """
    
    return x*1731.46

# Stumpff C function

def stumpff_c(z):
    
    """
    z - in radians 
    """
    
    if z > 0:
        
        value = (1 - np.cos(np.sqrt(z)))/z
        
    elif z < 0:
        
        value = (np.cosh(np.sqrt(-z)) - 1)/(-z)
        
    else:
        
        value = 0.5
        
    return value

# Stumpff S function

def stumpff_s(z):
    
    if z > 0:
        
        value = (np.sqrt(z) - np.sin(np.sqrt(z)))/(np.sqrt(z))**3
        
    elif z < 0:
        
        value = (np.sinh(np.sqrt(-z)) - np.sqrt(-z))/(np.sqrt(-z))**3
        
    else:
        
        value = 1/6
        
    return value

# F function

def f_func(z, r1_mag, r2_mag, A):
    
    """
    Inputs:
    

    """
    C = stumpff_c(z)
    S = stumpff_s(z)
    y = r1_mag + r2_mag + (A*(z*S - 1)/(np.sqrt(C)))
    
    value = 1 - (y/r1_mag) 
    
    return value

# G function

def g_func(z, r1_mag, r2_mag, A, obj='sun'):
    
    """
    Inputs:
    dt - time interval (sec)

    """

    C = stumpff_c(z)
    S = stumpff_s(z)
    y = r1_mag + r2_mag + (A*(z*S - 1)/(np.sqrt(C)))
    
    if obj=='earth':
        mu = 398600 # km^3s^-2
    elif obj == 'sun':
        mu = 1.327*(10**11)
        
    value = A*np.sqrt(y/mu)
    
    return value

# G dot

def gdot_func(z, r1_mag, r2_mag, A):
    
    C = stumpff_c(z)
    S = stumpff_s(z)
    y = r1_mag + r2_mag + (A*(z*S - 1)/(np.sqrt(C)))
    
    value = 1 - (y/r2_mag) 
    
    return value


# Curtis Algorithm 4.1

def orbital_elem(r,v,obj='sun'):
    
    
    """
    Inputs
    ----------
    r: Initial position vector in km
    v: Initial velocity vector in km/s
    
    """
    
    if obj=='earth':
        mu = 398600 # km^3s^-2
    elif obj == 'sun':
        mu = 1.327*(10**11)
    
    print('\n')
    print(f'r = {repr(r)} (km)\n')
    print(f'v = {repr(v)} (km/s)\n')
    
    
    r_mag = np.linalg.norm(r)
    v_mag = np.linalg.norm(v)
    
    # radial velocity
    # v_rdl > 0 --> sat flying away from perigee
    # v_rdl < 0 --> sat flying toward perigee

    v_rdl = np.dot(r,v)/r_mag
    print(f'Radial velocity = {v_rdl:.5e} (km/s)\n')
    
    # specific angular momentum
    h = np.cross(r, v)
    
    print(f'Sp. Angular Momentum = {repr(h)} (km^2/s)\n')
    
    """First Orbital Element"""
    h_mag = np.linalg.norm(h)
    
    print(f'Sp. Angular Momentum Magnitude = {h_mag:5e} (km^2/s)\n')
        
    """Second Orbital Element"""
    # inclination in radians
    i = np.arccos(h[-1]/h_mag)
    
    print(f'Inclination in deg = {i*180/np.pi:5f}\n')
    
    # Node line
    N = np.cross(np.array([0,0,1]), h)
    print(f'Node line = {repr(N)}\n')
    
    N_mag = np.linalg.norm(N)
    print(f'Node line Mag = {N_mag}\n')
    
    """Third Orbital Element"""
    # Right Ascension of ascending node
    if N[1] >= 0:
        RA = np.arccos(N[0]/N_mag)
        
    else:
        RA = 2*np.pi - np.arccos(N[0]/N_mag)
        
    print(f'(Right Ascension of ascending node) RA in deg = {RA*180/np.pi:.5f}\n')
    
    # Eccentricity Vector
    e = (1/mu)*((v_mag**2 - (mu/r_mag))*r - r_mag*v_rdl*v)
    print(f'Eccentricity vector = {repr(e)}\n')
    
    """Fourth Orbital Element"""
    # Eccentricity
    e_mag = np.linalg.norm(e)
    print(f'Eccentricity Magnitude = {e_mag}\n')
    
    if e_mag < 1:
        
        print('==> Elliptical Orbit\n')
        a = (h_mag**2/mu)/(1 - e_mag**2)
        
    if e_mag > 1:
        
        print('==> Hyperbolic Orbit\n')
        a = (h_mag**2/mu)/(e_mag**2 - 1)

    
    """Fifth Orbital Element"""
    # Argument of perigee
    if e[-1] >= 0:
        omega = np.arccos(np.dot(N, e)/(N_mag*e_mag))
    else:
        omega = 2*np.pi - np.arccos(np.dot(N, e)/(N_mag*e_mag))
        
    print(f'(Argument of perigee) Omega in deg = {omega*180/np.pi:.5f}\n')
        
    """Sixth Orbital Element"""
    # True anomaly
    
    if v_rdl >=0:
        theta = np.arccos(np.dot(e, r)/(e_mag*r_mag))
    else:
        theta = 2*np.pi - np.arccos(np.dot(e, r)/(e_mag*r_mag))
        
    print(f'(True Anomaly) Theta in deg = {theta*180/np.pi:.5f}\n')
    
    print(f'Semimajor axis = {a:.5f} (km)\n')
    
    # Orbit period
    T = (2*np.pi*np.power(a, 1.5))/np.sqrt(mu)
    
    print(f'Orbit period = {T/(3600):.5f} (hr)\n')
    
    return [h_mag, i, RA, e_mag, omega, theta]


# Testing Alg 4.1

"""Curtis Example 4.3"""

# r = np.array([-6045, -3490, 2500])
# v = np.array([-3.457, 6.618, 2.533])

# res = orbital_elem(r,v, obj='earth')
# print(res)

"""1I/’Oumouamoua - uncomment below lines to run"""

# x1 = np.array([3.5158e-2, -3.16204, 4.49398])
# x2 = np.array([-2.3175e-3, 9.8433e-3, -1.5418e-2])

# r = au_to_km(x1)
# v = au_per_day_to_km_per_s(x2)

# res = orbital_elem(r,v,obj='sun')
# print(res)


"""2I/Borisov - uncomment below lines to run"""

# x1 = np.array([7.24947, 14.61063, 14.24274])
# x2 = np.array([-8.241709e-3, -1.156219e-2, -1.3171359e-2])

# r = au_to_km(x1)
# v = au_per_day_to_km_per_s(x2)

# res = orbital_elem(r,v,obj='sun')
# print(res)

    
    
    
    
    
    
    

