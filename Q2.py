"""
Lambert solver - curtis 5.2
"""

import numpy as np
import utils as ul
import matplotlib.pyplot as plt

# see where F changes sign

def F_val(z, r1_mag, r2_mag, A, dt, obj='sun'):
    
    if obj=='earth':
        mu = 398600 # km^3s^-2
    elif obj == 'sun':
        mu = 1.327*(10**11)
    
    C = ul.stumpff_c(z)
    S = ul.stumpff_s(z)
    y = r1_mag + r2_mag + (A*(z*S - 1)/(np.sqrt(C)))
    
    F = ((y/C)**(3/2))*S + (A*np.sqrt(y)) - (np.sqrt(mu)*dt)
    
    return F

def Fdot_val(z, r1_mag, r2_mag, A, dt, obj='sun'):
    
    if obj=='earth':
        mu = 398600 # km^3s^-2
    elif obj == 'sun':
        mu = 1.327*(10**11)
    
    C = ul.stumpff_c(z)
    S = ul.stumpff_s(z)
    y = r1_mag + r2_mag + (A*(z*S - 1)/(np.sqrt(C)))
    
    if z != 0 :
    
        Fdot = ((y/C)**(3/2))*((0.5/z)*(C - (1.5*S/C)) + (0.75*(S**2)/C)) + \
                (A/8)*((3*S*np.sqrt(y)/C) + (A*np.sqrt(C/y))) 
        
    else:
        
        Fdot = ((np.sqrt(2)/40)*(y**(1.5))) + (A/8)*(np.sqrt(y) + A*np.sqrt(1/(2*y))) 
            
    
    return Fdot
    
    
# Curtis Alg 5.2

def lambert(r1, r2, dt, obj='sun'):
    
    """
    Inputs
    ----------
    dt - time interval (sec)
    r1 - position vector of point 1 [i, j] components array (km)
    r2 - position vector of point 2 [i, j] components array (km)
    """
    
    #print('Inputs to lambert solver:')
    #print('r1 = ', r1)
    #print('r2 = ', r2)
    #print('dt = ', dt)

    
    if obj=='earth':
        mu = 398600 # km^3s^-2
    elif obj == 'sun':
        mu = 1.327*(10**11)
        
    #print('mu = ', mu)
        
    tol = 1e-8
    
    r1_mag = np.linalg.norm(r1)
    r2_mag = np.linalg.norm(r2)
    
    #print('r1 Mag = ', r1_mag)
    #print('r2 Mag = ', r2_mag)

    
    # Suppose Prograde trajectory --> cos(i) > 0
    z_cmp = np.cross(r1, r2)[-1]
    
    # Change in true anomaly in radians
    if z_cmp >= 0:
        
        delta_theta = np.arccos(np.dot(r1, r2)/(r1_mag*r2_mag))
        
    else:
        
        delta_theta = (2*np.pi) - np.arccos(np.dot(r1, r2)/(r1_mag*r2_mag))
        
    #print('Change in True Anomaly in deg = ', delta_theta*180/np.pi)

    A = np.sin(delta_theta)*np.sqrt((r1_mag*r2_mag)/(1 - np.cos(delta_theta)))
    
    #print('A = ', A)
    
    # solving for initial value of z
    z_guess = np.linspace(0,10,101)
    val1 = F_val(0, r1_mag, r2_mag, A, dt, obj)

    for num, val in enumerate(z_guess):
        val2 = F_val(val, r1_mag, r2_mag, A, dt, obj)
        
        if val1*val2 < 0:
            break
            
        else:
            val1 = val2
        
    
    slope = (F_val(z_guess[num],r1_mag, r2_mag, A, dt, obj) - \
             F_val(z_guess[num-1],r1_mag, r2_mag, A, dt, obj))/\
    (z_guess[num] - z_guess[num-1])
    zcorr = z_guess[num] - (F_val(z_guess[num],r1_mag, r2_mag, A, dt, obj)/slope)
    
    # Solve for z
    z0 = zcorr #0 #1e-8
       
    #print('Initial guess: ', z0)
    
    zi = z0
    
    ratio_i = F_val(zi, r1_mag, r2_mag, A, dt, obj) / Fdot_val(zi, r1_mag, r2_mag, A, dt, obj)
    
    p = 0
    
    while (np.abs(ratio_i) >= tol):
        
        zip1 = zi - ratio_i
        ratio_i = F_val(zip1, r1_mag, r2_mag, A, dt, obj) / Fdot_val(zip1, r1_mag, r2_mag, A, dt, obj)
        zi = zip1
        p = p+1
        
        if p > 100:
            break
        
    print(f'Solution converged in itetaions: {p}, ratio: {ratio_i}, solution: {zi}')
    
    z = zi
    
    if z > 0:
        
        print('Elliptical orbit')
        
    elif z < 0:
        
        print('Hyperbolic orbit')
        
    else:
        
        print('Parabolic orbit')
        
    C = ul.stumpff_c(z)
    S = ul.stumpff_s(z)
    y = r1_mag + r2_mag + (A*(z*S - 1)/(np.sqrt(C)))
    
    # f, g and gdot
    
    f_val = ul.f_func(z, r1_mag, r2_mag, A)
    g_val = ul.g_func(z, r1_mag, r2_mag, A, obj)
    gdot_val = ul.gdot_func(z, r1_mag, r2_mag, A)
    
    # calculate v1 and v2 vectors
    
    v1_vec = (1/g_val)*(r2 - f_val*r1)
    v2_vec = (1/g_val)*(gdot_val*r2 - r1)
    
    # Calculate Orbital Elements
    #ul.orbital_elem(r1,v1_vec,obj)
    
    return v1_vec, v2_vec, p, ratio_i

# # Testing Lambert Solver - Example 5.2
# r1 = np.array([5000, 10000, 2100])
# r2 = np.array([-14600, 2500, 7000])
# dt = 60*60

# v1_vec, v2_vec, p, ratio_i = lambert(r1, r2, dt, obj ='earth')

# print(v1_vec)
# print(v2_vec)

    
    
    

    