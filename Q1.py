"""
Curtis - Algorithm 3.3 & 3.4 implementation
"""

import math
import numpy as np
import utils as ul


def func(Xi, dt, r0, vr0, alpha, obj='sun'):
    
    """Universal Kepler's Equation"""
    
    """
    Inputs
    --------
    Xi - Universal Anomaly 
    dt - time interval (sec)
    r0 - initial distance (km)
    vr0 - initial radial velocity (km/s)
    alpha - reciprocal of semi-major axis (km^-1)
    
    Outputs
    --------
    RHS - LHS of universal kepler's equation
    
    """
    
    if obj=='earth':
        mu = 398600 # km^3s^-2
    elif obj == 'sun':
        mu = 1.327*(10**11)
        
    zi = alpha * Xi**2
    value = ((r0*vr0/np.sqrt(mu))*(Xi**2)*ul.stumpff_c(zi)) + \
            ((1 - alpha*r0)*(Xi**3)*ul.stumpff_s(zi)) + \
            r0*Xi - (np.sqrt(mu)*dt)
    
    return value

def d_func(Xi, dt, r0, vr0, alpha, obj='sun'):
    
    """Derivative of Universal Kepler's Equation wrt X"""
    
    """
    Inputs
    --------
    Xi - Universal Anomaly 
    dt - time interval (sec)
    r0 - initial distance (km)
    vr0 - initial radial velocity (km/s)
    alpha - reciprocal of semi-major axis (km^-1)
    
    Outputs
    --------
    d(RHS - LHS)/dX of universal kepler's equation
    
    """
    
    if obj=='earth':
        mu = 398600 # km^3s^-2
    elif obj == 'sun':
        mu = 1.327*(10**11)
        
    zi = alpha * Xi**2
    
    value = (r0*vr0/np.sqrt(mu))*Xi*(1 - alpha*(Xi**2)*ul.stumpff_s(zi)) + \
            ((1 - alpha*r0)*(Xi**2)*ul.stumpff_c(zi)) + r0
    
    return value

def f_func(dt, r0, alpha, X):
    
    """F function"""
    
    """
    Inputs:
    --------
    dt - time interval (sec)
    r0 - initial distance (km)
    alpha - reciprocal of semi-major axis (km^-1)
    X - Universal Anomaly
    """
    
    value = 1 - (X**2/r0)*ul.stumpff_c(alpha*X**2)
    
    return value

def g_func(dt, r0, alpha, X, obj='sun'):
    
    """G function"""
    
    """
    Inputs:
    --------
    dt - time interval (sec)
    r0 - initial distance (km)
    alpha - reciprocal of semi-major axis (km^-1)
    X - Universal Anomaly
    """
    
    if obj=='earth':
        mu = 398600 # km^3s^-2
    elif obj == 'sun':
        mu = 1.327*(10**11)
        
    value = dt - (1/np.sqrt(mu))*(X**3)*ul.stumpff_s(alpha*X**2)
    
    return value

def fdot_func(r0, r, alpha, X, obj='sun'):
    
    """Time derivative of F function"""
    
    """
    Inputs
    -------
    r0 - magnitude of initial position vector (km)
    r - magnitude of position vector (km)
    alpha - reciprocal of semi-major axis (km^-1)
    X - Universal Anomaly
    
    """
    
    if obj=='earth':
        mu = 398600 # km^3s^-2
    elif obj == 'sun':
        mu = 1.327*(10**11)
        
    value = (np.sqrt(mu)/(r0*r)) * (alpha*(X**3)*ul.stumpff_s(alpha*X**2) - X)
    
    return value

def gdot_func(r0, r, alpha, X):
    
    """Time derivative of G function"""
    
    """
    Inputs
    -------
    r0 - magnitude of initial position vector (km)
    r - magnitude of position vector (km)
    alpha - reciprocal of semi-major axis (km^-1)
    X - Universal Anomaly
    
    """
    
    value = 1 - ((X**2)/r)*ul.stumpff_c(alpha*X**2)
    
    return value

    
# Curtis Algorithm 3.3
    
def uni_kepler_eq(dt, r0, vr0, alpha, obj='sun'):
    
    """
    Inputs:
    dt - time interval (sec)
    r0 - initial distance (km)
    vr0 - initial radial velocity (km/s)
    alpha - reciprocal of semi-major axis (km^-1)
    
    Other quantities:
    mu = G(m1 + m2) - gravitational parameter
    
    Outputs:
    Xi - Universal anomaly
    
    """
    
    #print('Universal Kepler Equation inputs:')
    #print('dt = ', dt)
    #print('r0 = ', r0)
    #print('vr0 = ', vr0)
    #print('alpha = ', alpha)
    
    tol = 1e-8
    if obj=='earth':
        mu = 398600 # km^3s^-2
    elif obj == 'sun':
        mu = 1.327*(10**11)
    
    X0 = np.sqrt(mu) * np.abs(alpha) * dt
    
    print('Initial guess: ', X0)
    
    Xi = X0
    
    ratio_i = func(X0, dt, r0, vr0, alpha, obj) / d_func(X0, dt, r0, vr0, alpha, obj)
    
    p = 0
    
    while (np.abs(ratio_i) >= tol):
        
        Xip1 = Xi - ratio_i
        ratio_i = func(Xip1, dt, r0, vr0, alpha, obj) / d_func(Xip1, dt, r0, vr0, alpha, obj)
        Xi = Xip1
        p = p+1
        
    print('Solution converged in itetaions: ', p, '  ratio: ', ratio_i)
        
        
    return Xi
    

# Curtis Algorithm 3.4

def alg_3p4(r0, v0, dt, obj='sun'):
    
    #print('Propagator inputs:')
    #print('dt = ', dt)
    #print('r0 = ', r0)
    #print('v0 = ', v0)
    
    """
    Inputs:
    dt - time interval (sec)
    r0 - initial distance vector [i, j] components array (km)
    v0 - initial velocity vector [i, j] components array (km/s)
    
    Other quantities:
    mu = G(m1 + m2) - gravitational parameter
    alpha - 
    > 0 --> ellipse
    < 0 --> hyperbola
    = 0 --> parabola
    
    Outputs:
    r vector - [x, y] array
    v vector - [x, y] array
    """
    
    if obj=='earth':
        mu = 398600 # km^3s^-2
    elif obj == 'sun':
        mu = 1.327*(10**11)
    
    r0_mag = np.linalg.norm(r0)
    v0_mag = np.linalg.norm(v0)
    
    vr0 = np.dot(r0, v0)/r0_mag
    
    alpha = (2/r0_mag) - (v0_mag**2 / mu)
    
    # Universal Anomaly determination
    X = uni_kepler_eq(dt, r0_mag, vr0, alpha, obj)
    
    # f and g functions
    f_val = f_func(dt, r0_mag, alpha, X)
    g_val = g_func(dt, r0_mag, alpha, X, obj)
    
    r_vector = (f_val*r0) + (g_val*v0)
    r_mag = np.linalg.norm(r_vector)
    
    # fdot and gdot
    fdot = fdot_func(r0_mag, r_mag, alpha, X, obj)
    gdot = gdot_func(r0_mag, r_mag, alpha, X)
    
    v_vector = (fdot*r0) + (gdot*v0)
    
    return r_vector, v_vector
    
    

# Test problems

# ALg 3.3  - Curtis Example 3.6
# dt = 3600 # sec
# r0 = 10000 # km
# vr0 = 3.0752 # km/s
# alpha = -5.0878 * 1e-5 # km^-1

# X = uni_kepler_eq(dt, r0, vr0, alpha, obj='earth')
# print(X)

# Alg 3.4 - Curtis Example 3.7

# r0 = np.array([7000, -12124])
# v0 = np.array([2.6679, 4.6210])
# dt = 3600
# r,v = alg_3p4(r0, v0, dt, obj='earth')
# print(r)
# print(v)
