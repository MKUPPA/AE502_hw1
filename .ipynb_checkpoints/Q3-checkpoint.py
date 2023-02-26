"""
Contour plots 
"""
import sys
import math
import Q1 as q1
import Q2 as q2
import numpy as np
import utils as ul
import pandas as pd
import matplotlib.pyplot as plt
from datetime import date, timedelta


"""
Obj1 - 1I/’Oumouamoua
Obj2 - 2I/Borisov
"""

def get_dv(date1, date2, r1I, v1I, rE, vE, mission='rendezvous'):
    
    """
    Inputs
    ----------
    date1 - Departure date
    date2 - Arrival date
    r1I - Initial state vector of interstellar object
    v1I - Initial velocity vector of interstellar object
    rE - Initial state vector of satellite
    """
    
    dt = ul.numOfsec(date1, date2)
    
    
    """Interstellar Object"""
    
    print('Interstellar Object Calculations')
    
    rI_final, vI_final = q1.alg_3p4(r0=r1I, v0=v1I, dt=dt)

    """Satellite"""
    
    print('Satellite Calculations')
    
    vsat_init, vsat_final, p, ratio = q2.lambert(r1=rE, r2=rI_final, dt=dt)
    
    
    if mission == 'rendezvous':
        
        dv = np.linalg.norm(vsat_init - vE) + np.linalg.norm(vsat_final - vI_final)
        
    elif mission == "flyby":
        
        dv = np.linalg.norm(vsat_init - vE) 
    
    return dv, p, ratio

"""1I/’Oumouamoua"""

r1I_O = np.array([3.5158e-2, -3.16204, 4.49398])
v1I_O = np.array([-2.3175e-3, 9.8433e-3, -1.5418e-2])

"""2I/Borisov"""

r1I_B = np.array([7.2494, 14.6106, 14.2427]) 
v1I_B = np.array([-8.2417e-3,-1.15621e-2,-1.31713e-2]) 

"""Earth"""
rE = np.array([-1.7961e-1, 9.6679e-1,-3.6686e-5])
vE = np.array([-1.7200e-2,-3.2111e-3, 7.9277e-7])

rE = ul.au_to_km(rE)
vE = ul.au_per_day_to_km_per_s(vE)              

def porkchop(interstellar_obj, mission):
    
    dep_date_list = []
    arr_date_list = []
    dv_arr = []

    if interstellar_obj == 'Oumouamoua':

        r1I = ul.au_to_km(r1I_O)
        v1I = ul.au_per_day_to_km_per_s(v1I_O)
        
        sdate = date(2017,1,1)
        edate = date(2017,12,1)

        dep_dates = ul.get_dates(sdate, edate)

        sdate = date(2017,8,1)
        edate = date(2019,1,1)

        arr_dates = ul.get_dates(sdate, edate)

    elif interstellar_obj == 'Borisov':

        r1I = ul.au_to_km(r1I_B)
        v1I = ul.au_per_day_to_km_per_s(v1I_B)
        
        sdate = date(2017,1,1)
        edate = date(2020,7,1)

        dep_dates = ul.get_dates(sdate, edate)

        sdate = date(2019,6,1)
        edate = date(2022,1,1)

        arr_dates = ul.get_dates(sdate, edate)
    
    
    for num1, iter1 in enumerate(dep_dates):

        dep_date = iter1 
        
        date_indices = np.where((arr_dates > dep_date))[0]

        for num2, iter2 in enumerate(arr_dates[date_indices]):

            arr_date = iter2 

            if num1 != 0:

                dt = ul.numOfsec(date1=date(2017, 1, 1), date2=dep_date)
                rE_sat, vE_sat = q1.alg_3p4(r0=rE, v0=vE, dt=dt)
                r1I_new, v1I_new = q1.alg_3p4(r0=r1I, v0=v1I, dt=dt)

            else:

                rE_sat=rE
                vE_sat=vE
                r1I_new = r1I
                v1I_new = v1I



            dv, p, ratio = get_dv(dep_date, arr_date, r1I_new, v1I_new, rE_sat, vE_sat, mission)


            if p > 100 or math.isnan(ratio):
                continue
            else:
                dv_arr.append(dv)
                dep_date_list.append(dep_date)
                arr_date_list.append(arr_date)
                
            print(f'departure date: {dep_date}, arrival date: {arr_date}, dv: {dv}')
                
                
    dv_arr = np.array(dv_arr)
    dep_date_list = np.array(dep_date_list)
    arr_date_list = np.array(arr_date_list)

    if mission == 'rendezvous' and interstellar_obj == 'Oumouamoua':
        indices = np.where((dv_arr<50))[0]
    
    elif mission == 'rendezvous' and interstellar_obj == 'Borisov':
        indices = np.where((dv_arr<60))[0]
        
    elif mission == 'flyby':
        indices = np.where((dv_arr<20))[0]

    plt.figure()
    plt.scatter(dep_date_list[indices], arr_date_list[indices], c=dv_arr[indices], cmap='plasma')
    cb = plt.colorbar()
    cb.set_label('DV')
    plt.xlabel('Departure Date')
    plt.ylabel('Arrival Date')
    plt.xticks(rotation = 30) 
    plt.savefig('plots/' + str(interstellar_obj) + '_' + str(mission) + \
                '.png',dpi=600,bbox_inches='tight',pad_inches=0)
    plt.show()  
    
    
    return 
    
    
mission = 'rendezvous'
interstellar_obj = 'Borisov'  
                
porkchop(interstellar_obj, mission)
        
                
                
            
            