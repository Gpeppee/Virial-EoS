import numpy as np
import src.thermodynamics as thermodynamics
import src.beta_equilibrium as beta_equilibrium
from src.constants import *
import src.virial as virial
import src.numericals as numericals
import matplotlib.pyplot as plt
import src.phases as phases
import os
import src.initialize as initialize
from scipy.interpolate import CubicSpline
import pickle

t_f = 50
n_t = int(t_f/5)+1
temperature = np.linspace(0,t_f,n_t)
temperature = temperature[1:]#[7:8]
print('The temperatures are: ', temperature)


#########IF YOU CAN, VECTORIZE THIS
energy, phase_shift_total_nn = phases.READ_PHASE_SHIFTS('nn')
energy, phase_shift_total_np = phases.READ_PHASE_SHIFTS('np')



#initialize.save_neutron(0, 0, 0, 'b2', 0, True)
from joblib import Parallel, delayed

def process_temperature(i):
    a, b, c, d = virial.virial(temperature[i], energy, phase_shift_total_nn, phase_shift_total_np, False)
    b3         = (-0.5022*(a + 2**(-5/2)) + 3**(-5/2) )
    db3        = -0.5022*b 
    #a = 0
    for xrho in rhos:
        #Calculate and save thermal EoS data with 2nd order virial coefficient 
        Out = beta_equilibrium.Calculate_neutron(temperature[i], xrho, a, b, 0, 0)
        initialize.save_neutron(Out, temperature[i], xrho, 'b2', True, False)
        print('T=',temperature[i], 'rho=',round(xrho,6), '\tb2')
        #Calculate and save thermal EoS data with 2nd order virial coefficient + approx of 3rd order
        Out = beta_equilibrium.Calculate_neutron(temperature[i], xrho, a, b, b3, db3)
        initialize.save_neutron(Out, temperature[i], xrho, 'b3+', False, False)
        print('T=',temperature[i], 'rho=',round(xrho,6), '\tb3+')

        #Calculate and save thermal EoS data with 2nd order virial coefficient - approx of 3rd order
        Out = beta_equilibrium.Calculate_neutron(temperature[i], xrho, a, b, -b3, -db3)
        initialize.save_neutron(Out, temperature[i], xrho, 'b3-', False, False)
        print('T=',temperature[i], 'rho=',round(xrho,6), '\tb3-')
        
        
    
    #'''
# Assuming temperature, b2_n, db2_n, be_np, db2_np, and other variables are already defined

# Specify the number of parallel jobs (adjust as needed)
num_jobs =1


# Parallelize the outer loop using joblib
Parallel(n_jobs=num_jobs)(delayed(process_temperature)(i) for i in (np.arange(len(temperature))))
exit()













