import numpy as np
import src.beta_equilibrium as beta_equilibrium
from src.constants import *
import src.virial as virial
import src.phases as phases
import os
import src.initialize as initialize
from joblib import Parallel, delayed

# Ask the user if the user want to calculate beta-equilibrium
print('Do you want beta-equilibrium? (y/n)')
alp = input()

# Define temperature range and number of temperature points
t_f = 50
n_t = int(t_f/5)+1
temperature = np.linspace(0.1, t_f, n_t)
temperature = temperature[1:]  # Exclude the first temperature value
temperature = np.linspace(1,10, 10)
#temperature = np.append([1.8], temperature)
print('The temperatures are: ', temperature)

# If beta-equilibrium is not required
if alp == 'n':
    print('Do you want to span different proton fractions? (y/n)')
    alp2 = input()
    if alp2 == 'n':
        print('Input the proton fraction:')
        xxp = [float(input())]
    if alp2 == 'y':
        n_xp = 21
        xxp = np.linspace(0.01, 0.5, n_xp)
        print('The code will calculate ', str(int(n_xp)), ' points equally spaced between 0 and 0.5')

# Remove old data files if they exist
try:
    os.remove('./data//phase_shifts.dat')
    os.remove('./data//virial_coeff.dat')
    os.remove('./data//thermal_barion.dat')
    os.remove('./data//thermal_electron.dat')
    print('Old data files removed')
except:
    try:
        os.remove('./data//data.npy')
    except:
        try:
            os.remove('./data//data.h5')
        except:
            print('No previous data: start calculations')


# Read phase shifts for neutron-neutron and neutron-proton interactions
energy, phase_shift_total_nn = phases.READ_PHASE_SHIFTS('nn')
energy, phase_shift_total_np = phases.READ_PHASE_SHIFTS('np')

# Save phase shift data to file
m = open(filename_phase, 'a')
m.write('#Energy\t phase shift nn \t phase shift np\n')
for i in range(len(energy)):
    m.write('{}\t{}\t{}\n'.format(round(energy[i], 5), round(phase_shift_total_nn[i], 5), round(phase_shift_total_np[i], 5)))
m.close()

# Initialize file for virial coefficients
m2 = open(filename_virial, 'a')
m2.write('#T\t b2_n\t T*db2_n\t b2_np\t T*db2_np\t b2_nuc\t T*db2_nuc \n')
m2.close()

# Function to process at each temperature
def process_temperature(i):
    #Calculate virial coefficients and write to file
    a, b, c, d = virial.virial(temperature[i], energy, phase_shift_total_nn, phase_shift_total_np, print_virial=True, interaction=True)
    if alp == 'y':
        #beta-equilibrium is required
        for xrho in rhos:
          #For each density, calculate and save thermal EoS data with 2nd order virial coefficient
          try:
            Out, xxpp = beta_equilibrium.Bisection(temperature[i], xrho, a, b, c, d, xxp_in, xxp_f)
            initialize.save(npy_name, Out, temperature[i], xxpp, xrho, 'beta')
          except:
            print('No beta equilibrium found for T = ', temperature[i], 'rho = ', xrho)    
    elif alp == 'n':
        #beta-equilibrium is not required
        for xp in xxp:
            for xrho in rhos:
                #For each density and proton fraction, calculate and save thermal EoS data with 2nd order virial coefficient
                Out = beta_equilibrium.Calculate(temperature[i], xrho, a, b, c, d, xp)
                initialize.save(npy_name, Out, temperature[i], xp, xrho, 'y_p')

# Number of parallel jobs, from 1 to max number of treads available
num_jobs = 18

# Parallelize the processing of temperatures using joblib
Parallel(n_jobs=num_jobs)(delayed(process_temperature)(i) for i in (np.arange(len(temperature))))

exit()









