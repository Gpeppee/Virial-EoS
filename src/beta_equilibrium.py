import numpy as np
import src.thermodynamics as thermodynamics
from src.constants import *
from time import sleep
np.seterr( invalid='ignore')

def  Bisection(T,rho,b2_n,db2_n,b2_np,db2_np,a,b):   
     """
    Bisection method to find the chemical potentials and proton fraction.
    Parameters:
    T (float): Temperature.
    rho (float): Density.
    b2_n (float): Second virial coefficient for neutrons.
    db2_n (float): Derivative of the second virial coefficient for neutrons.
    b2_np (float): Second virial coefficient for neutron-proton interactions.
    db2_np (float): Derivative of the second virial coefficient for neutron-proton interactions.
    a (float): Lower bound for bisection.
    b (float): Upper bound for bisection.
    Returns:
    dict: Dictionary containing the final values of the bisection for 'b' and 'e'.
    float: Final value of the proton fraction.
     """
     xx0 = a			#first value
     xx1 = b			#second value
     
     print('Rho entering bisection is ', rho)
     
     #First evaluation of the chemical potentials
     xmu_n0, xmu_p0         = thermodynamics.thermo(T,rho,xx0,b2_n,db2_n,b2_np,db2_np,False,False)   #
     xmu_e0                 = thermodynamics.thermo_e(T,rho,xx0,False) #
     #Second evaluation of the chemical potentials
     xmu_n1, xmu_p1         = thermodynamics.thermo(T,rho,xx1,b2_n,db2_n,b2_np,db2_np,False,False)   #
     xmu_e1                 = thermodynamics.thermo_e(T,rho,xx1,False) # 
     #Beta equilibrium condition evaluation
     func0      = beta_equilibrium(xmu_n0, xmu_p0, xmu_e0)   
     func1      = beta_equilibrium(xmu_n1, xmu_p1, xmu_e1)
 
     #Check if the initial conditions is correct in order to start the bisection
     if func1*func0>0:
        print('Error in choosing initial data')
        print('Error in T=', T, 'rho=', rho)
        print('Chem pot are xmu_e1 = ', xmu_e1, 'xmu_n1 = ', xmu_n1, 'xmu_p1 = ', xmu_p1, 'xmu_e0 = ', xmu_e0, 'xmu_n0 = ', xmu_n0, 'xmu_p0 = ', xmu_p0)
        print('func1 = ', func1, 'func0 = ', func0)
        return 0
     
     #Bisection loop
     #Initialize the counter
     k=0
     func2 = 10
     while np.abs(func2)>1e-7:
       #Evaluation at the middle point
       xx2 = (xx1+xx0)/2
       xmu_n2, xmu_p2         = thermodynamics.thermo(T,rho,xx2,b2_n,db2_n,b2_np,db2_np,False,False)   #
       xmu_e2                 = thermodynamics.thermo_e(T,rho,xx2,False) 
       func2 = beta_equilibrium(xmu_n2, xmu_p2, xmu_e2)
       
       #Check the sign third evaluation
       prod = func0*func2

       #Update the values
       if prod<0:               #Left side
           xx1=xx2
           func0 = func0
           k+=1
       elif prod>0:             #Right side
           xx0=xx2
           func0 = func2
           k+=1
       elif prod==0:            #Exact value
           print('Exact value found for rho=', rho, 'T=', T, 'xp=', xx2, 'iter=', k)
       if k >1000:              #Too many iterations
           print('Too many iterationw in T=', T, 'rho=', rho)
           print('Chem pot are xmu_e1 = ', xmu_e1, 'xmu_n1 = ', xmu_n1, 'xmu_p1 = ', xmu_p1, 'xmu_e0 = ', xmu_e0, 'xmu_n0 = ', xmu_n0, 'xmu_p0 = ', xmu_p0)
           print('func0 = ', func0)
           return 0
       xxp_final = xx2
     
     #Final evaluation of the chemical potentials with output of the thermodynamic properties
     out = {}
     xmu_n2, xmu_p2, out_b         = thermodynamics.thermo(T,rho,xxp_final,b2_n,db2_n,b2_np,db2_np,True,True)   #
     out['b']=out_b
     
     
     xmu_e2, out_e                 = thermodynamics.thermo_e(T,rho,xxp_final,True) #T,rho,xxp_final,xmu_e2,.True.
     out['e']=out_e
     
     #Print the initial values in order to take track of the overall calculation
     print('Final values of the bisection are: rho=', rho, 'T=', T, 'xp=', xxp_final, 'iter=', k, 'xmu_e = ', xmu_e2, 'xmu_n = ', xmu_n2, 'xmu_p = ', xmu_p2)
     return out, xxp_final

     
def  beta_equilibrium(xmu_n, xmu_p, xmu_e):   
     """
    Calculate the beta equilibrium condition.
    Parameters:
    xmu_n (float): Chemical potential of neutrons.
    xmu_p (float): Chemical potential of protons.
    xmu_e (float): Chemical potential of electrons.
    Returns:
    float: Beta equilibrium value.
     """    
     beta = xmu_n - xmu_p - xmu_e + xmassn - xmassp - xmasse                                 
     return beta
     
def Calculate(T,rho,b2_n,db2_n,b2_np,db2_np,xxp):
     """
    Calculate the chemical potentials and other thermodynamic properties.
    Parameters:
    T (float): Temperature.
    rho (float): Density.
    b2_n (float): Second virial coefficient for neutrons.
    db2_n (float): Derivative of the second virial coefficient for neutrons.
    b2_np (float): Second virial coefficient for neutron-proton interactions.
    db2_np (float): Derivative of the second virial coefficient for neutron-proton interactions.
    xxp (float): Proton fraction.
    Returns:
    dict: Dictionary containing the calculated thermodynamic properties for 'b' and 'e'.
     """
     out_b = {}
     xmu_n2, xmu_p2, out_b         = thermodynamics.thermo(T,rho,xxp,b2_n,db2_n,b2_np,db2_np,True,True)   #
     xmu_e2, out_e                 = thermodynamics.thermo_e(T,rho,xxp,True)
     TT = str(int(T))
     XXPP = str(int(xxp*100))
     out = {}
     out['b']=out_b
     out['e']=out_e
     print('T=',int(T),'\trho=',round(rho,5),'\txxp=',round(xxp,4),'\txmu_n=', round(xmu_n2,1),'\txmu_p=', round(xmu_p2,1), '\txmu_e=', round(xmu_e2,4))
     return out

def Calculate_neutron(T,rho,b2_n,db2_n,b3_n,db3_n):
     """
     Calculate the chemical potentials and other thermodynamic properties for neutrons.
     Parameters:
     T (float): Temperature.
     rho (float): Density.
     b2_n (float): Second virial coefficient for neutrons.
     db2_n (float): Derivative of the second virial coefficient for neutrons.
     b3_n (float): Third virial coefficient for neutrons.
     db3_n (float): Derivative of the third virial coefficient for neutrons.
     Returns:
     dict: Dictionary containing the calculated thermodynamic properties for neutrons.
     """
     xmu_n2, out_b         = thermodynamics.thermo_neutron(T,rho,b2_n,db2_n,b3_n,db3_n,True,True)
     return out_b
