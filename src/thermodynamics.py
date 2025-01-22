import numpy as np
from src.constants import *    
import src.numericals as numericals
import scipy
from scipy import integrate 
import matplotlib.pyplot as plt
import mpmath

from time import sleep

def safe_exp(x):
    # Set a threshold to avoid overflow or underflow
    threshold = 700  # You can adjust this based on your specific needs
    
    # Clip the input values to the threshold
    clipped_x = np.clip(x, -threshold, threshold)

    # Apply the exponential function to the clipped values
    result = np.exp(clipped_x)

    return result

def F(xz, xz0,b2_n,db2_n,b2_np,db2_np):
    """
    Computes equation B4 from Rivieccio et al 2025.
    Parameters:
        xz : The current values for the roots and weights.
        xz0 : The initial values for the roots and weights.
        b2_n : The second virial coefficient for neutrons.
        db2_n : The derivative of the second virial coefficient for neutrons.
        b2_np : The second virial coefficient for neutron-proton interactions.
        db2_np : The derivative of the second virial coefficient for neutron-proton interactions.
    """
    return np.array([
        xz[0] + 2.0 * b2_n * xz[0]**2 + 2.0 * xz[0] * xz[1] * b2_np - xz0[0],
        xz[1] + 2.0 * b2_n * xz[1]**2 + 2.0 * xz[0] * xz[1] * b2_np - xz0[1]
        ])

def zero_temperature_leptons(rho_lep, rho_barions):
        '''
        Compute the zero temperature pressure and energy density for a free fermi gas.
        Parameters:
            rho_lep: The lepton density.
            rho_barions: The barion density.
        '''
        #'''
        # Direct calculation of the zero temperature pressure and energy density
        #  for a free fermi gas where xmasse is the mass of the electron
        #  and rho_lep is the lepton density
        xkf   = (3*pi2*rho_lep)**(1/3)
        xx    = hbc*xkf/xmasse   
        x2    = xx**2

        ebar  =xn0*xmasse/8/np.sqrt(2)
        p0    = ebar/3*( xx*np.sqrt(1+x2)*(-3+2*x2) + 3*np.arcsinh(xx) )
        e0    = ebar*( xx*np.sqrt(1+x2)*(1+2*x2) - np.arcsinh(xx) ) / rho_barions
        #'''
        '''  FAUSSUIER 2016, Eq 129-130
        ebar   = xmasse**4/(8*pi2*hbc3)
        pf     = hbc*(3*pi**2*rho_lep)**(1/3)
        thetaf = np.sqrt((pf/m)**2 + 1) - 1
        
        e0 = ebar*(np.sqrt(thetaf*(thetaf + 2)) * (1-thetaf/3 + 10*thetaf**2/3 + 2*thetaf**3) - 2*np.arcsinh(np.sqrt(thetaf/2)) )
        p0 = ebar*(np.sqrt(thetaf*(thetaf + 2)) * (1+thetaf) *(-1+4*thetaf/3 + 2*thetaf**2/3) + 2*np.arcsinh(np.sqrt(thetaf/2)) )
        
        #'''
        return e0, p0

def Lattimer_eta(T, rho , xxp):
        '''
        Compute the chemical potential for electrons and positrons using the Lattimer formula.
        Parameters:
            T: The temperature.
            rho: The density.
            xxp: The proton fraction.
        Returns:
            mu/T :degeneracy parameter
        '''
        tt = 3*pi2*hbc3*rho*xxp/2                           #xxp is the proton fraction
        qq=(pi*T)**2/3 - xmasse**2*0.5		                #
        rr=(np.sqrt( qq**3 + tt**2 ) + tt)**(1/3)	        #
        xmu0=rr-qq/rr						                #
        return xmu0/T	
    
def elec_n_dens_eq(eta,theta, rho, xxp, ER):
    '''
    Compute the electron density for a given chemical potential.
    Parameters:
        eta: The degeneracy parameter.
        theta: The adimensional temperature.
        rho: The density (number density).
        xxp: The proton fraction.
        ER: The extreme relativistic flag.
    Returns:
        The value of the charge neutrality equation.
    '''
    if theta >= 1:# Theshhold Positron Calculation
        if ER==False: #if not extreme relativistic
            tmp12 = numericals.GLQ(prime = 0, idx = 0.5, theta = theta, eta = eta)
            tmp32 = numericals.GLQ(prime = 0, idx = 1.5, theta = theta, eta = eta)
            etap  = -eta-2/theta
            tmp12p = numericals.GLQ(prime = 0, idx = 0.5, theta = theta, eta = (etap))            #Positrons
            tmp32p = numericals.GLQ(prime = 0, idx = 1.5, theta = theta, eta = (etap))            #Positrons
        else: #if extreme relativistic
            Li2   =float(mpmath.polylog(2, -np.exp(eta)))                                             # Polylogarithm
            Li3   =float(mpmath.polylog(3, -np.exp(eta)))                                             # Polylogarithm
            tmp12 = np.sqrt(theta/2)*(-Li2)                                                            #E.R. limit
            tmp32 = np.sqrt(theta/2)*(-2 * Li3)                                                        #E.R. limit
            etap  = -eta-2/theta
            Li2   =float(mpmath.polylog(2, -np.exp(etap)))                                             # Polylogarithm
            Li3   =float(mpmath.polylog(3, -np.exp(etap)))                                             # Polylogarithm
            tmp12p = np.sqrt(theta/2)*(-Li2)                                                            #E.R. limit
            tmp32p = np.sqrt(theta/2)*(-2 * Li3)                                                        #E.R. limit
        n_e   = xn0*theta**1.5*(tmp12 + theta*tmp32)
        n_p   = xn0*theta**1.5*(tmp12p + theta*tmp32p)
        return  n_e-rho*xxp-n_p
    else: # NO positron calculation
        if ER==False:
            tmp12 = numericals.GLQ(prime = 0, idx = 0.5, theta = theta, eta = eta)
            tmp32 = numericals.GLQ(prime = 0, idx = 1.5, theta = theta, eta = eta)
        else:
            Li2   =float(mpmath.polylog(2, -safe_exp(eta)))                                             # Polylogarithm
            Li3   =float(mpmath.polylog(3, -safe_exp(eta)))                                             # Polylogarithm
            tmp12 = np.sqrt(theta/2)*(-Li2)                                                            #E.R. limit
            tmp32 = np.sqrt(theta/2)*(-2 * Li3)                                                        #E.R. limit
        n_e   = xn0*theta**1.5*(tmp12 + theta*tmp32)
        return  n_e-rho*xxp

def elec_n_dens_prime(eta,theta, rho, xxp):   
        #Derivative of the charge neutrality equation with positron for Newton-Raphson
        # This is not used in the final calculation due to the numerical instability
        tmp12 = numericals.GLQ(prime = 1, idx = 0.5, theta = theta, eta = eta)
        tmp32 = numericals.GLQ(prime = 1, idx = 1.5, theta = theta, eta = eta)
        tmp12p = numericals.GLQ(prime = 1, idx = 0.5, theta = theta, eta = (-eta-2/theta))            #Positrons
        tmp32p = numericals.GLQ(prime = 1, idx = 1.5, theta = theta, eta = (-eta-2/theta))            #Positrons
        ele = (tmp12 + theta*tmp32)
        pos = (tmp12p + theta*tmp32p)                                                       #Positrons
        return xn0*theta**1.5*(ele - pos)                                                   #Positrons
    
def  thermo_e(T, rho, xxp, final_calcul):

    theta=T/xmasse

    
    #RELATIVITY PAR. FROM LATTIMER
    eta     = [Lattimer_eta(T,rho,xxp)]
    eta_lat = eta[0]

    a = 500*rho**(0.01)
    b = 1100*rho**(0.01)
    #print(a, b, theta)

    #RELATIVITY PARAMETER FROM FERMI DIRAC INTEGRAL EQUATION
    #'''
    if theta < 1:
            eta     = scipy.optimize.newton(func = elec_n_dens_eq, x0 = a,x1 = b, args= (theta, rho, xxp, False), full_output = True, maxiter = 500, tol = 1e-10)
            eta_full  = eta[0]
            #eta_er = eta_full
    else:
            eta     = scipy.optimize.newton(func = elec_n_dens_eq, x0 = -70,x1 = 500, args= (theta, rho, xxp, False), full_output = True, maxiter = 500, tol = 1e-10)
            eta_full  = eta[0]
            #eta     = scipy.optimize.newton(func = elec_n_dens_eq, x0 = -10,x1 = 500, args= (theta, rho, xxp, True), full_output = True, maxiter = 500, tol = 1e-10)
            #eta_er  = eta[0]


    #xmu_e_er  = eta_er*T#-xmasse
    #xmu_e_lat = eta_lat*T#-xmasse
    xmu_e = eta_full*T

    

    if (final_calcul):
        '''
        ### MY FORMULAS ###
        ebar  = xn0*xmasse
        #press = (2/3)*ebar*theta**3/(np.sqrt(2))*(theta*eta**4/4 + eta**3 + eta**2*(theta*pi*0.5+1/theta) + eta*(pi2+  1/theta**2) + pi2/(3*theta) - 2/(3*theta**3) + 7*theta*pi2**2/60) ##MY FORMULA
        press = (2/3)*ebar*theta**3/(np.sqrt(2))*(theta*eta**4/4 + eta**3 + eta**2*(theta*pi*0.5+1/theta) + eta*(pi2) + pi2/(3*theta) + 7*theta*pi2**2/60) ##MY FORMULA approx 1/theta**2
        ener = ebar*theta**3/(np.sqrt(2))*(theta*eta**4*0.5 + 2*eta**3 + eta**2*(4/theta+pi2*theta) + eta*(2*pi2) + 4*(pi2)/(3*theta) + 7*pi2**2*theta/30)/rho
        
        s_e   = T*xmu_e**2 * (1 + xmu_e**(-2)*(7*pi2 * T**2/15 - 0.5*xmasse**2))/(3*rho*hbc3)
        n_e   = xn0*theta**1.5*np.sqrt(theta/2)* (eta**2/2 + pi2/6 + 2*theta*(eta**3/6 + pi2*eta/6))

        e0, p0 = zero_temperature_leptons(n_e, rho)
        press_th = press - p0
        ener_th  = ener - e0

        press_th_p = 0#press_th
        ener_th_p  = 0#ener_th
        s_p     = 0       #the electron formula should already take into account the positrons 
        n_p     = n_e - xxp*rho
        #print(press, ener, n_p)
        e0, p0 = zero_temperature_leptons(n_p, rho)
        press_th = press - p0
        ener_th  = ener - e0
        dir = 'my_form/'
        #'''
       
        '''
        #### LATTIMER ELECTRONS AND POSITRONS ####
        
        coff  = 2 * xmu_e * (xmu_e/hbc)**3 / (24 * pi2)
        
        press = coff * (1 + xmu_e**(-2) * (2 * pi2 * T**2 - 3 * xmasse**2) + pi2 * T**2 * (7*pi2 * T**2 / 15 - 0.5 * xmasse**2)/(xmu_e**4)) 
        ener  = coff * 3 * (1 + xmu_e**(-2) * (2 * pi2 * T**2 - xmasse**2) + pi2 * T**2 * (7*pi2 * T**2 / 15 - 0.5 * xmasse**2)/(xmu_e**4)) / rho  
        s_e   = T*xmu_e**2 * (1 + xmu_e**(-2)*(7*pi2 * T**2/15 - 0.5*xmasse**2))/(3*rho*hbc3)
        n_e   = xn0*theta**1.5*np.sqrt(theta/2)* (eta**2/2 + pi2/6 + theta*(eta**3 + pi2*eta)/3)
        
        e0, p0 = zero_temperature_leptons(n_e, rho)
        press_th = press - p0
        ener_th  = ener - e0

        press_th_p = 0#press_th
        ener_th_p  = 0#ener_th
        s_p     = 0       #the electron formula should already take into account the positrons 
        n_p     = n_e - xxp*rho
        #print(press, ener, n_p)
        e0, p0 = zero_temperature_leptons(n_p, rho)
        press_th = press - p0
        ener_th  = ener - e0
        press_p = 0
        ener_p = 0
        dir = 'Lattimer/'
         
 
        #'''

        ###              ELECTRONS           ###
                 #RELATIVISTIC ELECTRONS 
        '''
        #Li2   =-eta**2/2 - pi2/6 - float(mpmath.polylog(2, -np.exp(-eta)))                                 # Polylogarithm   APPROX
        #Li3   =-eta**3/6 - pi2*eta/6 + float(mpmath.polylog(3, -np.exp(-eta)))                             # Polylogarithm   APPROX
        #Li4   =-eta**4/24 -pi2*eta**2/12 - 7*pi2**2/(360) - float(mpmath.polylog(4, -np.exp(-eta)))        # Polylogarithm   APPROX
        Li2   =float(mpmath.polylog(2, -np.exp(eta)))                                             # Polylogarithm
        Li3   =float(mpmath.polylog(3, -np.exp(eta)))                                             # Polylogarithm
        Li4   =float(mpmath.polylog(4, -np.exp(eta)))                                             # Polylogarithm
        tmp12 = np.sqrt(theta/2)*(-Li2)                                                            #E.R. limit
        tmp32 = np.sqrt(theta/2)*(-2 * Li3)                                                        #E.R. limit
        tmp52 = np.sqrt(theta/2)*(-6 * Li4)                                                        #E.R. limit

        ebar  = xn0*xmasse
        n_e   = xn0*theta**1.5*(theta*tmp32)
        press = ebar*theta**2.5*(theta*tmp52)*(1/3)
        ener  = press*3/rho
        s_e   = xn0*theta**1.5*(4*theta*tmp52/3 - eta * (theta*tmp32))/rho 

        #! ZERO TEMPERATURE ELECTRONS
        e0, p0 = zero_temperature_leptons(rho*xxp, rho)
        press_th = press - p0
        ener_th  = ener - e0
        #'''

     
        ###          POSITRONS              ###
                #RELATIVISTIC POSITRONS
        ''' 
        etap  = -eta-2/theta                                                               #effective eta for positrons
        #Li2   =-etap**2/2 - pi2/6 - float(mpmath.polylog(2, -np.exp(-etap)))                                 # Polylogarithm   APPROX
        #Li3   =-etap**3/6 - pi2*etap/6 + float(mpmath.polylog(3, -np.exp(-etap)))                             # Polylogarithm   APPROX
        #Li4   =-etap**4/24 -pi2*etap**2/12 - 7*pi2**2/(360) - float(mpmath.polylog(4, -np.exp(-etap)))        # Polylogarithm   APPROX
        Li2p   =float(mpmath.polylog(2, -np.exp(etap)))#integrate.quad(numericals.integrand_ER,0, 700, args=(etap, 1)) #                                    # Polylogarithm
        Li3p   =float(mpmath.polylog(3, -np.exp(etap)))#integrate.quad(numericals.integrand_ER,0, 700, args=(etap, 2)) #                                    # Polylogarithm
        Li4p   =float(mpmath.polylog(4, -np.exp(etap)))#integrate.quad(numericals.integrand_ER,0, 700, args=(etap, 3)) #                                    # Polylogarithm
        #print(Li2, Li2p, Li3, Li3p, Li4,Li4p)
        tmp12p = np.sqrt(theta/2)*(-Li2p)#Li2p[0]#                                                            #E.R. limit
        tmp32p = np.sqrt(theta/2)*(-2 * Li3p)#Li3p[0]#                                                        #E.R. limit
        tmp52p = np.sqrt(theta/2)*(-6 * Li4p)#Li4p[0]#                                                        #E.R. limit
        dir = 'ER/'

        press_p = ebar*theta**2.5*(theta*tmp52p)*(1/3)
        print(press_p)
        ener_p  = press*3/rho    
        s_p     = xn0*theta**1.5*(4*theta*tmp52p/3 + eta * (theta*tmp32p))/rho     # similar to the previous but just s
        n_p     = xn0*theta**1.5*(theta*tmp32p)
        
        #! ZERO TEMPERATURE ELECTRONS
        #e0, p0 = zero_temperature_leptons(rho*xxp, rho)
        press_th_p = press_p 
        ener_th_p  = ener_p
        #'''




        #'''     #GENERAL ELECTRONS
        tmp12 = numericals.GLQ(prime = 0, idx = 0.5, theta = theta, eta = eta_full)     			   #Full fomulation
        tmp32 = numericals.GLQ(prime = 0, idx = 1.5, theta = theta, eta = eta_full)    			       #Full fomulation
        tmp52 = numericals.GLQ(prime = 0, idx = 2.5, theta = theta, eta = eta_full)				       #Full fomulation
        

        ebar  = xn0*xmasse
        n_e   = xn0*theta**1.5*(tmp12 + theta*tmp32)
        press = (2/3)*ebar*theta**2.5*(tmp32 + theta*0.5*tmp52)
        ener  = (ebar*theta**2.5*(tmp32 + theta*tmp52) + xmasse*n_e)/rho
        s_e   = xn0*theta**1.5*(5*tmp32/3 + 4*theta*tmp52/3 - eta_full * (tmp12 + theta*tmp32))/rho    

        
        # ZERO TEMPERATURE ELECTRONS
        e0, p0 = zero_temperature_leptons(rho*xxp, rho)
        press_th = press - p0
        ener_th  = ener - e0
        #'''
        
   

               #GENERAL POSITRONS
        #'''
        etap = -eta_full - 2/theta
        tmp12p =numericals.GLQ(prime = 0, idx = 0.5, theta = theta, eta = (etap))         #Full formulation
        tmp32p =numericals.GLQ(prime = 0, idx = 1.5, theta = theta, eta = (etap))         #Full formulation
        tmp52p =numericals.GLQ(prime = 0, idx = 2.5, theta = theta, eta = (etap))         #Full formulation
        
        n_p     = xn0*theta**1.5*(tmp12p + theta*tmp32p)
        press_p = (2/3)*ebar*theta**2.5*(tmp32p + theta*0.5*tmp52p)
        ener_p  = (ebar*theta**2.5*(tmp32p + theta*tmp52p) + xmasse*n_p)/rho 
        s_p     = xn0*theta**1.5*(5*tmp32p/3 + 4*theta*tmp52p/3 + etap * (tmp12p + theta*tmp32p))/rho     # similar to the previous but just s
        
        #positrons don't have a cold contribution!
        
        P_ph=pi2*T**4/(hbc3*45)   ##### taking k_b = c = 1
        E_ph=3*P_ph/rho

        '''
        m1 = open( filename_electrons_TD, 'a')
        m1.write('{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}\t{8}\t{9}\t{10}\t{11}\t{12}\t{13}\t{14}\t{15}\t{16}\t{17}\n'.format(rho,xxp,int(T),xmu_e,ener,ener_th,press, press_th, E_ph, P_ph, s_e, n_e, ener_th_p, press_th_p, s_p, n_p, press_p, ener_p))
        m1.close()
        '''
        out = {}
        out['rho']       =rho
        out['mu_e']      =xmu_e
        out['ener']      =ener
        out['press']     =press
        out['ener_th']   =ener_th
        out['press_th']  =press_th
        out['s_e']       =s_e
        out['n_e']       =n_e
        out['ener_p']    =ener_p
        out['press_p']   =press_p
        out['s_p']       =s_p
        out['n_p']       =n_p
        out['P_ph']      =P_ph
        out['E_ph']      =E_ph

        return xmu_e, out
    return xmu_e



### BARION THERMODYNAMICS
def  thermo(T,rho,xxp,b2_n,db2_n,b2_np,db2_np,calcul,escriure):
    #output is xmu_n, xmu_p
    Ratio= -0.5022
    rho_n=rho*(1-xxp)       #neutron density
    rho_p=rho*xxp           #proton density

    #define the DeBroglie wavelenght in fm
    xL=np.sqrt(2*pi/(xmass*T) )*hbc
    xL3=xL**3

    ###In order to find the fugacity z=e^{\mu/T}
    ###       we use the Newton-Raphson method (tangent root finding)
    ###       to solve equation 21 of the paper
    ### x_new = x_old + dx = x_old - J(x_old)^-1 * f(x_old)
    ###       where J^-1 is the inverse of the jacobian that is analog to f'
    ###In this case the function we use is the number density of proton and neutron
    
    #eq 21 of the paper, xz0 is the left hand
    xz0 = np.array([rho_n*xL3/2,rho_p*xL3/2])    #initialize lowest order for fugacity
    xz = xz0
    #print('Barions Fugacity =', xz)
    #print(xz[1])
    dz = np.array([10,10])                       #initialize the spacing 
    iter = 0  
    #exit()

    #eq 21 of the paper, F is the right hand
    #initialize
    J_adj = np.zeros((2,2))
    #'''
    xz = scipy.optimize.root(F, x0 = xz0, args=(xz0,b2_n,db2_n,b2_np,db2_np), tol = 1e-15)#, maxiter = 10000)
    xz = xz.x
    #'''
    '''
    while np.max([np.abs(dz[0]),np.abs(dz[1])] )> 1e-12:# and iter<itermax:# and iter <itermax:#14:# and iter < 10000:#14:# 
        F = np.array([
        xz[0] + 2.0 * b2_n * xz[0]**2 + 2.0 * xz[0] * xz[1] * b2_np - xz0[0],
        xz[1] + 2.0 * b2_n * xz[1]**2 + 2.0 * xz[0] * xz[1] * b2_np - xz0[1]
        ])

        J_adj = np.array([
        [1.0 + 4.0 * xz[1] * b2_n + 2.0 * xz[0] * b2_np, -2.0 * xz[1] * b2_np],
        [-2.0 * xz[0] * b2_np, 1.0 + 4.0 * xz[0] * b2_n + 2.0 * xz[1] * b2_np]
        ])

        J_det = 1.0 + (4.0 * b2_n + 2.0 * b2_np) * (xz[0] + xz[1]) + 8.0 * b2_n * b2_np * (xz[0]**2 + xz[1]**2) + 16.0 * xz[0] * xz[1] * (b2_n**2)

        dz = np.dot(J_adj,F)/J_det    #could use also np.matmul
        #dz = np.matmul(J_adj,F)/J_det    #could use also np.matmul
        
        xzp=xz - dz
        xz=xzp
        
        iter=iter+1
        if xz[1]<1e-150: xz[-1] = 0
    #'''    
    ###################
    
    if calcul == False:# and iter<itermax:
       escriure = False
       rhoz_p   = 2/xL3*( xz[1] + 2*xz[1]**2*b2_n + 2*xz[1]*xz[0]*b2_np )

       xmu_p    = T*np.log(xz[1])
       rhoz_n   = 2/xL3*( xz[0] + 2*xz[0]**2*b2_n + 2*xz[1]*xz[0]*b2_np )
       xmu_n    = T*np.log(xz[0])

    elif (calcul):# and iter < itermax:   # ! .and. np.max(xz(1),xz(2)) < 0.5 ####we are excluding the iteration that don't converge, mh interesting
       escriure=True
    
       rhoz_p   = 2/xL3*( xz[1] + 2*xz[1]**2*b2_n + 2*xz[1]*xz[0]*b2_np )
       xmu_p    = T*np.log(xz[1])
       rhoz_n   = 2/xL3*( xz[0] + 2*xz[0]**2*b2_n + 2*xz[1]*xz[0]*b2_np )
       xmu_n    = T*np.log(xz[0])

       #! ... PRESSURE
       press    = (2*T/xL3)*( xz[0]+xz[1] + (xz[0]**2+xz[1]**2)*b2_n + 2*xz[0]*xz[1]*b2_np )

       #! ... ENTROPY
       entro    = (5/2)*(press/T) - rhoz_n*np.log(xz[0]) - rhoz_p*np.log(xz[1]) \
         + (2*T/xL3)*(  (xz[0]**2+xz[1]**2)*db2_n + 2*xz[0]*xz[1]*db2_np )


       #! ... ENERGY
       ener     = (3/2)*press+(2*T**2/xL3)*((xz[0]**2+xz[1]**2)*db2_n + 2*xz[0]*xz[1]*db2_np)    

       #! ... FREE
       free     = ener-T*entro

       '''
    #! ... THERMAL INDEX
       b2p_n    = T*db2_n
       b2p_np   = T*db2_np
       th_c0    = 4/3
       th_c1    = T/xL3
       th_c21   = -3*T**2/xL3**2
       th_c22   = T/xL3 *(b2_n+(2/3)*db2_n*T)
       th_c31   = -6*T**2/xL3**2
       th_c32   = T/xL3 *(2*b2_np+(4/3)*db2_np*T)

       S_E      =((3/8)*T*xL3*(b2_n-b2_np) + (xL3*T**2)*(db2_n-db2_np)/4)*rho
       S_F      = T*np.log(2)+ T*xL3*(b2_n-b2_np)*rho
       '''

       '''
       m_bar = open(filename_barion_TD, 'a')
       m_bar.write('{:.2E}\t{:.7E}\t{:.7E}\t{:.7E}\t{:.7E}\t{:.7E}\t{:.7E}\t{:.7E}\t{:.7E}\t{:.7E}\t{:.7E}\n'.format(rho,xxp,int(T),xmu_n,xmu_p,ener/rho,entro/rho,free/rho,press,S_E/rho,S_F/rho))
       m_bar.close()
       '''

       out = {}
       out['mu_p']      =xmu_p
       out['rho']       =rho
       out['mu_n']      =xmu_n
       out['ener']      =ener/rho
       out['press']     =press
       out['entro']       =entro/rho
       out['free']       =free/rho
       #if iter<itermax:
       return xmu_n, xmu_p , out
    #print('finished barions')
    if iter<itermax and calcul == False:
        return xmu_n, xmu_p
    else:
        return 0.1,0.1







def F_neutron(xz, xz0,b2,b3):
    return xz + 2.0 * b2 * xz**2 + 3*xz**3*b3 - xz0


### BARION THERMODYNAMICS
def  thermo_neutron(T,rho,b2_n,db2_n,b3_n,db3_n,calcul,escriure):
    #output is xmu_n, xmu_p
    Ratio= -0.5022

    #define the DeBroglie wavelenght in fm
    xL=np.sqrt(2*pi/(xmassn*T) )*hbc
    xL3=xL**3

    ###In order to find the fugacity z=e^{\mu/T}
    ###       we use the Newton-Raphson method (tangent root finding)
    ###       to solve equation 21 of the paper
    ### x_new = x_old + dx = x_old - J(x_old)^-1 * f(x_old)
    ###       where J^-1 is the inverse of the jacobian that is analog to f'
    ###In this case the function we use is the number density of proton and neutron
    
    #eq 21 of the paper, xz0 is the left hand
    xz0 = rho*xL3/2    #initialize lowest order for fugacity
    xz = xz0
    
    xz = scipy.optimize.root(F_neutron, x0 = xz0, args=(xz0,b2_n,b3_n), tol = 1e-15)#, maxiter = 10000)
    xz = xz.x[0]
      
    ###################
    
    if calcul == False:# and iter<itermax:
       escriure = False
       rhoz_n   = 2/xL3*( xz + 2*xz**2*b2_n + 2*xz**3*b3_n )
       xmu_n    = T*np.log(xz)

    elif (calcul):# and iter < itermax:   # ! .and. np.max(xz(1),xz(2)) < 0.5 ####we are excluding the iteration that don't converge, mh interesting
       escriure=True
       rhoz_n   = 2/xL3*( xz + 2*xz**2*b2_n + 2*xz**3*b3_n )
       xmu_n    = T*np.log(xz)
       
       #! ... PRESSURE
       press    = (2*T/xL3)*( xz + xz**2*b2_n + xz**3*b3_n )

       #! ... ENTROPY
       entro    = (5/2)*(press/T) - rhoz_n*np.log(xz) + (2*T/xL3)*(  (xz**2)*db2_n + 2*xz**3*db3_n )


       #! ... ENERGY
       ener     = (3/2)*press+(2*T**2/xL3)*((xz**2)*db2_n + xz**3*db3_n)    

       #! ... FREE
       free     = ener-T*entro

       '''
    #! ... THERMAL INDEX
       b2p_n    = T*db2_n
       b2p_np   = T*db2_np
       th_c0    = 4/3
       th_c1    = T/xL3
       th_c21   = -3*T**2/xL3**2
       th_c22   = T/xL3 *(b2_n+(2/3)*db2_n*T)
       th_c31   = -6*T**2/xL3**2
       th_c32   = T/xL3 *(2*b2_np+(4/3)*db2_np*T)

       S_E      =((3/8)*T*xL3*(b2_n-b2_np) + (xL3*T**2)*(db2_n-db2_np)/4)*rho
       S_F      = T*np.log(2)+ T*xL3*(b2_n-b2_np)*rho
       '''


       out = {}
       out['mu_n']      =xmu_n
       out['ener']      =ener/rho
       out['press']     =press
       out['entro']     =entro/rho
       out['free']      =free/rho
       out['b2n']       = b2_n
       out['b2n_p']     = db2_n
       #if iter<itermax:
       return xmu_n, out
    #print('finished barions')
    if iter<itermax and calcul == False:
        return xmu_n
    else:
        return 0.1,0.1








