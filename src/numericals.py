import numpy as np
from scipy import integrate 
import matplotlib.pyplot as plt


"""
This module provides numerical integration functions and utilities for scientific computations.
Functions:
    integrand_ER(x, eta, idx):
        Computes the integrand for the ER function.
    integrand(x, eta, theta, idx):
        Computes the integrand for the given parameters.
    integrand_prime(x, eta, theta, idx):
        Computes the derivative of the integrand for the given parameters.
    GLQ(prime, idx, theta, eta):
        Performs Gaussian-Laguerre Quadrature (GLQ) integration.
    safe_exp(x):
        Computes the exponential function with safety checks to avoid overflow or underflow.
"""

def integrand_ER(x,eta,idx):
    return x**idx/(1+np.exp(x-eta))

def integrand(x,eta,theta,idx):
    return x**idx*(1+0.5*theta*x)**0.5/(1+np.exp(x-eta, dtype='longdouble'))

def integrand_prime(x,eta,theta,idx):
    return x**idx*(1+0.5*theta*x)**0.5/((2+safe_exp(x-eta)+safe_exp(eta-x)))

def GLQ(prime,idx,theta,eta):
    if prime == 0:
        result = integrate.quad(integrand,0, 5000, args=(eta, theta, idx))
        return result[0]
    if prime == 1:
        result = integrate.quad(integrand_prime,0, 5000,args=(eta, theta, idx))
        return result[0]



def safe_exp(x):
    # Set a threshold to avoid overflow or underflow
    threshold = 700  # You can adjust this based on your specific needs
    
    # Clip the input values to the threshold
    clipped_x = np.clip(x, -threshold, threshold)

    # Apply the exponential function to the clipped values
    result = np.exp(clipped_x)

    return result

'''
def GAUSS_LAGUERRE(n,alf):
    max_it = 100
    eps    = 1e-14
    x      = np.zeros(n, dtype=np.float64)
    w      = np.zeros(n, dtype=np.float64)
    for i in range(n):
        if i == 0:
            #GUESS FOR THE SMALLEST ROOT
            z = (1+alf)*(3+0.92*alf)/(1+2.4*float(n)+1.8*alf)
        elif i == 1:
            #GUESS FOR THE SECOND ROOT
            z = z + (15+6.25*alf)/(1+0.9*alf+2.5*float(n))
        else:
            #GUESS FOR ROOTS 3 TO N
            ai = float(i-1)
            z = z + ((1+2.55*ai)/(1.9*ai) + 1.26*ai*alf/(1+3.5*ai)) \
            *( z-x[i-2] )/(1+0.3*alf)

        #Newthon-Raphson for roots of Laguerre polynomial
        for its in range(max_it):
            p1 = 1
            p0 = 0
            ## RECURREENCE LOOP FOR LAGUERRE POLYNOMIAL AT Z
            for j in range(n):
                j  +=1
                p2 = p1
                p3 = p2
                p1 = ((2*float(j)-1+alf-z)*p2-(float(j-1)+alf)*p3)/float(j)
            #p1 is now the desired Laguerre polynomial. 
            #      We then compute pp, its derivative, from recurrence
            pp=(float(n)*p1-(float(n)+alf)*p2)/z
            z1=z
            z=z1-p1/pp # Newton formula
            if (abs(z-z1) <= eps): break
        # End N-R loop
        if(its >= max_it):
            print('too many iterations for Gauss-Laguerre')
            exit()
        # Store the root z and the weight w.
        x[i] = z 
        w[i] = -np.exp(gammln(alf+n)-gammln(float(n)))/(pp*float(n)*p2)


    return x, w
'''

### I trust for now this formulation although I did not go into details
def gammln(xx):
    cof = [76.18009172947146, -86.50532032941677, 24.01409824083091, -1.231739572450155, 0.01208650973866179, -.5395239384953E-5]
    stp =  2.5066282746310005

    x=xx
    y=x
    tmp=x+5.5
    tmp=(x+0.5)*np.log(tmp)-tmp
    ser=1.000000000190015
    for j in range(6):
        y=y+1
        ser=ser+cof[j]/y
 
    gamma=tmp+np.log(stp*ser/x)


    return gamma

