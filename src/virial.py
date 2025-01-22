import numpy as np
import src.numericals as numericals
from scipy.interpolate import CubicSpline
from scipy.integrate import quad
from src.constants import *

def func0(x,spl,T):
    if 2*T*x<350:
        return np.exp(-x)*spl(2*T*x)*pi/180
    else:
        return np.exp(-x)*spl(350)*pi/180
def func1(x,spl,T):
    if 2*T*x<350:
        return np.exp(-x)*spl(2*T*x)*x*pi/180
    else:
        return np.exp(-x)*spl(350)*(350/(2*T))*pi/180
    

def  virial(T, energy, phase_shift_total_nn, phase_shift_total_np, print_virial, interaction = True):
    # Calculation of the virial coefficients
    #filename_virial = './data/virial_coeff.dat'
  if interaction:
    #Interpolation of the phase shifts
    spl_np         = CubicSpline(energy,phase_shift_total_np)
    spl_nn         = CubicSpline(energy,phase_shift_total_nn)

    # Calculation of b2_n
    b2d_n=np.sqrt(2)/pi*quad(func0, 0, 1000, args=(spl_nn,T), limit=250)[0]           
    b2_n=b2d_n - 2**(-2.5)

    # Calculation of db2_n, derivative of b2_n
    db2d_n=np.sqrt(2)/pi*quad(func1, 0, 1000, args=(spl_nn,T), limit=250)[0]          
    db2_n=(db2d_n-b2d_n)/T
      
    b2d_nuc=1/(np.sqrt(8)*pi)*quad(func0, 0, 1000, args=(spl_np,T), limit=250)[0]     
    b2_nuc=b2d_nuc - 2**(-2.5) + 3/np.sqrt(2)*(np.exp(ener_deutero/T)-1) 
      
    db2d_nuc=1/(np.sqrt(8)*pi)*quad(func1, 0, 1000, args=(spl_np,T), limit=250)[0]    
    db2_nuc=(db2d_nuc-b2d_nuc)/T - 3*ener_deutero/(np.sqrt(2)*T**2)*np.exp(ener_deutero/T)
    
    # Calculation of b2_np and its derivative
    b2_np= b2_nuc - b2_n
    db2_np= db2_nuc - db2_n

    # 
    Gamma2 = T*db2_n*(T*db2_n*2/3 + b2_n)
    if print_virial:
        m = open(filename_virial, 'a')
        m.write('\t{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}\n'.format(round(T,4),round(b2_n,4),round(T*db2_n,4),round(b2_np,4),round(T*db2_np,4), round(b2_nuc,4), round(T*db2_nuc,4), round(Gamma2, 4)))
        #m2.write('{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\n'.format(t, b2_n[counter], t*db2_n[counter], b2_np[counter], t*db2_np[counter], Gamma2[counter], 0))
    return b2_n, db2_n, b2_np, db2_np
  else:
    ### Free fermi gas
    #b2_n= -2**(-2.5)
    #db2_n = 0
    #b2_np = 0
    #db2_np = 0
    #b2_nuc = 0
    #db2_nuc = 0
    #Gamma2 = 0

    ## Only deuteron contribution
    b2_n= -2**(-2.5)
    db2_n = 0
    b2_nuc = 3/np.sqrt(2)*(np.exp(ener_deutero/T)-1) - 2**(-2.5)
    db2_nuc = -3*ener_deutero/(np.sqrt(2)*T**2)*np.exp(ener_deutero/T)
    b2_np = b2_nuc - b2_n
    db2_np = db2_nuc
    Gamma2 = 0
    return b2_n, db2_n, b2_np, db2_np