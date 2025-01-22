import numpy as np

hbc   =197.3269718   #h*c  [MeV*fm] where h is the slashed planck constant   
hbc2  =hbc*hbc
hbc3  =hbc2*hbc
pi    = np.pi    #3.14         #pi greco
pi2   = pi*pi   
c     = 1


xmassn=939.56536     #proton mass
xmassp=938.27203     #neutron mass
xmasse=0.5109989     #electron mass
xmass =xmassn        #should be (xmassn+xmassp)/2

htm   =[hbc2/xmassn, hbc2/xmassp]   #
htmm  =(htm[0]+htm[1])/2
zi    =[0,1]

ener_deutero=2.22   #deuteron binding energy
deg   =2            #isospin degeneracy

xn0   = np.sqrt(2) * (xmasse/hbc)**3/pi2

#Scattering lengths
a_np=-23.768  # fm
reff_np=2.68  #
a_nn=-18.5    # fm
reff_nn=2.86  #

########### INITIAL CONTIDIONS ############
###    take in input rho_i, drho, N samples of rho. 		###
rho_i = 1e-5
nrho  = 101 
rho_f= 0.1
#rhos = np.logspace(np.log10(rho_i), np.log10(rho_f), nrho, endpoint=True)
rhos = np.logspace(-7, -1, 101, endpoint=True)

###    same for proton fraction xxp				###
xxp_in = 0.001
xxp_f  = 0.999


###    same for temperature t					###
t_i = 5
t_f = 50
###    read phase shifts, if yes than read it from the file,	###
read_ps = True
###	  otherwise use only scattering length and reff		###

###    Read J_max
J_max = 8
Nch=4


itermax = 10000



#Global quantities




# OUTPUT FILENAMES
h5_name               = 'data//data.h5'
npy_name              = 'data//data'
npy_name_neutron      = 'data//data_neutron'
filename_barion_TD    = 'data//thermal_barion.dat' 
filename_electrons_TD = 'thermal_electron.dat'
filename_virial       = 'data//virial_coeff.dat'
filename_phase        = 'data//phase_shifts.dat'


