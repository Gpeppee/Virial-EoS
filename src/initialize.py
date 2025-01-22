from src.constants import *
import os
import h5py
import pickle

def initialize_files():
    m_b = open(filename_barion_TD,'a')
    m_e = open(filename_electrons_TD, 'a')
    m_e.write('#rho,xxp,T,xmu_e,ener,ener_th,press, press_th, E_ph, P_ph, s_e, n_e, ener_th_p, press_th_p, s_p, n_p\n')
    m_b.write('#rho\t xxp\t T\t xmu_n\t xmu_p\t ener/rho\t entro/rho\t free/rho\t press\t S_E/rho\t S_F/rho\n')
    m_b.close()
    m_e.close()
    #m = {}
    #np.save(npy_name, m)
    

def interrupt_dataline(filename):
    m = open(filename,'a')
    m.write('\n\n\n')
    m.close()

def rename(xxp):
    filename_barion_TD0 = filename_barion_TD.replace('.dat', str(xxp) + '.dat')
    filename_electrons_TD0 = filename_electrons_TD.replace('.dat', str(xxp) + '.dat')
    
    os.rename(filename_barion_TD, filename_barion_TD0)
    os.rename(filename_electrons_TD, filename_electrons_TD0)

def WriteDictionaryVar(h5f,data) :
    for varname in sorted(list(data)) :
        x_data = data[varname]
        x_type = type(x_data)
        if x_type == dict :
            grp = h5f.create_group(varname)
            WriteDictionaryVar(grp,x_data) 
        else:
            if x_type != np.ndarray :
                values = np.array(x_data)
            else :
                values = x_data
            try :
                h5f.create_dataset(varname, data=values)
            except:
                print ("Error creating dataset ",varname)

def WriteScalarHDF5(file,data,group='',mode='a') :
    f = h5py.File(file,mode)
    if len (group) > 0 :
       try:
           grp = f.create_group(group)
       except : 
           grp = f[group]
       grp = f[group]
    else :
       grp = f
    WriteDictionaryVar(grp,data)
    f.close()
    return 1

def store_vectorized(filename, h5f=False):
    if type(filename).__name__=='str':
        with open(filename, 'rb') as handle:
            file = pickle.load(handle)
        TT   = filename.replace('.npy','').replace('data\\y_p\\data','')
        T    = int(TT)
        
        XXPP = list(file.keys())
        RHO  = list(file[XXPP[0]].keys())
        rho  = np.array(RHO).astype(float)
        xp   = np.array(XXPP).astype(float)/1000
        dim = np.size(rho)

        list_var_b = list(file[XXPP[0]][list(file[XXPP[0]].keys())[0]]['b'].keys())
        list_var_e = list(file[XXPP[0]][list(file[XXPP[0]].keys())[0]]['e'].keys())
        list_var_b = list_var_b
        list_var_e = list_var_e    

        
        varr={}
        for xxp in XXPP:
            varr[xxp]={}
            varr[xxp]['e']={}
            varr[xxp]['b']={}
            for var in list_var_e:
                varr[xxp]['e'][var] = np.zeros(dim)
            for var in list_var_b:
                varr[xxp]['b'][var] = np.zeros(dim)

        for xxp in XXPP:
            for count,r in enumerate(RHO):
                for var in list_var_e:
                    varr[xxp]['e'][var][count] = file[xxp][r]['e'][var]
                for var in list_var_b:
                    varr[xxp]['b'][var][count] = file[xxp][r]['b'][var]
        if h5f==False:
            new_name = 'data/vec_data' + TT +'.npy'
            with open(new_name, 'wb') as handle:
                pickle.dump(varr, handle)
        else:
            WriteScalarHDF5(h5f,varr,TT)


    elif type(filename).__name__=='list':
      for fname in filename:
        with open(fname, 'rb') as handle:
            file = pickle.load(handle)
        TT   = fname.replace('.npy','').replace('data','')
        T    = int(TT)
        XXPP = list(file.keys())
        RHO  = list(file[XXPP[0]].keys())
        rho  = np.array(RHO).astype(float)
        xp   = np.array(XXPP).astype(float)/1000
        dim = np.size(rho)

        list_var_b = list(file[XXPP[0]][list(file[XXPP[0]].keys())[0]]['b'].keys())
        list_var_e = list(file[XXPP[0]][list(file[XXPP[0]].keys())[0]]['e'].keys())
        list_var_b = list_var_b
        list_var_e = list_var_e    

        
        varr={}
        varr[TT]={}
        for xxp in XXPP:
            varr[TT][xxp]={}
            varr[TT][xxp]['e']={}
            varr[TT][xxp]['b']={}
            for var in list_var_e:
                varr[TT][xxp]['e'][var] = np.zeros(dim)
            for var in list_var_b:
                varr[TT][xxp]['b'][var] = np.zeros(dim)

        for xxp in XXPP:
            for count,r in enumerate(RHO):
                for var in list_var_e:
                    varr[TT][xxp]['e'][var][count] = file[xxp][r]['e'][var]
                for var in list_var_b:
                    varr[TT][xxp]['b'][var][count] = file[xxp][r]['b'][var]
        new_name = 'data/vec_data.npy'
        with open(new_name, 'wb') as handle:
            pickle.dump(varr, handle)
    #AA=varr[XXPP[0]]['b']['mu_n']
    #try:
    #    print(list(AA.keys()))
    #except:
    #    print(AA)
    
            
            
def save(file,data, T, xxp, rho, type):
    if type == 'y_p':
        TT = str(int(T))
        XXPP = str(int(xxp*1000))
        RHO = str(rho)
        file = 'data\\y_p\\data'+TT+'.npy'
        if rho==rhos[0] and xxp == 0.01:
            old = {}
            old[XXPP] = {}
            old[XXPP][RHO]=data
            with open(file, 'wb') as handle:
                pickle.dump(old, handle)
        elif rho==rhos[0] and xxp!=0.01:
            
            with open(file, 'rb') as handle:
                #print('AAAAAAAAAAAAAA',TT, XXPP, RHO)
                old = pickle.load(handle)
            old[XXPP] = {}
            old[XXPP][RHO]=data
        else:
            with open(file, 'rb') as handle:
                #print('AAAAAAAAAAAAAA',TT, XXPP, RHO)
                old = pickle.load(handle)
            old[XXPP][RHO]=data#{}
        with open(file, 'wb') as handle:
            pickle.dump(old, handle)
    elif type == 'beta':
        TT = str(int(T))
        XXPP = str(int(xxp*1000))
        RHO = str(rho)
        file = 'data\\beta\\data'+TT+'.npy'
    
        if rho==rhos[0]:
            old = {}
            old[RHO] = {}
            old[RHO][XXPP]=data
            with open(file, 'wb') as handle:
                pickle.dump(old, handle)
        elif rho!=rhos[0]:
            with open(file, 'rb') as handle:
                old = pickle.load(handle)
            old[RHO] = {}
            old[RHO][XXPP]=data
        #else:
        #    with open(file, 'rb') as handle:
        #        old = pickle.load(handle)
        #    old[RHO][XXPP]=data
        with open(file, 'wb') as handle:
            pickle.dump(old, handle)


def save_neutron(data, T, rho, type, new, Newfile):
        file = 'data\\data_neutron.npy'
        
        if Newfile:
            old = {}#{'' : 0}
            with open(file, 'wb') as handle:
                pickle.dump(old, handle)
            return
        TT = str(int(T))
        RHO = str(rho)
        '''
        #PARALLELIZED
        if rho==rhos[0] and new:
            with open(file, 'rb') as handle:
                old = pickle.load(handle)
            old[TT] = {}
            old[TT][RHO] = {}
            old[TT][RHO][type]=data###

        #elif rho==rhos[0] and T!=5 and new:
        #    with open(file, 'rb') as handle:
        #        old = pickle.load(handle)
        #    old[TT] = {}
        #    old[TT][RHO] = {}
        #    old[TT][RHO][type]=data

        elif rho!=rhos[0] and new:
            with open(file, 'rb') as handle:
                old = pickle.load(handle)
            old[TT][RHO] = {}
            old[TT][RHO][type]=data#{}

        else:
            with open(file, 'rb') as handle:
                old = pickle.load(handle)
            old[TT][RHO][type]=data
        #'''
        #'''
        #NOT PARALLELIZED
        if rho==rhos[0] and T == 5 and new:
            old = {}
            old[TT] = {}
            old[TT][RHO] = {}
            old[TT][RHO][type]=data

        elif rho==rhos[0] and T!=5 and new:
            with open(file, 'rb') as handle:
                old = pickle.load(handle)
            old[TT] = {}
            old[TT][RHO] = {}
            old[TT][RHO][type]=data

        elif rho!=rhos[0] and new:
            with open(file, 'rb') as handle:
                old = pickle.load(handle)
            old[TT][RHO] = {}
            old[TT][RHO][type]=data#{}

        else:
            with open(file, 'rb') as handle:
                old = pickle.load(handle)
            old[TT][RHO][type]=data#{}
        #'''
        with open(file, 'wb') as handle:
            pickle.dump(old, handle)
    
        
#, protocol=pickle.HIGHEST_PROTOCOL)
    

'''
with open(npy_name, 'wb') as handle:
    pickle.dump(old, handle, protocol=pickle.HIGHEST_PROTOCOL)

# Load data (deserialize)
with open('filename.pickle', 'rb') as handle:
    unserialized_data = pickle.load(handle)
    '''