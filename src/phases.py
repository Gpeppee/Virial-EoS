import numpy as np
import src.numericals as numericals
from src.constants import *
import os

### I did not check this
def READ_PHASE_SHIFTS(isospin_pairs):
    base = os.getcwd()
    #!folder="./data_phase_shifts/"//isospin_pairs//"/"
    folder='data_phase_shifts/np/'
    if (isospin_pairs == 'np'):
        itz1=-1
        itz2=1
    elif (isospin_pairs == 'nn'):
        itz1=-1
        itz2=-1
    hhhh = 0
    if( read_ps ):
        #! ... THERE ARE 6 AVAILABLE PARTIAL WAVES PER VALUE OF J
        #! ... BUT WE ONLY CONSIDER THE 4 DIAGONAL CHANNELS
        Ddim = 728
        phase_shifts = np.zeros((Ddim,J_max+1,Nch))
        global energy       
        energy = np.zeros((Ddim))
        for j in range(J_max+1):
            for ich in range (Nch):
                ich +=1
                ss, l1,l2,cc = PARTIAL_WAVE_ANGMOM(j,ich)
                # ! ... ANTISIMMETRY CONSTRAINTS FOR T=1
                ittot=itz1+itz2
                if( isospin_pairs == 'nn' ): #np.abs(ittot) == 2):
                   isum=1+l1+ss
                   if((isum%2)==0): 
                       ss=2

                if( ss != 2):
                    pwname=PARTIAL_WAVE_NAME(j,ich,l1,l2,ss,cc)
                    hhhh+=1
                else:
                    pwname = 'nope'
            #! ... FIND ISOSPIN
                if(j==0):
                    tt=1
                else:
                    #! ... If the channel is the singlet or the coupled triplet,
                    #! ... the isospin is given by J+1 modulus 2
                    tt = ((ss+l1 % 2) + 1) % 2 
                    #! ... SAFE CHECK
                    icheck=(l1+ss+tt)%2
                    if(icheck==0): 
                        print('Error in the pws?')#write(*,*) 'error in pws?'
                        break
                    
                #! ... READING BEGINS HERE
                #! ... WE DO NOT ASSUME ANYTHING ABOUT FILE LENGTH
                filename_er=base +'/'+ folder + pwname + '.dat'
                try:
                    data = np.loadtxt(filename_er)
                    energy[:] = data[:,0]
                    print('YOOO', pwname)
                except:
                    data = np.zeros(data.shape)            
                nlines = np.shape(data[:,0])[0]
                #if hhhh == 1: energy[:] = data[:,0]#np.reshape(data[:,0],(nlines))
                if ss != 2:
                    phase_shifts[:, j, ich-1] = data[:,1]
 
        #    endif ! GOOD PARTIAL WAVE, ss \= 2

            #enddo ! LOOP OVER CHANNELS
        #enddo ! LOOP OVER J

        #! ... THE PHASE SHIFTS HAVE BEEN OBTAINED HERE
        #! ... COMPUTE THEIR SUM AND STORE IN FILE
        #exit()
        Nener=nlines
        phase_shift_total= np.zeros(nlines)
        for iener in range(Nener):
            ps=0
            for J in range(J_max+1):# J=0,Jmax
              for ich in range(Nch):#=1,4
                
                ittot=itz1+itz2
                if( np.abs(ittot) == 2):
                    ps = ps + float(2*J+1)*phase_shifts[iener,J,ich]
                elif ( np.abs(ittot) == 0):
                    #!calculate the isospin
                    if ((ich==0) or (ich==2) or (ich==3)):
                        tt=((J%2)+1)%2
                    elif (ich==1):
                        tt=((((J+1)%2)+1)%2)
            	    #endif
                    
                    ps = ps + float(2*J+1)*float(2*tt+1)*phase_shifts[iener,J,ich]
               #endif
             #enddo
            #enddo
            #!        write(*,*)
            phase_shift_total[iener]=ps
        #enddo
        #folder="output\\"
        #filename=folder+"pstot_J"+str(J_max)+"_"+isospin_pairs+".dat"

    else: #! IF WE DO NOT READ, WE COMPUTE THE PS FROM NP OR NN SCATTERING LENGTH/EFF RANGE
        folder="output/"
        filename=folder//"ps_a"//isospin_pairs//".dat"
        Nener=700
        energy = np.zeros(Nener)
        waux   = np.zeros(Nener)
        phase_shift_total = np.zeros(Nener)
        #if( .not.ALLOCATED(energy) ) ALLOCATE( energy(Nener) )
        #if( .not.ALLOCATED(waux) ) ALLOCATE( waux(Nener) )
        #if( .not.ALLOCATED(phase_shift_total) ) ALLOCATE( phase_shift_total(Nener) )
        
        e_in=0
        e_fi=350
        energy, waux = np.gauss(e_in,e_fi,Nener)
               
    #endif

    #! ... OPEN OUTPUT PHASESHIFT FILE
    #open(11,file=trim(filename))
    #! ... WRITE FILE HEADER
    #if( not read_ps): write(11,211) isospin_pairs
        #  write(11,111)

        #! ... OUTPUT DATA
    #'''
    for iener in range(Nener):
            if( read_ps ):
                if( isospin_pairs == "nn"):
                    # !  CIB CORRECTION FOR 1S0 WAVE ONLY
                    pss= PHASE_SHIFT_SWAVE(energy[iener], a_nn,reff_nn,xmassn,xmassn) - PHASE_SHIFT_SWAVE(energy[iener], a_np,reff_np,xmassn,xmassp)
                    phase_shift_total[iener] += pss
                #endif
            else:
                if( isospin_pairs == "np"):
                    phase_shift_total[iener]=PHASE_SHIFT_SWAVE(energy(iener), a_np,reff_np,xmassn,xmassp)
                elif( isospin_pairs == "nn"):
                    phase_shift_total[iener]=PHASE_SHIFT_SWAVE(energy(iener), a_nn,reff_nn,xmassn,xmassn)
               
    #'''
  #END SUBROUTINE
    #print(hhhh)
    return energy, phase_shift_total                





def PHASE_SHIFT_SWAVE(E,a,reff,xm1,xm2):
    xm=(xm1*xm2)/2#(xm1+xm2)
    k=np.sqrt(E*xmassn)/hbc

    if( reff != 0 ):
        #    ! ... NP SCATTERING LENGTH + EFFECTIVE RANGE FORMULA
        ps= np.arctan(k*a/(-1 + k**2*a*reff/2))#*180/pi
    else:
        #    ! ... NP SCATTERING LENGTH ONLY
        #    ! ... arccot(-1/x)=arctan(x)
        ps=-np.arctan(k*a)#*180/pi
    
    return ps



def  PARTIAL_WAVE_NAME(J,ich,l1,l2,ss,cc):
    #LPW="SPDFGHIKLM"
    LPW="spdfghiklm"
    JPW="0123456789"
    SPW="13"
    if ss != 2:
        SPW = str(int(2*ss+1))

        if(J < 9): 
            Ja=JPW[J]
        else:
            Ja=str(int(J))

        if(not cc):
            if(l1<9):
                La=LPW[l1]
            else:
                La=str(int(l1))
            return SPW + La + Ja
        elif (cc):
            if(l1<9):
              La=LPW[J-1]     #####check the indexes
              Lb=LPW[J+1]     #####check the indexes
            else:
              La = str(int(J-1))
              Lb = str(int(J+1))
              Lc = str(int(J-2))
              Ld = str(int(J+2))
            

            #!PARTIAL_WAVE_NAME=SPW//trim(La)//trim(Lb)//trim(Ja)
            if( ich == 3 ): return SPW + La + Ja
            if( ich == 4 ): return SPW + Lb + Ja
            #if( ich == 5 and J>2 ): PARTIAL_WAVE_NAME=SPW + Lc + Ja
            #if( ich == 6 and J>2 ): 
            #    PARTIAL_WAVE_NAME=SPW + Ld + Ja
            #else:
            #    PARTIAL_WAVE_NAME=""
    else:
        return ""
      

    return PARTIAL_WAVE_NAME


def PARTIAL_WAVE_ANGMOM(jc, lsj):
    s=2
    if(jc == 0):            # ...  J0
        coupled= False
        
        if(lsj == 1):       # .... 1S0
            l1=0
            l2=0
            s=0
            return  s, l1, l2, coupled
        elif(lsj == 4):     # .... 1P0

            l1=1
            l2=1
            s=1
            return  s, l1, l2, coupled
        else:
            return  s, False, False, coupled   ################
        
    else:
  
        if(lsj==1):         # ... Singlet
            coupled=False
            l1=jc
            l2=jc
            s=0
            
  
        elif(lsj==2):       # ... Uncoupled triplet
            coupled=False
            l1=jc
            l2=jc
            s=1
  
        elif(lsj==3):       # ... Coupled triplet - V--
            coupled=True
            l1=jc-1
            l2=jc-1
            s=1

        elif(lsj==4):       # ... Coupled triplet - V++
            coupled=True
            l1=jc+1
            l2=jc+1
            s=1
  
        elif(lsj==5):       # ... Coupled triplet - V-+
            coupled=True
            l1=jc-1
            l2=jc+1
            s=1
  
        elif(lsj==6):       # ... Coupled triplet - V+-
            coupled=True
            l1=jc+1
            l2=jc-1
            s=1

    return  s, l1, l2, coupled
