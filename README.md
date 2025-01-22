# Virial-EoS

This code calculates thermodynamic properties of the matter in the Virial approximation as in CITATION. 

## Requirement
The code is in Python 3.11.9 and uses the following packages that may need installation:
 - joblib  ( parallel calculations )
 - scipy   ( scientific computations )
 - mpmath  ( for calculation in extreme relativistic case using polylogarithm)
 - numpy   ( standard library for handling variables )
 - h5py / pickle  ( print and read the data )

Alert: The coding has been done on windows so the user may need to change the format of the paths for reading and printing files

## Running
To run the code: 
    -   python main.py
It will ask you whether to compute $\beta$-equilibrium. If no, it will ask you whether to fix a specific proton fraction. If no, it will run automatically what is written in the constant.py file as xp0 and xp1

## Output
The output is divided in different files depending from the temperature. Every file is built in hierarchical format: every entry runs  
