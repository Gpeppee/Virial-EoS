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
It will ask you whether to compute $\beta$-equilibrium. If no, it will ask you whether to fix a specific proton fraction. If no, it will run automatically what is written in the constant.py file as xp0 and xp1.

After running the calculation make sure to run
    -   python organize.py
To organize the data in a more memory efficient way (still to improve).

## Output
The output is divided in different files depending from the temperature. Every file is built in hierarchical format to be read by python as a dictionary so all the quantities will be clearly stated. 

All the files will be printed in a folder named 'data', so either you change this in the constants.py file or you create an empty folder named 'data'.

## Refernces
To know more about the approach and the formalism go to Rivieccio et al. 2025 (arXiv:2501.16795). 

Have fun!
