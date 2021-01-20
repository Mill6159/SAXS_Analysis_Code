## Description ##
# This script runs Crysol
# i.e. Models experimental SAXS data given a known high resolution structure
# NOTE: The calculations can be done with just a high resolution structure NOT fit to experimental data
# It then reports the determined Radius of Gyration (Rg) & if an experimental profile is provided
# the script will report the calculated chi-squared

## Imports

from PlotClass import *
import os

## Terminal Header

print('#'*50)
print('#'*13,'Crysol Fitting Package','#'*13)
print('#'*50)
print('#'*5, 'Package requires the following modules', '#'*5)
print('NumPy','\nmatplotlib','\nitertools','\nSAXS_Calcs.py - Written by Rob Miller \n->(https://github.com/Mill6159/SAXS_Analysis_Code)')
# print('')
print('#'*50)
print('#'*15,'Default Parameters','#'*15)
print('#'*50)
print('Maximum number of harmonics: ', 50)
print('Mode #2: Prediction Mode')
print('#'*50)
print('#'*18,'User Inputs', '#'*19)
print('#'*50)

## User Inputs

user_Just_PDB = input('Are you inputing a high-res structure & \nexperimental data? (y/n): ')

trueList = ('y','Y','yes','Yes','YES', '',' ')
falseList = ('n','N','No','NO')

if user_Just_PDB in trueList:
  justPDB_Mode=True
elif user_Just_PDB in falseList:
  justPDB_Mode = False
else:
  print('*'*4,'Uninterpretable response.. Terminating script')
  print('#'*50)
  quit()




