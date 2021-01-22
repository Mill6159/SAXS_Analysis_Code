## Description ##
# This script runs Crysol
# i.e. Models experimental SAXS data given a known high resolution structure
# NOTE: The calculations can be done with just a high resolution structure NOT fit to experimental data
# It then reports the determined Radius of Gyration (Rg) & if an experimental profile is provided
# the script will report the calculated chi-squared

## Imports

from PlotClass import *
from FileParser import *
from SAXS_Calcs import *
import os
from itertools import islice

## Classes

figs = PlotClass(notify=False)
parser = FileParser(notify=False)
calcs = SAXSCalcs()

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
  justPDB_Mode=False
elif user_Just_PDB in falseList:
  justPDB_Mode = True
else:
  print('*'*4,'Uninterpretable response.. Terminating script')
  print('#'*50)
  quit()



if justPDB_Mode == True:
  print('Running just PDB mode')





else:
  user_exptFile = input('Experimental Profile (.dat): ')
  fileList=user_exptFile
  user_PDBFile = input('PDB File: ')

  Q, IQ, ERROR = calcs.runCrysol(datFileList=fileList,
                  PDB_File=user_PDBFile,
                  qmax=0.3)

  # dat_files = parser.load_datFiles(fileList=fileList)

  # # Crysol Modeling

  # for i in dat_files:
  #     n=len(dat_files[str(i)]['ERROR']) # scale the size of the profile for downstream residuals calculations
  #     cmd = 'crysol %s %s -lm 50 -sm 0.3 -kp=ii -p %s'%(str(user_PDBFile),str(fileList),str(fileList))
  #     proc = subprocess.Popen([cmd], stdout=subprocess.PIPE, shell=True)
  #     (out, err) = proc.communicate()

  # crysol_List=[]
  # crysol_List.append(fileList + '.log')

      
  # crysol_fits = []
  # crysol_fits.append(fileList+'.fit')

  # crysol_files = parser.load_fitFiles(fileList=crysol_fits)
      
  # for fn,k in zip(crysol_List,crysol_files):
  #     crysol_parse=[]
  #     crysol_parse2=[]
  #     with open(fn, "r") as f: # important to use with command when dealing with files
  #         counter = 0
  #         print('File: %s' %str(fileList))
  #         c=0
  #         for line in f:
  #             c+=1
  #             if 'Chi^2' in line:
  #                 crysol_parse.append(''.join(islice(f, 1)))
  #             if 'Rg from the slope of net intensity' in line:
  #                 print(line)
  #                 crysol_parse2.append(line)
      
  #     Rg=[i for j in crysol_parse2[0].split() for i in (j, ' ')][:-1][18]
  #     print('Final Rg (A): ', Rg)
  #     chi_sq=[i for j in crysol_parse[0].split() for i in (j, ' ')][:-1][18]
  #     print('Chi-Squared: ',chi_sq)

      
  #     pairList=[[crysol_files[str(k)]['Q'],crysol_files[str(k)]['I(Q)']],
  #               [crysol_files[str(k)]['Q'],crysol_files[str(k)]['I(Q)_mod']]]
      

  #     labelList=['%s'%str(fileList),'Model']

  #     colorList=['#9B9B9B','#00B2AF']

  # #     figs.nPlot_variX_and_Color(pairList=pairList,labelList=labelList,colorList=colorList,
  # #                               xlabel='q ($\\AA^{-1}$)',
  # #                               ylabel='Intensity, I(Q)',
  # #                               savelabel='%s_crysol'%str(z))

  #     X1=crysol_files[str(k)]['Q']
  #     Y1=crysol_files[str(k)]['I(Q)']
  #     X2=crysol_files[str(k)]['Q']
  #     Y2=crysol_files[str(k)]['I(Q)_mod']
  #     Y1err=[1] * len(X1)
      
  #     skip=next((i for i, x in enumerate(Y1) if x), None) # x!= 0 for strict match

  #     figs.vertical_stackPlot(X1=X1[skip:],Y1=Y1[skip:],Y1err=Y1err[skip:],X2=X2[skip:],Y2=Y2[skip:],ylabel1='Intensity, I(Q)',
  #                             ylabel2='Residuals',xlabel='Q ($\mathring{A}^{-1}$)',
  #                             Label1='%s'%fileList, Label2='Crysol Model',saveLabel='%s_Crysol_fit_EM10mer'%str(fileList),
  #                             bottomPlot_yLabel='$ln(\\frac{I_{expt}(q)}{I_{model}(q)}) \cdot (\\frac{1}{\sigma_{expt}})$',
  #                             LinLin=False,
  #                             linewidth=4,labelSize=20,darkmode=False,plot=True)










