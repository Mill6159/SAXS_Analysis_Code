# SAXS Analysis Pipeline
The automated pipeline for analysis of SAXS data

Requires the following classes located in the root respository:  
*(1)* Basic_SAXS_Calcs.py  
*(2)* FileParser.py  
*(3)* PlotClass.py  
*(4)* SAXS_Calcs.py  

## General flow of data analysis

**(1)** Import .dat files  
* Should contain three columns with any non-numeric data contained on lines starting with #  
  

    Q     I(Q)    Error  

* Place all .dat files into the sub-directory *SAXS_Pipeline/dat_files/*  
	* The script will then read in every single file in that directory that contains the .dat suffix  


