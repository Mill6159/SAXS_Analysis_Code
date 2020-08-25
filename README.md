# Basic to semi-advanced SAXS analysis #

## _Welcome_  ##
### _This repository contains a series of scripts used for processing SAXS data. Each script/class should have a description within the script as well as here in the README file. Note, these functions are continually updated and improved for bug fixes. Therefore, if issues arise, please post them to the GitHub issues page here or contact Robert Miller (rcm347@cornell.edu). Also, note that because they are updated the README.md file may become slightly outdated, but will be updated entirely on a semi-regular basis._ ###

### _FileParser.py_ ###

A script that builds the class FileParser(). This class reads in file types from subtracted SAXS profiles (.dat)
to GNOM and DATGNOM output files (.out). I will describe each of the functions below:

**(1)** 

```python
loadOutFile(filename)
```
This function reads the output from the runGNOM() or runDatgnom() functions available in Basic_SAXS_Calcs.py.
In fact, it just reads the output file (which is a bit messy) from GNOM and DATGNOM and does some preliminary processing
for downstream use. If the .out file is not in the current working directory, the filepath must be specified along with the file name.

### _PlotClass.py_ ###

A script that builds the class PlotClass(). This class contains many functions for various different plots
all built with matplotlib. I will describe each of the functions below:

**(1)**
```python
basicPlot(X,Y,plotlabel='',savelabel='',xlabel='',ylabel='NOT PROVIDED')
``` 
This function takes in an X/Y pair, these of course must be of equal length. 
It also requires a plotlabel, which is placed into the legend.
savelabel is the label of the .png file (i.e. savelabel.png) output into the current working directory. xlabel/ylabel are self explanatory.

**(2)** 

```python
semilogyPlot(X,Y,plotlabel='',savelabel='',xlabel='',ylabel='',linewidth=4)
```

Just like the basicPlot() function but plots the X/Y pair on a Log-linear scale.

### _Basic_SAXS_Calcs.py_ ###

A script that actually performs the SAXS calculations (P(r), Guiner, etc) and generates the class
BasicSAXS().

**(1)** 
```python
lineModel(x,m,b)
```

This function is used in other functions and simply defines a line.


**(2)**
```python
lsq_w_sigma(X, Y, SIG)
```
This is a special version of linear least squares that uses 
known sigma values of Y to calculate standard error in slope and intercept. 
The actual fit is also SIGMA-weighted

Outputs: slope(my), intercept(by), sigma slope(smy), sigma intercept (sby)


### _FoxS.py_ ###

A script that runs FoxS calculations using a local installation of the FoxS
profile generator 
url: https://modbase.compbio.ucsf.edu/foxs/download.html

**(1)**
```python
FoxS_simple(pdb1=pdb1,expt1=expt1,fast_mode=fast_mode,nq=nq,maxq=maxq, exH=exH,offset=offset,
	plot=True)
```


## I hope these are useful to more people than myself ##
## Best of luck ##