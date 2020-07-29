# Basic to semi-advanced SAXS analysis #

## _Welcome_  ##
### _This repository contains a series of scripts used for processing SAXS data. Each script/class should have a description within the script as well as here in the README file. Note, these functions are continually updated and improved for bug fixes. Therefore, if issues arise, please post them to the GitHub issues page here or contact Robert Miller (rcm347@cornell.edu). Also, note that because they are updated the README.md file may become slightly outdated, but will be updated entirely on a semi-regular basis._ ###

### _FileParser.py_ ###

A script that builds the class FileParser(). This class reads in file types from subtracted SAXS profiles (.dat)
to GNOM and DATGNOM output files (.out). I will describe each of the functions below:

**(1)** loadOutFile()

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

### _Basic_SAXS_Calcs.py_ ###
