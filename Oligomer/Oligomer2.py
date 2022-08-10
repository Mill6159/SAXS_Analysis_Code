#!/usr/bin/env python3

############################################
# Describe script here
# Runs ATSAS OLIGOMER software
# Url:
# & then generates a nice quality figure
# from the output
#
# It also makes the input "interactive"

############################################################
# Module Inputs
import sys
import os
import subprocess
from subprocess import PIPE, Popen
import re
from termcolor import colored,cprint
# import numpy as np
import matplotlib as mpl
from matplotlib import pyplot as plt

############################################################

############################################################
# User Inputs


############################################################

############################################################
# Figure Inputs/Building

def nPlot_variColor(pairList,labelList,colorList,savelabel,xlabel='No Label Provided',ylabel='No Label Provided',
              LogLin=True,LinLin=False,LogLog=False,linewidth=3,
              set_ylim=False,ylow=0.0001,yhigh=1,darkmode=False):
        '''
        :param pairList: list of lists (tuple), must be [[x1,y1],...[xn,yn]]
        :param labelList: list of length n, labeling the sets of tuples in pairList
        :param savelabel:
        :param xlabel:
        :param ylabel:
        :param linewidth:
        :return:
        '''

        if darkmode==True:
            plt.style.use('dark_background')
        else:
            mpl.rcParams.update(mpl.rcParamsDefault)

        fig=plt.figure(figsize=(10,8)) # set figure dimensions
        ax1=fig.add_subplot(1,1,1) # allows us to build more complex plots
        for tick in ax1.xaxis.get_major_ticks():
            tick.label1.set_fontsize(20) # scale for publication needs
            tick.label1.set_fontname('Helvetica')
        for tick in ax1.yaxis.get_major_ticks():
            tick.label1.set_fontsize(20) # scale for publication needs
            tick.label1.set_fontname('Helvetica')

        if LogLin == True and LinLin == True and LogLog == True: # kicks you out of the function if you set more then one mode to true
            print("Cannot set more than one mode equal to True")
            return
        elif LogLin == True and LinLin == True:
            print("Cannot set more than one mode equal to True")
            return
        elif LogLin == True and LogLog == True:
            print("Cannot set more than one mode equal to True")
            return
        elif LinLin == True and LogLog == True:
            print("Cannot set more than one mode equal to True")
            return

        n=0
        if LogLin==True:
            for i,j in zip(pairList,colorList):
                plt.semilogy(i[0],i[1],
                            label=labelList[n],
                            color=j,
                            linewidth=linewidth,linestyle='dashed')
                n+=1
        elif LinLin==True:
            for i in pairList:
                plt.plot(i[0],i[1],
                            label=labelList[n],
                            linewidth=linewidth,linestyle='dashed')
                n+=1
        elif LogLog==True:
            for i in pairList:
                plt.plot(i[0],i[1],
                            label=labelList[n],
                            linewidth=linewidth,linestyle='dashed')
                n+=1
                ax1.set_yscale('log')
                ax1.set_xscale('log')


        plt.ylabel(ylabel,size=22)
        plt.xlabel(xlabel,size=22)
        plt.legend(numpoints=1,fontsize=18,loc='best')


        if set_ylim==True:
            ax1.set_ylim(ylow,yhigh)
        else:
            print('Using default y-limit range for the plot: %s'%savelabel)
        fig.tight_layout()

        plt.savefig(savelabel+'.png',format='png',bbox_inches='tight',dpi=300)
        plt.show()


############################################################
# Build a few functions to disable and enable print output
# to the terminal window

# Disable
def blockPrint():
    sys.stdout = open(os.devnull, 'w')

# Restore
def enablePrint():
    sys.stdout = sys.__stdout__


############################################################


cprint('-----> Beginning Oligomer calculations','green',attrs=['bold'])
cprint('Note: Assumes ATSAS is installed in $USER directory \n& you also (of course) have access to FFMAKER','green',attrs=['bold'])


############################################################

def cmdline(command):
    process = Popen(
        args=command,
        stdout=PIPE,
        shell=True
    )
    return process.communicate()[0]


proc = subprocess.Popen(["echo $USER"], stdout=subprocess.PIPE, shell=True)
(out, err) = proc.communicate()
localUser = str(out.decode('ascii')).rstrip() # rstrip() removes line break

############################################################


OligomerDir= '/Users/' + localUser + '/ATSAS-3.0.1-1/bin/Oligomer'
# print(OligomerDir)

if os.path.exists(OligomerDir):
	cprint('--> Local Oligomer install located','green',attrs=['bold'])
else:
	cprint('--> Local Oligomer install was NOT able to be located..\n--> This script assumes the install is located at: %s'%OligomerDir,'green',attrs=['bold'])
	cprint('--> Terminating Oligomer calculations','green',attrs=['bold'])
	sys.exit()

############################################################
# Let user input the PDB files themselves

pdbNo_text=colored("How many different structures as input? ",'cyan',
					  attrs=['bold'])

pdbNo = int(input(pdbNo_text))

pdbList = []
for i in range(pdbNo):
	pdbinput=colored("PDB File #%i: "%(i+1),'cyan',
				  attrs=['bold'])
	pdbList.append(input(pdbinput))
	# pdbList[i] = (str(input(pdbinput)))

# print(pdbList)	

############################################################
############################################################
# Must build the form factor file

ff_cmd = 'ffmaker'

for i in pdbList:
	ff_cmd = ff_cmd + ' ' + i

tmp_lst=[]
for i in pdbList:
	# print(i)
	tmp_lst.append(i.rsplit('.',1)[0])


# print(tmp_lst)

outfile='FF_'
for i in tmp_lst:
	outfile = outfile + i + '_'
outfile = outfile + '.dat'

# print(outfile)

ff_cmd = ff_cmd + ' ' + '-out ' + outfile

qmax_text=colored("Max q-value (A)? (Default=0.3): ",'cyan',
				  attrs=['bold'])

qmax = input(qmax_text)
if qmax == '':
	qmax = str(0.3)

ff_cmd = ff_cmd	+ ' ' + '-smax=' + qmax + ' ' + '-lmmax=32'

# print(ff_cmd)
cprint(ff_cmd,'yellow')

proc = subprocess.Popen([ff_cmd], stdout=subprocess.PIPE, shell=True)
(out, err) = proc.communicate()

############################################################
############################################################
expt_txt=colored("Expt File: ",'cyan',
				  attrs=['bold'])
expt = input(expt_txt)
expt_tmp=expt.rsplit('.',1)[0]
# print('regular expression outfile:', expt_tmp)

out_tmp = outfile.rsplit('.',1)[0]

ol_cmd = 'Oligomer -ff ' + outfile + ' -dat ' + expt + ' -smax=' + qmax + \
' -out=' + 'Oligo_' + out_tmp + '.log' + ' -fit=' + 'Oligo_' + out_tmp + '.fit' + ' --ws'

# print(ol_cmd)

fit_file = 'Oligo_' + out_tmp + '.fit'
# print(fit_file)

cprint(ol_cmd,'yellow')

proc = subprocess.Popen([ol_cmd], stdout=subprocess.PIPE, shell=True)
(out, err) = proc.communicate()

############################################################
############################################################
### NEED TO READ IN THE .log FILE!! RM

N=5
# Method 1:
out2=[]
filename='Oligo_' + out_tmp + '.log'
with open(filename, "r") as f: # important to use with command when dealing with files
    counter = 0
    # print('File: %s' % filename)
    for line in f:
    	counter += 1
    	out2.append(line)
        # print (line)


# print(out2)
n=len(expt_tmp)
# pattern_chi2 = re.compile(r'C[a-zA-Z]+\^2...([0-9]+).([)0-9]+)?')
pattern_chi2 = re.compile(r'%s(..)([0-9]+)\.([0-9]{2})'%expt_tmp) # not working... 
# pattern_chi2 = re.compile(r'%s'%expt_tmp) # not working... 
# pattern_chi2 = re.compile(r'Number of oligomers')
# match_chi2 = pattern_chi2.finditer(out2.decode('ascii')) # reading out
match_chi2 = pattern_chi2.finditer(str(out2)) 

for match in match_chi2:
    il,ij=match.span()[0],match.span()[1]
    chi2=str(out2)[il+n+2:ij]
    print('chi-squared:',chi2) # return Chi^2 value from fit
    # approx_chi2=match.group(1)

os.system('more %s'%filename)
############################################################
############################################################


data=np.loadtxt(fit_file,
	dtype={'names': ('q', 'I(q)','I(q)_Error','I(q)_Fit'), 'formats': (np.float,np.float,np.float,np.float)},skiprows=1)


pairList=[[data['q'],data['I(q)']],[data['q'],data['I(q)_Fit']]]
labelList=['Raw Data','Model']
colorList=['#D6610B','#01A9BD']
savelabel='Oligo_' + out_tmp
xlabel='q ($\\AA^{-1}$)'
ylabel='I(q)'

blockPrint()
nPlot_variColor(pairList=pairList,labelList=labelList,colorList=colorList,savelabel=savelabel,xlabel=xlabel,ylabel=ylabel,
              LogLin=True,LinLin=False,LogLog=False,linewidth=3,
              set_ylim=False,ylow=0.0001,yhigh=1,darkmode=False)
enablePrint()




