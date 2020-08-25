############################################################
# Author: Robert Miller (rcm347@cornell.edu)
# 

# --->
# This script runs a local download of FoxS
# So, of course.. You must first have FoxS installed
# https://modbase.compbio.ucsf.edu/foxs/download.html
# Make sure it installs to /usr/local/bin/foxs
# which should be the default

# --->
# This script also requires the two scripts 
# (1) PlotClass.py and (2) SAXSCalcs.py
# They must be located in the same directory as the FoxS.py
# script and the expt/pdb files!



############################################################
# Import all of the dependencies

from PlotClass import *
import os
import subprocess
import re

############################################################
# Call the necessary classes()
visuals=PlotClass(notify=False) # This class requires SAXSCalcs.py


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
# Check if the user has FoxS installed

print('-----> Beginning FoxS calculations')

foxsDir='/usr/local/bin/foxs'

if os.path.exists(foxsDir):
	print('--> Local FoxS install located')
else:
	print('--> Local FoxS install was NOT able to be located..\n--> This script assumes the install is located at: %s'%foxsDir)
	print('--> Terminating FoxS calculations')
	sys.exit()

############################################################
# Let user input the files themselves
# make this more sophisticated for n inputs

pdb1 = input("PDB File #1: ")
# print(pdb1)
expt1 = input("Expt File #1: ")
# print(expt1)

fast_mode=input('Run fast mode? (default=False), (Y/N): ')
if fast_mode=='Y' or fast_mode=='y':
	fast_mode='True'
elif fast_mode=='N' or fast_mode=='n':
	fast_mode='False'
elif fast_mode=='':
	fast_mode='False'
# print(fast_mode)

nq=input('Number of points in the profile? (default=500), (Y/#): ')
if nq=='Y' or nq=='y':
	nq='True' # true means use default of 500
elif nq=='':
	nq='True'
else:
	nq=nq
# print(nq)

maxq=input('Maximum q-value (%s). (default=length of dataframe=Y), (Y/#): '%str(u'\u212B'))
if maxq=='Y' or maxq=='y':
	maxq='True' # true means use default of 500
elif maxq=='':
	maxq='True'
else:
	maxq=maxq

exH=input('Explicitly consider hydrogens in PDB files? (default = false), (Y/N): ')
if exH=='Y' or exH=='y':
	exH='True' # true means use default of 500
elif exH=='':
	exH='False'
else:
	exH='False'

offset=input('Use offset in fitting? (default = false), (Y/N): ')
if offset=='Y' or offset=='y':
	offset='True' # true means use default of 500
elif offset=='':
	offset='False'
else:
	offset='False'


############################################################
# FoxS function

def FoxS_simple(pdb1=pdb1,expt1=expt1,fast_mode=fast_mode,nq=nq,maxq=maxq, exH=exH,offset=offset,
	plot=True):
	'''
	This function runs a local FoxS download
	It takes one experimental file and one PDB file
	It asks the user what modes it wants to access (almost all are available.. at least the ones we want)
	It then generates the FoxS command, passes it to your terminal, parses the output and 
	returns a few statistics on the fit as well as a beautiful plot.

	--> input:
	pdb1:
	expt1:

	--> output:


	'''

	print('\n----> Now building FoxS command....\n')

	cmd = 'foxs' + ' ' + pdb1 + ' ' + expt1 #+ ' ' + '-v True'

	if fast_mode=='True':
		cmd = cmd + ' ' + '-r True'
	else:
		cmd = cmd

	if nq=='True':
		cmd = cmd
	else:
		cmd = cmd + ' ' + '-s %s'%str(nq)

	if maxq=='True':
		cmd = cmd
	else:
		cmd = cmd + ' ' + '-q %s'%str(maxq)

	if exH=='True':
		cmd = cmd + ' ' + '-h True'
	else:
		cmd = cmd

	if offset=='True':
		cmd = cmd + ' ' + '-o True'
	else:
		cmd = cmd

	# print(cmd)
	# # run command
	# os.system(cmd)

	proc = subprocess.Popen([cmd], stdout=subprocess.PIPE, shell=True)
	(out, err) = proc.communicate()
	# print ("program output:", out.decode('ascii')) # use re module to parse the output
	print('\n')

	# ---->
	# search for chi-squared value

	pattern_chi2 = re.compile(r'C[a-zA-Z]+\^2...([0-9]+).([)0-9]+)?') # finds chi-squared value from FoxS output

	# ---->
	# search for coefficient values (c1 and c2)

	pattern_c1 = re.compile(r'c1...(-)?[0-9]+.([)0-9]+)?')
	pattern_c2 = re.compile(r'c2...(-)?([0-9]+).([)0-9]+)?')

	# += 1 or more until we hit @

	# ---->
	match_chi2 = pattern_chi2.finditer(out.decode('ascii'))
	match_c1 = pattern_c1.finditer(out.decode('ascii'))
	match_c2 = pattern_c2.finditer(out.decode('ascii'))


	# ---->

	# print('############################')
	# print('FoxS OUTPUT')
	# print(out.decode('ascii'))

	print('####################################')
	print('Parsed output from FoxS calculations')
	print('####################################')
	# ---->	

	for match in match_chi2:
	    il,ij=match.span()[0],match.span()[1]
	    chi2=out.decode('ascii')[il:ij]
	    print(chi2) # return Chi^2 value from fit
	    approx_chi2=match.group(1)

	for match in match_c1:
	    il,ij=match.span()[0],match.span()[1]
	    c1=out.decode('ascii')[il:ij]
	    print(c1) # return c1 value
	    if c1=='c1 = 0.99':
	    	print('----> Minimum c1 value used. Quality of the fit is likely poor')
	    elif c1=='c1 = 1.05':
	    	print('----> Maximum c1 value used. Quality of the fit is likely poor')

	for match in match_c2:
	    il,ij=match.span()[0],match.span()[1]
	    c2=out.decode('ascii')[il:ij]
	    print(c2) # return c2 value
	    # print(match.group(1),match.group(2))
	    if '-2' in c2:
	    	print('----> Minimum c2 value used. Quality of the fit is likely poor')
	    elif '4' in c2:
	    	print('----> Maximum c2 value used. Quality of the fit is likely poor')


	#### fetch output file for plotting

	output_name=pdb1.replace('.pdb','') + '_' + expt1.replace('.dat','') + '.fit'

	plot_data=np.loadtxt(output_name,
	dtype={'names': ('Q', 'I(Q)','ERROR','Model'), 'formats': (np.float,np.float,np.float,np.float)},
	comments='#')

	pairList=[[plot_data['Q'],plot_data['I(Q)']],[plot_data['Q'],plot_data['Model']]]
	labelList=['Expt - $\\chi^{2} \\approx$%s'%approx_chi2,'FoxS Model']
	colorList=['k','#4BC8C8']

	blockPrint() # suppressing plot function print outs.. not necessary here.
	visuals.nPlot_variColor(pairList=pairList,labelList=labelList,colorList=colorList,savelabel=output_name+'_ScatteringProfile',
		xlabel='q ($\\AA^{-1}$)',ylabel='I(q)')
	enablePrint()



FoxS_simple()


############################################################















