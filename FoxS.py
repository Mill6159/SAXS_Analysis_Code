########################################################
# This function runs FoxS from "imp" as installed from:
# https://integrativemodeling.org/download-mac.html
#
# I mainly developed this so I could have more control
# of how many input structures there were! But..
# now we'll have much better plots and more control
# in general.
########################################################


########################################################
# Imports
import os 
import sys

########################################################

imp_dir='/usr/local/Cellar/imp/2.13.0_2/bin' # build an automated way to find this
# OR
imp_dir='/usr/local/bin/foxs'


class FoxS:

	def __init__(self, imp_dir):
		''' 
		Explain class here
		'''
		self.imp_dir=imp_dir

	def runFoxS(self, pdb_list,expt_prof_list,ns,maxq):
		'''
		explain function here.
		pdb_list:
		expt_prof_list:
		ns:
		maxq:
		'''



