# Main Script #

import os

class SAXS:

	'''
	This should take in a "SAXS Object"
	which should be a list of paths to .dat files
	'''

	def __init__(self,
		file_directory=None,
		ATSAS_directory=None,
		Rg=0.0, # ! maybe remove as input variable
		MW = 0.0,
		IO = 0.0): # ! maybe remove as input variable
		'''
		Doc String
		Class "SAXS"

		This is my MASTER class, i.e. it controls everything else.

		Basic idea:


		Inputs:

		Functions:


		'''
		# TODO: should make this a super class (i.e. the parent)


		print('#'*51)
		print('---> Class "SAXS" has been called')
		print('#'*51)

		# ! deal with file(s) location

		if file_directory == None:
			print('#'*51)
			print('File directory NOT provided ... Terminating script.')
			print('#'*51)
			quit()
		else:
			self.file_directory = file_directory
			print('.dat files from the following directory will be parsed: ',file_directory)

		# ! Find ATSAS directory for GNOM, etc

		if ATSAS_directory == None:
			print('lets find it...')
		else:
			self.ATSAS_directory = ATSAS_directory

		# ! Declaring class variables
		# - Need to do this on a per file basis
		# - i.e. if you import 10 files there is a dictionary key (self.Rg['DICTIONARY KEY']) that you can reference.

		# TODO: Make a dictionary of Rg values per input file
		# - Also from the dictionary generate a list of the values for plotting, etc later on.
		

		self.Rg = {}

		def RgList(self):
			'''
			Create a list of Rg values from the dictionary
			'''
			# ! How to deal with this if dictionary is empty

			RgList = []

			if bool(self.Rg) == True: # only run if dictionary is not empty
				for key,value in self.Rg:
					RgList.append(value)

				self.RgList = RgList # angstrom
			else:
				print('Dictionary is empty.. No list created.')

			return RgList 

		print(bool(self.Rg))

		self.MW = MW # kDa
		self.IO = IO # should be absolute scale - match units of input files






test = SAXS(file_directory = '../datFiles/')