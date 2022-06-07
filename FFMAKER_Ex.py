		def run_FFMAKER(self,
										FFMaker_outname='FFout'):
			'''
			Another simple function.


			All experimental data files with "pdb" extension from the working directory 
			will be taken and CRYSOL will run for these files with the parameters:

			-smax=0.3, -ns=201 and -lmmax=20

			for FFMAKER command:
			-smax sets 0.3 A^-1 as max 'q' value
			-ns sets 256 points for the Crysol curve (max)
			-lmmax sets 50 harmonics (max)

			The formfactors will be taken from the corresponding intensity values
			'''

			ffmaker_cmd = 'FFMAKER '
			for i in self.PDBList:
				ffmaker_cmd = ffmaker_cmd + (str(i) + ' ')

			outname = '%s.dat'%str(FFMaker_outname)
			ffmaker_cmd = ffmaker_cmd + '-out %s -smax 0.3 -ns 256 -lmmax 50'%str(outname)
			print('Running FFMAKER . . . .') # ! Make a live updater b/c this is kinda slow
			print('Output filename: %s'%str(outname))
			proc = subprocess.Popen([ffmaker_cmd], 
															stdout=subprocess.PIPE, 
															shell=True)
			(out, err) = proc.communicate()
			print('FFMAKER finished')
			print('#'*50)
