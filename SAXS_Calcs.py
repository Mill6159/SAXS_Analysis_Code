#!/Users/robmiller/opt/anaconda3/bin/python3
#######################################################
## Describe overall goal of the code here
#######################################################


#######################################################
# Imports
#######################################################
import numpy as np
from scipy.optimize import curve_fit as cf
import os
import sys
import csv
import subprocess
import platform
from numba import jit
from numba import types, typed
from numba.experimental import jitclass
import math
import warnings
import traceback
import fabio
# from PlotClass import *
from FileParser import *
from itertools import islice

#######################################################
#######################################################
# Working notes:

# Create a new class to read HDF5 files directly if a config file is given!
# raw_settings are important!
# deal with this next!!

#######################################################
#######################################################


class SAXSCalcs:

    def __init__(self,notify=True):
        self.notify=notify
        if self.notify==True:
            print('--------------------------------------------------------------')
            print('SAXSCalcs class was called..')
            print('--------------------------------------------------------------')

        # other stuff.. atsas?

    def superimpose(self,ref_q,ref_i,qmin,qmax,scale_list,choice='Scale'):
        """
        Find the scale and/or offset factor between a reference curve and the
        curves of interest.
        The reference curve need not be sampled at the same q-space points.

        """

        q_star = ref_q
        i_star = ref_i
        # err_star = sasm_star.err

        q_star_qrange_min,q_star_qrange_max = qmin,qmax

        q_star = q_star[q_star_qrange_min:q_star_qrange_max]
        i_star = i_star[q_star_qrange_min:q_star_qrange_max]

        for each_scale in scale_list:

            each_q = each_scale[0]
            each_i = each_scale[1]

            each_q_qrange_min,each_q_qrange_max = (0,len(each_scale[0]))

            # resample standard curve on the data q vector
            min_q_each = each_q[each_q_qrange_min]
            max_q_each = each_q[each_q_qrange_max - 1]

            min_q_idx = np.where(q_star >= min_q_each)[0][0]
            max_q_idx = np.where(q_star <= max_q_each)[0][-1]

            if np.all(q_star[min_q_idx:max_q_idx + 1] != each_q[each_q_qrange_min:each_q_qrange_max]):
                I_resamp = np.interp(q_star[min_q_idx:max_q_idx + 1],
                                     each_q[each_q_qrange_min:each_q_qrange_max],
                                     each_i[each_q_qrange_min:each_q_qrange_max])
            else:
                I_resamp = each_i[each_q_qrange_min:each_q_qrange_max]

            if not np.all(I_resamp == i_star):
                if choice == 'Scale and Offset':
                    A = np.column_stack([I_resamp,np.ones_like(I_resamp)])
                    scale,offset = np.linalg.lstsq(A,i_star[min_q_idx:max_q_idx + 1])[0]
                elif choice == 'Scale':
                    A = np.column_stack([I_resamp,np.zeros_like(I_resamp)])
                    scale,offset = np.linalg.lstsq(A,i_star[min_q_idx:max_q_idx + 1],rcond=None)[0]
                    offset = 0
                elif choice == 'Offset':
                    A = np.column_stack([np.zeros_like(I_resamp),np.ones_like(I_resamp)])
                    scale,offset = np.linalg.lstsq(A,i_star[min_q_idx:max_q_idx + 1] - I_resamp)[0]
                    scale = 1

                return scale,offset
                # each_scale.scale(scale)
                # each_scale.offset(offset)

    def quickNormalize(self,sig):
        '''
        This function just normalizes a signal from 0-1
        but DOES NOT normalize by total area or with respect to a reference signal
        '''
        n = len(sig)
        sig = np.absolute(sig)
        maxY = np.nanmax(sig)
        minY = min(sig)
        if np.isnan(minY):
            minY = 0
        Ynew = np.empty(n)
        for i in range(len(sig)):
            Ynew[i] = (sig[i] - minY) / (maxY - minY)
        return Ynew

    def load_fitFiles(self,fileList):
      '''
      Describe function here:


      '''
      c=1
      diction={}
      for i in fileList:
        diction['data_%s'%str(c)]=np.loadtxt(i, dtype={'names': ('Q', 'I(Q)','Q_mod','I(Q)_mod'), 'formats': (np.float,np.float,np.float,np.float)}, comments='#',
          skiprows=1)
        c+=1

      return diction


    def runCrysol(self, datFileList,
                  PDB_File,
                  qmax=0.3):
        '''
        Description
        '''

        # Import data files

        
        def load_datFiles(fileList):
            '''
            Describe function here:


            '''
            c=1
            diction={}
            for i in fileList:
              # if len(fileList)==1: # RM! Some issues when fileList only contains one string entry
              #   diction['data_1']=np.loadtxt(fileList, dtype={'names': ('Q', 'I(Q)','ERROR'), 'formats': (np.float,np.float,np.float)}, comments='#')
              # else:
              # print('Loading multiple .dat files....')
              # print(i)
              diction['data_%s'%str(c)]=np.loadtxt(i, dtype={'names': ('Q', 'I(Q)','ERROR'), 'formats': (np.float,np.float,np.float)}, comments='#',
                                                   skiprows=1)
              c+=1
            
            return diction

        fileList=datFileList
        dat_files = load_datFiles(fileList=fileList)

        ## Crysol Modeling

        for i,j in zip(dat_files,fileList):
            print('Starting Crysol Run')
            print('')
            n=len(dat_files[str(i)]['ERROR']) # scale the size of the profile for downstream residuals calculations
            cmd = 'crysol %s %s -lm 50 -sm %s -kp=ii -p %s'%(str(PDB_File),str(j),str(qmax),str(j))
            print('Final command line input')
            print(cmd)
            print('')
            proc = subprocess.Popen([cmd], stdout=subprocess.PIPE, shell=True)
            (out, err) = proc.communicate()
            # print(out)

        crysol_List=[]
        for i in fileList:
            crysol_List.append(i + '.log')

            
        crysol_fits = []
        for k in fileList:
            crysol_fits.append(k+'.fit')

        def load_fitFiles(fileList):
          '''
          Describe function here:


          '''
          c=1
          diction={}
          for i in fileList:
            diction['data_%s'%str(c)]=np.loadtxt(i, dtype={'names': ('Q', 'I(Q)','Q_mod','I(Q)_mod'), 'formats': (np.float,np.float,np.float,np.float)}, comments='#',
              skiprows=1)
            c+=1
    
          return diction

        crysol_files = load_fitFiles(fileList=crysol_fits)

        export_dictionary={}
            
        for fn,k,z in zip(crysol_List,crysol_files,fileList):
            crysol_parse=[]
            crysol_parse2=[]
            with open(fn, "r") as f: # important to use with command when dealing with files
                counter = 0
                print('File: %s' %str(z))
                c=0
                for line in f:
                    c+=1
                    if 'Chi^2' in line:
                        crysol_parse.append(''.join(islice(f, 1)))
                    if 'Rg from the slope of net intensity' in line:
                        # print(line)
                        crysol_parse2.append(line)
            
            Rg=[i for j in crysol_parse2[0].split() for i in (j, ' ')][:-1][18]
            print('Final Rg (A): ', Rg)
            chi_sq=[i for j in crysol_parse[0].split() for i in (j, ' ')][:-1][18]
            print('Chi-Squared: ',chi_sq)
            print('')

            Q = crysol_files[str(k)]['Q']
            expt_IQ = crysol_files[str(k)]['I(Q)']
            model_IQ = crysol_files[str(k)]['I(Q)_mod']

            export_dictionary['%s'%str(z)] = {'Q': Q,'exptIQ': expt_IQ, 'modelIQ': model_IQ}

            # pairList=[[crysol_files[str(k)]['Q'],crysol_files[str(k)]['I(Q)']],
            #           [crysol_files[str(k)]['Q'],crysol_files[str(k)]['I(Q)_mod']]]
            # labelList=['%s'%str(z),'Model']
            # colorList=['#9B9B9B','#00B2AF']
            # X1=crysol_files[str(k)]['Q']
            # Y1=crysol_files[str(k)]['I(Q)']
            # X2=crysol_files[str(k)]['Q']
            # Y2=crysol_files[str(k)]['I(Q)_mod']
            # Y1err=[1] * len(X1)
          
            # skip=next((i for i, x in enumerate(Y1) if x), None) # x!= 0 for strict match

            # figs.vertical_stackPlot(X1=X1[skip:],Y1=Y1[skip:],Y1err=Y1err[skip:],X2=X2[skip:],Y2=Y2[skip:],ylabel1='Intensity, I(Q)',
            #                         ylabel2='Residuals',xlabel='Q ($\mathring{A}^{-1}$)',
            #                         Label1='%s'%z, Label2='Crysol Model',saveLabel='%s_Crysol_fit_EM10mer'%str(z),
            #                         bottomPlot_yLabel='$ln(\\frac{I_{expt}(q)}{I_{model}(q)}) \cdot (\\frac{1}{\sigma_{expt}})$',
            #                         LinLin=False,
            #                         linewidth=4,labelSize=20,darkmode=False,plot=True)

        return export_dictionary


    def runCrysol_single(self,
                         filename,
                         PDB_File,
                         qmax=0.3):
        '''
        Description
        '''

        cwd = str(os.getcwd())
        crysol_dir = '/Crysol_Modeling' # create folder for modeling files
        save_dir = cwd + crysol_dir

        if not os.path.exists(save_dir):
          os.mkdir(save_dir) # create it if it doesn't already exist

        os.chdir(save_dir) # switch to that directory

        print(save_dir)

        # Import data files

        fileList = filename

        test = os.path.basename(fileList)
        print(test)

        output_name = cwd + crysol_dir + '/fileList'

        ## Crysol Modeling

        print('CRYSOL RUNNING')
        # n=len(datFile['ERROR']) # scale the size of the profile for downstream residuals calculations
        cmd = 'crysol %s %s -lm 50 -sm %s -kp=ii -p %s'%(str(PDB_File),str(filename),str(qmax),str(output_name))
        print(cmd)
        proc = subprocess.Popen([cmd], stdout=subprocess.PIPE, shell=True)
        (out, err) = proc.communicate()
        # print(out)

        crysol_List=[]
        crysol_List.append(output_name + '.log')

        print(crysol_List)

            
        crysol_fits = []
        crysol_fits.append(output_name + '.fit')

        print(crysol_fits)

        ##### Fancier parsing..


        crysol_files = self.load_fitFiles(fileList=crysol_fits)

        export_dictionary={}
            
        for fn,k,z in zip(crysol_List,crysol_files,fileList):
            crysol_parse=[]
            crysol_parse2=[]
            with open(fn, "r") as f: # important to use with command when dealing with files
                counter = 0
                print('File: %s' %str(z))
                c=0
                for line in f:
                    c+=1
                    if 'Chi^2' in line:
                        crysol_parse.append(''.join(islice(f, 1)))
                    if 'Rg from the slope of net intensity' in line:
                        # print(line)
                        crysol_parse2.append(line)
            
            Rg=[i for j in crysol_parse2[0].split() for i in (j, ' ')][:-1][18]
            print('Final Rg (A): ', Rg)
            chi_sq=[i for j in crysol_parse[0].split() for i in (j, ' ')][:-1][18]
            print('Chi-Squared: ',chi_sq)
            print('')

            Q = crysol_files[str(k)]['Q']
            expt_IQ = crysol_files[str(k)]['I(Q)']
            model_IQ = crysol_files[str(k)]['I(Q)_mod']

            export_dictionary['%s'%str(z)] = {'Q': Q,'exptIQ': expt_IQ, 'modelIQ': model_IQ}

        return export_dictionary


        def load_fitFiles(file):
          '''
          Just loads in the intensity file generated by crysol.
          '''
          c=1
          diction={}
          diction['data_%s'%str(c)]=np.loadtxt(file, dtype={'names': ('Q', 'I(Q)','Q_mod','I(Q)_mod'), 
                                               'formats': (np.float,np.float,np.float,np.float)}, comments='#',
                                               skiprows=1)
          return diction

        # crysol_files = load_fitFiles(file=crysol_fits)
            
        # crysol_parse=[]
        # crysol_parse2=[]
        # with open(crysol_List, "r") as f: # important to use with command when dealing with files
        #     counter = 0
        #     print('File: %s' %str(filename))
        #     c=0
        #     for line in f:
        #         c+=1
        #         if 'Chi^2' in line:
        #             crysol_parse.append(''.join(islice(f, 1)))
        #         if 'Rg from the slope of net intensity' in line:
        #             print(line)
        #             crysol_parse2.append(line)
        
        # Rg=[i for j in crysol_parse2[0].split() for i in (j, ' ')][:-1][18]
        # print('Final Rg (A): ', Rg)
        # chi_sq=[i for j in crysol_parse[0].split() for i in (j, ' ')][:-1][18]
        # print('Chi-Squared: ',chi_sq)

        # Q = crysol_files[str(crysol_files)]['Q']
        # expt_IQ = crysol_files[str(crysol_files)]['I(Q)']
        # model_IQ = crysol_files[str(crysol_files)]['I(Q)_mod']

        # pairList=[[crysol_files[str(crysol_files)]['Q'],crysol_files[str(crysol_files)]['I(Q)']],
        #           [crysol_files[str(crysol_files)]['Q'],crysol_files[str(crysol_files)]['I(Q)_mod']]]
        # labelList=['%s'%str(filename),'Model']
        # colorList=['#9B9B9B','#00B2AF']
        # X1=crysol_files[str(crysol_files)]['Q']
        # Y1=crysol_files[str(crysol_files)]['I(Q)']
        # X2=crysol_files[str(crysol_files)]['Q']
        # Y2=crysol_files[str(crysol_files)]['I(Q)_mod']
        # Y1err=[1] * len(X1)
      
        # skip=next((i for i, x in enumerate(Y1) if x), None) # x!= 0 for strict match

        # figs.vertical_stackPlot(X1=X1[skip:],Y1=Y1[skip:],Y1err=Y1err[skip:],X2=X2[skip:],Y2=Y2[skip:],ylabel1='Intensity, I(Q)',
        #                         ylabel2='Residuals',xlabel='Q ($\mathring{A}^{-1}$)',
        #                         Label1='%s'%filename, Label2='Crysol Model',saveLabel='%s_Crysol_fit_EM10mer'%str(filename),
        #                         bottomPlot_yLabel='$ln(\\frac{I_{expt}(q)}{I_{model}(q)}) \cdot (\\frac{1}{\sigma_{expt}})$',
        #                         LinLin=False,
        #                         linewidth=4,labelSize=20,darkmode=False,plot=True)

        return Q, expt_IQ, model_IQ

class IFTM(object):
    """
    Inverse fourier tranform measurement (IFTM) object. Contains the P(r), r
    and error vectors, as well as the original data, the fit of the P(r) to
    the data, and all associated metadata.

    Attributes
    ----------
    r: numpy.array
        The r vector of the P(r) function.
    p: numpy.array
        The values of the P(r) function.
    err: numpy.array
        The errors of the P(r) function.
    q_orig: numpy.array
        The q vector of the input data.
    i_orig: numpy.array
        The intensity vector of the input data.
    err_orig: numpy.array
        The error vector of the input data.
    i_fit: numpy.array
        The intensity vector of the fit to the input data.
    q_extrap: numpy.array
        The q vector of the input data extrapolated to q=0.
    i_extrap: numpy.array
        The intensity vector of the fit to the input data extrapolated to q=0.
    """

    def __init__(self, p, r, err, i_orig, q_orig, err_orig, i_fit, parameters, i_extrap = [], q_extrap = []):
        """
        Constructor

        Parameters
        ----------
        p: numpy.array
            The input P(r) values.
        r:  numpy.array
            The input r values for the P(r) function.
        err: numpy.array
            The input error values for the P(r) function.
        i_orig: numpy.array
            The intensity values of the data used to do the IFT.
        q_orig: numpy.array
            The q values of the data used to do the IFT.
        err_orig: numpy.array
            The error values of the data used to do the IFT.
        i_fit: numpy.array
            The intensity values of the fit of the P(r) function to the data.
        parameters: dict
            A dictionary of the metadata. Should ontain at least {'filename':
            filename_with_no_path}
        i_extrap: numpy.array, optional
            The intensity values of the fit of the P(r) function to the data
            extrapolated to q=0. If not provided, an empty array is used.
        q_extrap: numpy.array, optional
            The q values of the input data extrapolated to q=0. If not
            provided, an empty array is used.
        """

        #Raw intensity variables
        self._r_raw = np.array(r)
        self._p_raw = np.array(p)
        self._err_raw = np.array(err)
        self._i_orig_raw = np.array(i_orig)
        self._q_orig_raw = np.array(q_orig)
        self._err_orig_raw = np.array(err_orig)
        self._i_fit_raw = np.array(i_fit)
        self._i_extrap_raw = np.array(i_extrap)
        self._q_extrap_raw = np.array(q_extrap)
        self._parameters = parameters

        # Make an entry for analysis parameters i.e. Rg, I(0) etc:
        # if 'analysis' not in self._parameters:
        #     self._parameters['analysis'] = {}
        # if 'history' not in self._parameters:
        #     self._parameters['history'] = {}

        #Modified intensity variables
        self.r = self._r_raw.copy()
        self.p = self._p_raw.copy()
        self.err = self._err_raw.copy()

        self.i_orig = self._i_orig_raw.copy()
        self.q_orig = self._q_orig_raw.copy()
        self.err_orig = self._err_orig_raw.copy()

        self.i_fit = self._i_fit_raw.copy()

        self.i_extrap = self._i_extrap_raw.copy()
        self.q_extrap = self._q_extrap_raw.copy()

        #variables used for plot management
        self.item_panel = None

        self.plot_panel = None

        self.r_line = None
        self.qo_line = None
        self.qf_line = None

        self.r_origline = None
        self.qo_origline = None
        self.qf_origline = None

        self.r_err_line = None
        self.qo_err_line = None

        self.r_axes = None
        self.qo_axes = None
        self.qf_axes = None

        self.canvas = None

        self.is_plotted = False

    def __deepcopy__(self, memo):
        p = copy.deepcopy(self._p_raw, memo)
        r = copy.deepcopy(self._r_raw, memo)
        err = copy.deepcopy(self._err_raw, memo)

        i_orig = copy.deepcopy(self._i_orig_raw, memo)
        q_orig = copy.deepcopy(self._q_orig_raw, memo)
        err_orig = copy.deepcopy(self._err_orig_raw, memo)

        i_fit = copy.deepcopy(self._i_fit_raw, memo)

        i_extrap = copy.deepcopy(self._i_extrap_raw, memo)
        q_extrap = copy.deepcopy(self._q_extrap_raw, memo)

        parameters = copy.deepcopy(self._parameters, memo)

        new_iftm = IFTM(p, r, err, i_orig, q_orig, err_orig, i_fit, parameters,
        i_extrap, q_extrap)

        return new_iftm

    def getScale(self):
        """
        Returns the scale factor for the P(r) function.

        Returns
        -------
        scale: float
            The scale factor.
        """
        return self._scale_factor

    def getOffset(self):
        """
        Returns the offset for the P(r) function.

        Returns
        -------
        offset: float
            The offset.
        """
        return self._offset_value

    def getLine(self):
        """
        Returns the plotted line for the P(r) function. Only used in the RAW GUI.

        Returns
        -------
        line: matplotlib.lines.Line2D
            The plotted line.
        """
        return self.line

    def setAllParameters(self, new_parameters):
        """
        Sets the parameters dictionary, which contains the IFT metadata,
        to the new input value.

        Parameters
        ----------
        new_parameters: dict
            A dictionary containing the new parameters.
        """
        self._parameters = new_parameters

    def getAllParameters(self):
        """
        Returns all of the metadata parameters associated with the IFT as
        a dictionary.

        Returns
        -------
        parameters: dict
            The metadata associated with the IFT.
        """
        return self._parameters

    def getParameter(self, key):
        """
        Gets a particular metadata parameter based on the provided key.

        Parameters
        ----------
        key: str
            A string that is a key in the parameters metadata dictionary.

        Returns
        -------
        parameter
            The parameter associated with the specified key. If the key is not
            in the parameter dictionary, None is returned.
        """

        if key in self._parameters:
            return self._parameters[key]
        else:
            return None

    def setParameter(self, key, value):
        """
        Sets a particular metadata parameter based on the provided key and value.

        Parameters
        ----------
        key: str
            The name of the new bit of metadata.
        value: object
            The value of the new bit of metadata. Could be anything that is
            an acceptable value for a dictionary.
        """
        self._parameters[key] = value

    def extractAll(self):
        """
        Extracts the raw and scaled q, intensity, and error, the scale and
        offset values, the selected q range, and the parameters in a single
        dictionary.

        Returns
        -------
        all_data: dict
            A dictionary with keys r_raw, p_raw, err_raw, i_orig_raw, q_orig_raw,
            err_orig_raw, i_fit_raw, i_extrap_raw, q_extrap_raw, and parameters,
            which correspond to those values from the IFTM.
        """

        all_data = {}

        all_data['r_raw'] = self._r_raw
        all_data['p_raw'] = self._p_raw
        all_data['err_raw'] = self._err_raw

        all_data['i_orig_raw'] = self._i_orig_raw
        all_data['q_orig_raw'] = self._q_orig_raw
        all_data['err_orig_raw'] = self._err_orig_raw

        all_data['i_fit_raw'] = self._i_fit_raw
        all_data['i_extrap_raw'] = self._i_extrap_raw
        all_data['q_extrap_raw'] = self._q_extrap_raw

        all_data['parameters'] = self._parameters

        return all_data

        pass

    def copy(self):
        """
        Creates a copy of the IFT.

        Returns
        -------
        ift: bioxtasraw.SASM.IFTM
            The copied IFTM
        """

        iftm_copy = IFTM(copy.deepcopy(self._p_raw), copy.deepcopy(self._r_raw),
            copy.deepcopy(self._err_raw), copy.deepcopy(self._i_orig_raw),
            copy.deepcopy(self._q_orig_raw), copy.deepcopy(self._err_orig_raw),
            copy.deepcopy(self._i_fit_raw), copy.deepcopy(self._parameters),
            copy.deepcopy(self._i_extrap_raw), copy.deepcopy(self._q_extrap_raw))

        return iftm_copy


def postProcessSasm(sasm, raw_settings):

    if raw_settings.get('ZingerRemoval'):
        std = raw_settings.get('ZingerRemoveSTD')
        winlen = raw_settings.get('ZingerRemoveWinLen')
        start_idx = raw_settings.get('ZingerRemoveIdx')

        sasm.removeZingers(start_idx, winlen, std



                           )


class SAXS_FileReader(object):
    '''
    Must pass in an object with all of the __init__ stuff ready to go!
    '''

    def __init__(self, i, q, err, parameters, settings = None):
        """
        Constructor

        Parameters
        ----------
        i: numpy.array
            The intensity vector.
        q: numpy.array
            The q vector.
        err: numpy.array
            The error vector.
        parameters: dict
            A dictionary of metadata for the object. This should contain at
            least {'filename': filename_with_no_path}. Other reserved keys are:
            'counters' : [(countername, value),...] Info from counterfiles
            'fileHeader' : [(label, value),...] Info from the header in the
            loaded file
        """

        #Raw intensity variables
        self._i_raw = np.array(i)
        self._q_raw = np.array(q)
        self._err_raw = np.array(err)
        self._parameters = parameters

        # Make an entry for analysis parameters i.e. Rg, I(0) etc:
        if 'analysis' not in self._parameters:
            self._parameters['analysis'] = {}
        if 'history' not in self._parameters:
            self._parameters['history'] = {}

        #Modified intensity variables
        self.i = self._i_raw.copy()
        self.q = self._q_raw.copy()
        self.err = self._err_raw.copy()

        self._scale_factor = 1
        self._offset_value = 0
        self._q_scale_factor = 1

        #variables used for plot management
        self.item_panel = None
        self.plot_panel = None
        self.line = None
        self.err_line = None
        self.axes = None
        self.is_plotted = False
        self._selected_q_range = (0, len(self._q_raw))

        #Calculated values
        try:
            if len(self.q)>0:
                self.total_intensity = integrate.trapz(self.getI(), self.getQ())
                self.mean_intensity = self.getI().mean()
            else:
                self.total_intensity = -1
                self.mean_intensity = -1

        except Exception as e:
            print(e)
            self.total_intensity = -1
            self.mean_intensity = -1
    # def SafeUnpickler(pickle.Unpickler):
    #     find_class = staticmethod(find_global)



    def readSettings(filename):
        '''
        Describe here...
        '''
        try:
            with open(filename, 'r') as f:
                settings = f.read()
            settings = dict(json.loads(settings))
        except Exception as e:
            try:
                with open(filename, 'rb') as f:
                    if six.PY3:
                        pickle_obj = self.SafeUnpickler(f, encoding='latin-1')
                    else:
                        pickle_obj = pickle.Unpickler(f)
                        pickle_obj.find_global = SASUtils.find_global
                    settings = pickle_obj.load()
            except Exception as e:
                print('Error type: %s, error: %s' %(type(e).__name__, e))
                return None

        return settings

    def checkFileType(filename):
        ''' Tries to find out what file type it is and reports it back '''

        path, ext = os.path.splitext(filename)

        if ext == '.fit':
            return 'fit'
        elif ext == '.fir':
            return 'fir'
        elif ext == '.out':
            return 'out'
        elif ext == '.nxs': #Nexus file
            return 'image'
        elif ext == '.edf':
            return 'image'
        elif ext == '.ccdraw':
            return 'image'
        elif ext == '.int':
            return 'int'
        elif ext == '.img' or ext == '.imx_0' or ext == '.dkx_0' or ext == '.dkx_1' or ext == '.png' or ext == '.mpa':
            return 'image'
        elif ext == '.dat' or ext == '.sub' or ext =='.txt':
            return 'primus'
        elif ext == '.mar1200' or ext == '.mar2400' or ext == '.mar2300' or ext == '.mar3600':
            return 'image'
        elif (ext == '.img' or ext == '.sfrm' or ext == '.dm3' or ext == '.edf' or ext == '.xml' or ext == '.cbf' or ext == '.kccd' or
            ext == '.msk' or ext == '.spr' or ext == '.tif' or ext == '.mccd' or ext == '.mar3450' or ext =='.npy' or
            ext == '.pnm' or ext == '.No'):
            return 'image'
        elif ext == '.ift':
            return 'ift'
        elif ext == '.csv':
            return 'csv'
        elif ext == '.h5':
            return 'hdf5'
        else:
            try:
                fabio.open(filename)
                return 'image'
            except Exception:
                try:
                    float(ext.strip('.'))
                except Exception:
                    return 'rad'
                return 'csv'

    def loadFabio(filename, hdf5_file=None):
        if hdf5_file is None:
            fabio_img = fabio.open(filename)
        else:
            fabio_img = hdf5_file

        if fabio_img.nframes == 1:
            data = fabio_img.data
            hdr = fabio_img.getheader()

            img = [data]
            img_hdr = [hdr]

        else:
            img = [None for i in range(fabio_img.nframes)]
            img_hdr = [None for i in range(fabio_img.nframes)]

            img[0] = fabio_img.data
            img_hdr[0] = fabio_img.getheader()

            for i in range(1,fabio_img.nframes):
                fabio_img = fabio_img.next()
                img[i] = fabio_img.data
                img_hdr[i] = fabio_img.getheader()

        fabio_img.close()

        return img, img_hdr

    def parseCHESSEIGER4MCountFile(filename):
        ''' Loads information from the counter file at CHESS, id7a from
        the image filename. EIGER .h5 files with 1-based frame numbers '''
        dir, file = os.path.split(filename)
        underscores = file.split('_')

        countFile = underscores[0]

        filenumber = int(underscores[-3])

        try:
            frame_number = int(underscores[-1].split('.')[0])
        except Exception:
            frame_number = 0

        # REG: if user root file name contains underscores, include those
        # note: must start at -3 to leave out "data" in image name

        if len(underscores)>3:
            for each in underscores[1:-3]:
                countFile += '_' + each

        countFilename = os.path.join(dir, countFile)

        with open(countFilename,'rU') as f:
            allLines = f.readlines()

        line_num = 0
        start_found = False
        start_idx = None
        label_idx = None
        date_idx = None

        for eachLine in allLines:
            splitline = eachLine.split()

            if len(splitline) > 1:
                if splitline[0] == '#S' and splitline[1] == str(filenumber):
                    start_found = True
                    start_idx = line_num

                if splitline[0] == '#D' and start_found:
                    date_idx = line_num

                if splitline[0] == '#L' and start_found:
                    label_idx = line_num
                    break

            line_num = line_num + 1

        counters = {}
        try:
            if start_idx and label_idx:
                labels = allLines[label_idx].split()
                # REG: hdf5 indices start at 1 not 0 as was our Pilatus convention!
                vals = allLines[label_idx+0+frame_number].split()

            for idx in range(0,len(vals)):
                counters[labels[idx+1]] = vals[idx]

            if date_idx:
                counters['date'] = allLines[date_idx][3:-1]

        except:
            print('Error loading CHESS id7a counter file')

        return counters

    def loadImage(filename, raw_settings, hdf5_file=None):
        ''' returns the loaded image based on the image filename
        and image type. '''
        image_type = raw_settings.get('ImageFormat')
        fliplr = raw_settings.get('DetectorFlipLR')
        flipud = raw_settings.get('DetectorFlipUD')

        try:
            if all_image_types[image_type] == loadFabio:
                img, imghdr = all_image_types[image_type](filename, hdf5_file)
            else:
                img, imghdr = all_image_types[image_type](filename)
        except (ValueError, TypeError, KeyError, fabio.fabioutils.NotGoodReader, Exception) as msg:
            raise SASExceptions.WrongImageFormat('Error loading image, ' + str(msg))

        if not isinstance(img, list):
            img = [img]
        if not isinstance(imghdr, list):
            imghdr = [imghdr]

        #Clean up headers by removing spaces in header names and non-unicode characters)
        for hdr in imghdr:
            if hdr is not None:
                hdr = {key.replace(' ', '_').translate(str.maketrans('', '', '()[]'))
                    if isinstance(key, str) else key: hdr[key] for key in hdr}
                # hdr = { key : str(hdr[key], errors='ignore') if isinstance(hdr[key], str)
                #     else hdr[key] for key in hdr}


        if image_type != 'SAXSLab300':
            for i in range(len(img)):
                if fliplr:
                    img[i] = np.fliplr(img[i])
                if flipud:
                    img[i] = np.flipud(img[i])
        return img, imghdr

#################################
#--- ** MAIN LOADING FUNCTION **
#################################

    def loadFile(filename, raw_settings, no_processing = False):
        ''' Loads a file an returns a SAS Measurement Object (SASM) and the full image if the
            selected file was an Image file

             NB: This is the function used to load any type of file in RAW
        '''
        try:
            file_type = self.checkFileType(filename)
            # print file_type
        except IOError:
            raise
        except Exception as msg:
            print(str(msg))
            file_type = None

        if file_type == 'hdf5':
            try:
                hdf5_file = fabio.open(filename)
                file_type = 'image'
            except Exception:
                pass
        else:
            hdf5_file = None

        if file_type == 'image':
            try:
                sasm, img = self.loadImageFile(filename, raw_settings, hdf5_file)
            except (ValueError, AttributeError) as msg:
                # raise SASExceptions.UnrecognizedDataFormat('No data could be retrieved from the file, unknown format.')
                print('No data could be retrieved.. unknown format')

            #Always do some post processing for image files
            if not isinstance(sasm, list):
                sasm = [sasm]

            for current_sasm in sasm:
                self.postProcessProfile(current_sasm, raw_settings, no_processing)

        elif file_type == 'hdf5':
            sasm = loadHdf5File(filename, raw_settings)
            img = None
        else:
            print('FAILURE!!!')
            # sasm = loadAsciiFile(filename, file_type)
            # img = None

            #If you don't want to post process asci files, return them as a list
            if not isinstance(sasm, list):
                SASM.postProcessSasm(sasm, raw_settings)

        if not isinstance(sasm, list) and (sasm is None or len(sasm.i) == 0):
            # raise SASExceptions.UnrecognizedDataFormat('No data could be retrieved from the file, unknown format.')
            print('unknown format..')

        return sasm, img

    def postProcessSasm(sasm, raw_settings):
        '''
        Description..
        '''
        if raw_settings.get('ZingerRemoval'):
            std = raw_settings.get('ZingerRemoveSTD')
            winlen = raw_settings.get('ZingerRemoveWinLen')
            start_idx = raw_settings.get('ZingerRemoveIdx')

            sasm.removeZingers(start_idx, winlen, std)


    def postProcessProfile(sasm, raw_settings, no_processing):
        """
        Does post-processing on profiles created from images.
        """
        self.postProcessSasm(sasm, raw_settings) #Does dezingering

        if not no_processing:
            #Need to do a little work before we can do glassy carbon normalization
            if raw_settings.get('NormAbsCarbon') and not raw_settings.get('NormAbsCarbonIgnoreBkg'):
                bkg_filename = raw_settings.get('NormAbsCarbonSamEmptyFile')
                bkg_sasm = raw_settings.get('NormAbsCarbonSamEmptySASM')
                if bkg_sasm is None or bkg_sasm.getParameter('filename') != os.path.split(bkg_filename)[1]:
                    bkg_sasm, junk_img = self.loadFile(bkg_filename, raw_settings, no_processing=True)
                    if isinstance(bkg_sasm,list):
                        if len(bkg_sasm) > 1:
                            bkg_sasm = SASProc.average(bkg_sasm) # RM!
                        else:
                            bkg_sasm = bkg_sasm[0]
                    raw_settings.set('NormAbsCarbonSamEmptySASM', bkg_sasm)

            try:
                #Does fully glassy carbon abs scale
                SASCalib.postProcessImageSasm(sasm, raw_settings)
            except SASExceptions.AbsScaleNormFailed:
                raise

    def postProcessImageSasm(sasm, raw_settings):
        if (raw_settings.get('NormAbsCarbon') and
            not raw_settings.get('NormAbsCarbonIgnoreBkg')):
            self.normalizeAbsoluteScaleCarbon(sasm, raw_settings)

    def normalizeAbsoluteScaleCarbon(sasm, raw_settings):
        abs_scale_constant = float(raw_settings.get('NormAbsCarbonConst'))
        sam_thick = float(raw_settings.get('NormAbsCarbonSamThick'))

        bkg_sasm = raw_settings.get('NormAbsCarbonSamEmptySASM')

        ctr_ups = raw_settings.get('NormAbsCarbonUpstreamCtr')
        ctr_dns = raw_settings.get('NormAbsCarbonDownstreamCtr')

        sample_ctrs = sasm.getParameter('imageHeader')
        sample_file_hdr = sasm.getParameter('counters')
        sample_ctrs.update(sample_file_hdr)

        bkg_ctrs = bkg_sasm.getParameter('imageHeader')
        bkg_file_hdr = bkg_sasm.getParameter('counters')
        bkg_ctrs.update(bkg_file_hdr)

        sample_ctr_ups_val = float(sample_ctrs[ctr_ups])
        sample_ctr_dns_val = float(sample_ctrs[ctr_dns])
        bkg_ctr_ups_val = float(bkg_ctrs[ctr_ups])
        bkg_ctr_dns_val = float(bkg_ctrs[ctr_dns])

        sample_trans = (sample_ctr_dns_val/sample_ctr_ups_val)/(bkg_ctr_dns_val/bkg_ctr_ups_val)
        sasm.scaleRawIntensity(1./sample_ctr_ups_val)
        bkg_sasm.scale((1./bkg_ctr_ups_val)*sample_trans)

        try:
            sub_sasm = SASProc.subtract(sasm, bkg_sasm, forced = True, full = True)
        except SASExceptions.DataNotCompatible:
            sasm.scaleRawIntensity(sample_ctr_ups_val)
            raise SASExceptions.AbsScaleNormFailed('Absolute scale failed because empty scattering could not be subtracted')

        sub_sasm.scaleRawIntensity(1./(sample_trans)/sam_thick)
        sub_sasm.scaleRawIntensity(abs_scale_constant)

        sasm.setRawQ(sub_sasm.getRawQ())
        sasm.setRawI(sub_sasm.getRawI())
        sasm.setRawErr(sub_sasm.getRawErr())
        sasm.scale(1.)
        sasm.setQrange((0,len(sasm.q)))

        bkg_sasm.scale(1.)

        abs_scale_params = {
            'Method'    : 'Glassy_carbon',
            'Absolute_scale_factor': abs_scale_constant,
            'Sample_thickness_[mm]': sam_thick,
            'Ignore_background': False,
            }

        abs_scale_params['Background_file'] = raw_settings.get('NormAbsCarbonSamEmptyFile')
        abs_scale_params['Upstream_counter_name'] = ctr_ups
        abs_scale_params['Downstream_counter_name'] = ctr_dns
        abs_scale_params['Upstream_counter_value_sample'] = sample_ctr_ups_val
        abs_scale_params['Downstream_counter_value_sample'] = sample_ctr_dns_val
        abs_scale_params['Upstream_counter_value_background'] = bkg_ctr_ups_val
        abs_scale_params['Downstream_counter_value_background'] = bkg_ctr_dns_val
        abs_scale_params['Sample_transmission'] = sample_trans

        norm_parameter = self.getParameter('normalizations')

        norm_parameter['Absolute_scale'] = abs_scale_params

        self.setParameter('normalizations', norm_parameter)

        return sasm, abs_scale_constant

    def getParameter(self, key):
        """
        Gets a particular metadata parameter based on the provided key.

        Parameters
        ----------
        key: str
            A string that is a key in the parameters metadata dictionary.

        Returns
        -------
        parameter
            The parameter associated with the specified key. If the key is not
            in the parameter dictionary, None is returned.
        """

        if key in self._parameters:
            return self._parameters[key]
        else:
            return None
            
    def setParameter(self, key, value):
        """
        Sets a particular metadata parameter based on the provided key and value.

        Parameters
        ----------
        key: str
            The name of the new bit of metadata.
        value: object
            The value of the new bit of metadata. Could be anything that is
            an acceptable value for a dictionary.
        """
        self._parameters[key] = value


    def removeZingers(self, start_idx = 0, window_length = 10, stds = 4.0):
        """
        Removes spikes (zingers) from radially averaged data based on smoothing
        of the intensity profile.

        Parameters
        ----------
        start_idx: int
            The initial index in the intensity to start the dezingering process.
        window_length: int
            The size of the window used to search for an replace spikes, as the
            number of intensity points.
        stds: The standard deviation threshold used to detect spikes.
        """

        intensity = self._i_raw

        for i in range(window_length + start_idx, len(intensity)):

            averaging_window = intensity[i - window_length : i]
            averaging_window_std = np.std(averaging_window)
            averging_window_mean = np.mean(averaging_window)

            threshold = averging_window_mean + (stds * averaging_window_std)

            if intensity[i] > threshold:
                intensity[i] = averging_window_mean

        print('Ran removeZingers')

    def loadImageFile(filename, raw_settings, hdf5_file=None):
        hdr_fmt = raw_settings.get('ImageHdrFormat')

        loaded_data, loaded_hdr = self.loadImage(filename, raw_settings, hdf5_file)

        sasm_list = [None for i in range(len(loaded_data))]

        #Process all loaded images into sasms
        for i in range(len(loaded_data)):
            img = loaded_data[i]
            img_hdr = loaded_hdr[i]

            if len(loaded_data) > 1:
                temp_filename = os.path.split(filename)[1].split('.')
                if len(temp_filename) > 1:
                    temp_filename[-2] = temp_filename[-2] + '_%05i' %(i+1)
                else:
                    temp_filename[0] = temp_filename[0] + '_%05i' %(i+1)

                new_filename = '.'.join(temp_filename)
            else:
                new_filename = os.path.split(filename)[1]

            hdrfile_info = self.loadHeader(filename, new_filename, hdr_fmt)

            parameters = {'imageHeader' : img_hdr,
                          'counters'    : hdrfile_info,
                          'filename'    : new_filename,
                          'load_path'   : filename}

            sasm = self.processImage(img, parameters, raw_settings)

            sasm_list[i] = sasm


        return sasm_list, loaded_data












