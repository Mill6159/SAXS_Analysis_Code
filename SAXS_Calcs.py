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
from numba import jitclass, types, typed
import math
import warnings
import traceback
# from PlotClass import *
# from FileParser import *

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

        sasm.removeZingers(start_idx, winlen, std)
