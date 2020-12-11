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