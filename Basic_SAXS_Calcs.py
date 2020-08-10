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
from PlotClass import *
from SAXS_Calcs import *
# import PlotClass as pc
from FileParser import *

#######################################################
# Define Class
#######################################################


# @jitclass([('d', types.DictType(types.float64)),
#            ('l', types.ListType(types.float64))])
class BasicSAXS:

    def __init__(self,notify=True,atsas_dir=' '):
        '''
        Describe the class here.
        Ideally, this function could read and determine whether or not it is a .dat file or OTHER..
        and deal with an input file accordingly... hmmm. Perhaps, I need another class (i.e. data reader?)
        '''
        self.notify = notify
        if self.notify == True:
            print('--------------------------------------------------------------')
            print('Basic SAXS Calculations Class was called')


        self.atsas_dir=atsas_dir
        if atsas_dir==' ':
            print('##############################################################')
            print('No GNOM directory provided, attempting to find it...')
            get_dammif_dir=subprocess.Popen('cd; which dammif', shell=True, stdout=subprocess.PIPE).stdout
            dammif_dir=get_dammif_dir.read()
            dammif_dir_s2=dammif_dir.decode()
            dammif_dir_s3=re.sub(r'\w*dammif\w*', '', dammif_dir_s2)
            gnom_dir=dammif_dir_s3.strip()+'gnom'
            self.atsas_dir=gnom_dir # setting gnom directory finally
            if '/bin/' not in self.atsas_dir:
                print('We were NOT able to locate the ATSAS-GNOM library on your local computer. Define this manually when calling the BasicSAXS() class.')
                print('We will terminate the script...')
                sys.exit('Analysis terminated!!!')
            else:
                print('The GNOM library was found in the directory: %s. \nIf this is the incorrect library, the calculation may fail'%gnom_dir)
            print('##############################################################')
        else:
            print('#######################')
            print('GNOM directory provided')
            print('#######################')
            print('--------------------------------------------------------------')

        self.plots=PlotClass(notify=False)
        self.FileParser=FileParser(notify=False)
        self.calcs=SAXSCalcs(notify=False)
        # self.d=typed.List.empty_list(types.float64)
        # self.l=typed.List.empty_list(types.float64)

    def lineModel(self,x,m,b):
        return m*x+b

    '''
    Input some functions necessary for AUTORG calculations
    These require numba, which speeds things up significantly.
    '''

    # @jit(nopython=True,cache=False,parallel=False)
    def linear_func(self,x,a,b):
        return a + b * x

    # @jit(nopython=True,cache=True,parallel=False)
    def lin_reg(self,x,y):
        x_sum = x.sum()
        xsq_sum = (x ** 2).sum()
        y_sum = y.sum()
        xy_sum = (x * y).sum()
        n = len(x)

        delta = n * xsq_sum - x_sum ** 2.

        if delta != 0:
            a = (xsq_sum * y_sum - x_sum * xy_sum) / delta
            b = (n * xy_sum - x_sum * y_sum) / delta

            cov_y = (1. / (n - 2.)) * ((y - a - b * x) ** 2.).sum()
            cov_a = cov_y * (xsq_sum / delta)
            cov_b = cov_y * (n / delta)
        else:
            a = -1
            b = -1
            cov_a = -1
            cov_b = -1

        return a,b,cov_a,cov_b

    # @jit(nopython=True,cache=True,parallel=False)
    def weighted_lin_reg(self,x,y,err):
        weights = 1. / (err) ** 2.

        w_sum = weights.sum()
        wy_sum = (weights * y).sum()
        wx_sum = (weights * x).sum()
        wxsq_sum = (weights * x ** 2.).sum()
        wxy_sum = (weights * x * y).sum()

        delta = weights.sum() * wxsq_sum - (wx_sum) ** 2.

        if delta != 0:
            a = (wxsq_sum * wy_sum - wx_sum * wxy_sum) / delta
            b = (w_sum * wxy_sum - wx_sum * wy_sum) / delta

            cov_a = wxsq_sum / delta
            cov_b = w_sum / delta
        else:
            a = -1
            b = -1
            cov_a = -1
            cov_b = -1

        return a,b,cov_a,cov_b

    # @jit(nopython=True,cache=False,parallel=False)
    def calcRg(self,q,i,err,transform=True,error_weight=True):
        if transform:
            # Start out by transforming as usual.
            x = np.square(q)
            y = np.log(i)
            yerr = np.absolute(err / i)  # I know it looks odd, but it's correct for a natural log
        else:
            x = q
            y = i
            yerr = err

        if error_weight:
            a,b,cov_a,cov_b = self.weighted_lin_reg(x,y,yerr)
        else:
            a,b,cov_a,cov_b = self.lin_reg(x,y)

        if b < 0:
            RG = np.sqrt(-3. * b)
            I0 = np.exp(a)

            # error in rg and i0 is calculated by noting that q(x)+/-Dq has Dq=abs(dq/dx)Dx, where q(x) is your function you're using
            # on the quantity x+/-Dx, with Dq and Dx as the uncertainties and dq/dx the derviative of q with respect to x.
            RGer = np.absolute(0.5 * (np.sqrt(-3. / b))) * np.sqrt(np.absolute(cov_b))
            I0er = I0 * np.sqrt(np.absolute(cov_a))

        else:
            RG = -1
            I0 = -1
            RGer = -1
            I0er = -1

        return RG,I0,RGer,I0er,a,b

    '''
    Functions for basic Guiner analysis
    '''

    def lsq_w_sigma(self,X, Y, SIG):
        """ This is a special version of linear least squares that uses
        known sigma values of Y to calculate standard error in slope and
        intercept. The actual fit is also SIGMA-weighted"""
        XSIG2 = np.sum(X / (SIG ** 2))
        X2SIG2 = np.sum((X ** 2) / (SIG ** 2))
        INVSIG2 = np.sum(1.0 / (SIG ** 2))

        XYSIG2 = np.sum(X * Y / (SIG ** 2))
        YSIG2 = np.sum(Y / (SIG ** 2))

        delta = INVSIG2 * X2SIG2 - XSIG2 ** 2

        sby = X2SIG2 / delta
        smy = INVSIG2 / delta

        by = (X2SIG2 * YSIG2 - XSIG2 * XYSIG2) / delta
        my = (INVSIG2 * XYSIG2 - XSIG2 * YSIG2) / delta
        # outputs slope(my), intercept(by), sigma slope(smy), sigma intercept (sby)
        return my, by, smy, sby

    def Guiner_Error(self,q,I,I_Err,nmin,nmax,file='No File Description Provided', plot=True,
                     darkmode=False):
        '''
        Describe basic function here
        '''
        print('--------------------------------------------------------------------------')
        hb_slope,hb_inter,hb_sigslope,hb_siginter = self.lsq_w_sigma(q[nmin:nmax] ** 2,
                                                                     np.log(I[nmin:nmax]),
                                                       I_Err[nmin:nmax] /I[nmin:nmax])
        hbI0=np.exp(hb_inter)
        hbRg = np.sqrt(-3*hb_slope)
        hbRg_Err = np.absolute(hbRg* np.sqrt(hb_sigslope)/(2*hb_slope))
        hb_qminRg = hbRg*q[nmin]
        hb_qmaxRg = hbRg*q[nmax]

        model=self.lineModel(q**2,hb_slope,hb_inter)

        print('Guiner Analysis:')
        print('The input file was: %s'%file)
        print('The Rg is approximately %.2f'%hbRg + ' ' + '+/- %.2f'%(hbRg_Err) + 'Angstroms')
        print('The guiner range is approximately: (qminRg, qmaxRg) - %.2f, %.2f' % (hb_qminRg,hb_qmaxRg))

        if plot==True:
            self.plots.vertical_stackPlot(X1=q[nmin:nmax] ** 2,Y1=np.log(I[nmin:nmax]),Y1err=np.log(I_Err[nmin:nmax]),X2=q[nmin:nmax]**2,Y2=model[nmin:nmax],
                                          ylabel1='ln(I(q))',ylabel2='Residuals (a.u.)',xlabel='q = $\\frac{4 \pi sin(\\theta)}{\\lambda}$ ($\\AA^{-1}$)',
                                          Label1='Expt',Label2='Model',darkmode=darkmode)

        return hbI0,hbRg,hbRg_Err,hb_qminRg,hb_qmaxRg,model


    '''
    More advanced "autoRG" calculations
    '''

    # @jit(nopython=True,cache=False,parallel=False)
    def autoRg(self,q,i,err,single_fit=False,error_weight=True,plot=True,
               lowQclip=0,output_suppress=False,saveLabel='AutoRg'):
        '''This function automatically calculates the radius of gyration and scattering intensity at zero angle
        from a given scattering profile. It roughly follows the method used by the autorg function in the atsas package'''

        if output_suppress is not True:
            print('--------------------------------------------------------------------------')
            print('################################################################')
            print('AutoRg calculations beginning..')
            print('We will cutoff %.2f points at the beginning of the profile'%lowQclip)

        # RM!
        warnings.filterwarnings("ignore",category=RuntimeWarning) # deals with dividing through by zero issues, which will generate a RuntimeWarning

        q = q[lowQclip:]
        i = i[lowQclip:]
        err = err[lowQclip:]
        qmin,qmax = (0,len(q))

        q = q[qmin:qmax]
        i = i[qmin:qmax]
        err = err[qmin:qmax]
        # rg,rger,i0,i0er,idx_min,idx_max = self.autoRg_inner(q,i,err,qmin,single_fit,error_weight)

        try:
            rg,rger,i0,i0er,idx_min,idx_max = self.autoRg_inner(q,i,err,qmin,single_fit,error_weight)
        except Exception:  # Catches unexpected numba errors, I hope
            print('Darn it, AUTORG did not work!')
            traceback.print_exc()
            rg = -1
            rger = -1
            i0 = -1
            i0er = -1
            idx_min = -1
            idx_max = -1

        '''
        Other possible errors..
        '''

        if rg==-1.00 and rger==-1.00:
            print('AutoRg failed... Take a look at the profile and try manually fitting the GuinerError() function.')
        else:
            if output_suppress is not True:
                print('autoRg predicts the Rg to be: %.2f +/- %.2f' % (rg,rger))
                print('with a Guiner region of (qminRg, qmaxRg): %.2f, %.2f' % (q[idx_min] * rg,q[idx_max] * rg))

        slope=(-rg**2)/3
        inter=np.log(i0)


        if plot==True:
            self.plots.vertical_stackPlot(X1=q[idx_min:idx_max] ** 2,Y1=np.log(i[idx_min:idx_max]),Y1err=np.log(err[idx_min:idx_max]),X2=q[idx_min:idx_max]**2,Y2=self.lineModel(q[idx_min:idx_max]**2,slope,inter),
                                          ylabel1='ln(I(q))',ylabel2='Residuals (a.u.)',xlabel='q$^{2}$ = $(\\frac{4 \pi sin(\\theta)}{\\lambda})^{2}$ ($\\AA^{-2}$)',
                                          Label1='Expt',Label2='Model',saveLabel=saveLabel)
        else:
            if output_suppress is not True:
                print('No plots were generated. Set plot=True to visualize the result.')

        if output_suppress is not True:
            print('################################################################')
            print('--------------------------------------------------------------------------')

        return rg,rger,i0,i0er,idx_min,idx_max

    # @jit(nopython=True,cache=False,parallel=False)
    def autoRg_inner(self,q,i,err,qmin,single_fit,error_weight):
        # Pick the start of the RG fitting range. Note that in autorg, this is done
        # by looking for strong deviations at low q from aggregation or structure factor
        # or instrumental scattering, and ignoring those. This function isn't that advanced
        # so we start at 0.

        # Note, in order to speed this up using numba, I had to do some unpythonic things
        # with declaring lists ahead of time, and making sure lists didn't have multiple
        # object types in them. It makes the code a bit more messy than the original
        # version, but numba provides a significant speedup.
        data_start = 0

        # Following the atsas package, the end point of our search space is the q value
        # where the intensity has droped by an order of magnitude from the initial value.
        data_end = np.abs(i - i[data_start] / 10.).argmin()

        # This makes sure we're not getting some weird fluke at the end of the scattering profile.
        if data_end > len(q) / 2.:
            found = False
            idx = 0
            while not found:
                idx = idx + 1
                if i[idx] < i[0] / 10.:
                    found = True
                elif idx == len(q) - 1:
                    found = True
            data_end = idx

        # Start out by transforming as usual.
        qs = np.square(q)
        il = np.log(i)
        iler = np.absolute(err / i)

        # Pick a minimum fitting window size. 10 is consistent with atsas autorg.
        min_window = 10

        max_window = data_end - data_start

        if max_window < min_window:
            max_window = min_window

        # It is very time consuming to search every possible window size and every possible starting point.
        # Here we define a subset to search.
        tot_points = max_window
        window_step = tot_points // 10
        data_step = tot_points // 50

        if window_step == 0:
            window_step = 1
        if data_step == 0:
            data_step = 1

        window_list = [0 for k in range(int(math.ceil((max_window - min_window) / float(window_step))) + 1)]

        for k in range(int(math.ceil((max_window - min_window) / float(window_step)))):
            window_list[k] = min_window + k * window_step

        window_list[-1] = max_window

        num_fits = 0

        for w in window_list:
            num_fits = num_fits + int(math.ceil((data_end - w - data_start) / float(data_step)))

        if num_fits < 0:
            num_fits = 1

        ''' 
        Initializing a bunch of lists..
        '''

        start_list = [0 for k in range(num_fits)]
        w_list = [0 for k in range(num_fits)]
        q_start_list = [0. for k in range(num_fits)]
        q_end_list = [0. for k in range(num_fits)]
        rg_list = [0. for k in range(num_fits)]
        rger_list = [0. for k in range(num_fits)]
        i0_list = [0. for k in range(num_fits)]
        i0er_list = [0. for k in range(num_fits)]
        qrg_start_list = [0. for k in range(num_fits)]
        qrg_end_list = [0. for k in range(num_fits)]
        rsqr_list = [0. for k in range(num_fits)]
        chi_sqr_list = [0. for k in range(num_fits)]
        reduced_chi_sqr_list = [0. for k in range(num_fits)]

        success = np.zeros(num_fits)

        current_fit = 0
        # This function takes every window size in the window list, stepts it through the data range, and
        # fits it to get the RG and I0. If basic conditions are met, qmin*RG<1 and qmax*RG<1.35, and RG>0.1,
        # We keep the fit.
        for w in window_list:
            for start in range(data_start,data_end - w,data_step):
                x = qs[start:start + w]
                y = il[start:start + w]
                yerr = iler[start:start + w]

                # Remove NaN and Inf values:
                x = x[np.where(np.isfinite(y))]
                yerr = yerr[np.where(np.isfinite(y))]
                y = y[np.where(np.isfinite(y))]

                qmaxRgLimit=1.55

                RG,I0,RGer,I0er,a,b = self.calcRg(x,y,yerr,transform=False,error_weight=error_weight)

                if RG > 0.1 and q[start] * RG < 1 and q[start + w - 1] * RG < qmaxRgLimit and RGer / RG <= 1:

                    r_sqr = 1 - np.square(il[start:start + w] - self.linear_func(qs[start:start + w],a,b)).sum() / np.square(
                        il[start:start + w] - il[start:start + w].mean()).sum()

                    if r_sqr > .15:
                        chi_sqr = np.square(
                            (il[start:start + w] - self.linear_func(qs[start:start + w],a,b)) / iler[start:start + w]).sum()

                        # All of my reduced chi_squared values are too small, so I suspect something isn't right with that.
                        # Values less than one tend to indicate either a wrong degree of freedom, or a serious overestimate
                        # of the error bars for the system.
                        dof = w - 2.
                        reduced_chi_sqr = chi_sqr / dof

                        start_list[current_fit] = start
                        w_list[current_fit] = w
                        q_start_list[current_fit] = q[start]
                        q_end_list[current_fit] = q[start + w - 1]
                        rg_list[current_fit] = RG
                        rger_list[current_fit] = RGer
                        i0_list[current_fit] = I0
                        i0er_list[current_fit] = I0er
                        qrg_start_list[current_fit] = q[start] * RG
                        qrg_end_list[current_fit] = q[start + w - 1] * RG
                        rsqr_list[current_fit] = r_sqr
                        chi_sqr_list[current_fit] = chi_sqr
                        reduced_chi_sqr_list[current_fit] = reduced_chi_sqr

                        # fit_list[current_fit] = [start, w, q[start], q[start+w-1], RG, RGer, I0, I0er, q[start]*RG, q[start+w-1]*RG, r_sqr, chi_sqr, reduced_chi_sqr]
                        success[current_fit] = 1 # success was an empty array of zeroes

                current_fit = current_fit + 1 # starts at 0, adds a value each time we iterate through the loop.
        # Extreme cases: may need to relax the parameters.
        if np.sum(success) == 0:
            # Stuff goes here
            pass

        if np.sum(success) > 0:

            fit_array = np.array([[start_list[k],w_list[k],q_start_list[k],
                                   q_end_list[k],rg_list[k],rger_list[k],i0_list[k],i0er_list[k],
                                   qrg_start_list[k],qrg_end_list[k],rsqr_list[k],chi_sqr_list[k],
                                   reduced_chi_sqr_list[k]] for k in range(num_fits) if success[k] == 1])

            # Now we evaluate the quality of the fits based both on fitting data and on other criteria.

            # Choice of weights is pretty arbitrary. This set seems to yield results similar to the atsas autorg
            # for the few things I've tested.
            qmaxrg_weight = 1
            qminrg_weight = 1
            rg_frac_err_weight = 1
            i0_frac_err_weight = 1
            r_sqr_weight = 4
            reduced_chi_sqr_weight = 0
            window_size_weight = 4

            weights = np.array([qmaxrg_weight,qminrg_weight,rg_frac_err_weight,i0_frac_err_weight,r_sqr_weight,
                                reduced_chi_sqr_weight,window_size_weight])

            quality = np.zeros(len(fit_array))

            max_window_real = float(window_list[-1])

            # all_scores = [np.array([]) for k in range(len(fit_array))]

            # This iterates through all the fits, and calculates a score. The score is out of 1, 1 being the best, 0 being the worst.
            indices = list(range(len(fit_array)))
            for a in indices:
                k = int(a)  # This is stupid and should not be necessary. Numba bug?

                # Scores all should be 1 based. Reduced chi_square score is not, hence it not being weighted.
                qmaxrg_score = 1 - abs((fit_array[k,9] - 1.3) / 1.3)
                qminrg_score = 1 - fit_array[k,8]
                rg_frac_err_score = 1 - fit_array[k,5] / fit_array[k,4]
                i0_frac_err_score = 1 - fit_array[k,7] / fit_array[k,6]
                r_sqr_score = fit_array[k,10]
                reduced_chi_sqr_score = 1 / fit_array[k,12]  # Not right
                window_size_score = fit_array[k,1] / max_window_real

                scores = np.array([qmaxrg_score,qminrg_score,rg_frac_err_score,i0_frac_err_score,r_sqr_score,
                                   reduced_chi_sqr_score,window_size_score])

                total_score = (weights * scores).sum() / weights.sum()

                quality[k] = total_score

                # all_scores[k] = scores

            # I have picked an aribtrary threshold here. Not sure if 0.6 is a good quality cutoff or not.
            if quality.max() > 0.50: # RM!
                if not single_fit:
                    idx = quality.argmax()
                    rg = fit_array[:,4][quality > quality[idx] - .1].mean()
                    rger = fit_array[:,5][quality > quality[idx] - .1].std()
                    i0 = fit_array[:,6][quality > quality[idx] - .1].mean()
                    i0er = fit_array[:,7][quality > quality[idx] - .1].std()
                    idx_min = int(fit_array[idx,0])
                    idx_max = int(fit_array[idx,0] + fit_array[idx,1] - 1)
                else:
                    idx = quality.argmax()
                    rg = fit_array[idx,4]
                    rger = fit_array[idx,5]
                    i0 = fit_array[idx,6]
                    i0er = fit_array[idx,7]
                    idx_min = int(fit_array[idx,0])
                    idx_max = int(fit_array[idx,0] + fit_array[idx,1] - 1)

            else:
                rg = -1
                rger = -1
                i0 = -1
                i0er = -1
                idx_min = -1
                idx_max = -1

        else:
            rg = -1
            rger = -1
            i0 = -1
            i0er = -1
            idx_min = -1
            idx_max = -1
            # quality = []
            # all_scores = []

        idx_min = idx_min + qmin
        idx_max = idx_max + qmin

        # We could add another function here, if not good quality fits are found, either reiterate through the
        # the data and refit with looser criteria, or accept lower scores, possibly with larger error bars.

        # returns Rg, Rg error, I0, I0 error, the index of the first q point of the fit and the index of the last q point of the fit
        return rg,rger,i0,i0er,idx_min,idx_max

    def runGNOM(self,file_path,output_name, gnom_nmin=50,
        rmax='Auto',force_zero_rmin=None,force_zero_rmax=None,system=0,radius56=None,alpha=None,
        plot=True,darkmode=False):
        '''
        Need to fetch the following:
        (1) ATSAS - GNOM directory
        (2) Input file directory
        (3) GNOM Parameters (A dictionary of additional arguments to provide)
            i.e. if an entry is empty, skip it, otherwise append it to the string of commands.
            output_name = output. Must be input as a string
            rmax: must be input as an integer or float. Default 100 angstroms.
            force_zero_rmax: Must be set to 'yes' or 'no'

        Output file goes into the current working directory:
        To determine what that is run:
            > import os
            > print(os.getcwd())

        rmax: options -- (1) 'Auto': AutoRg tries to determine Rg

        '''

        '''
        Consider that the command
        > which dammif SORT OF tells us where ATSAS is located.. I could probably manipulate that to find the directory automatically.
        '''

        print('--------------------------------------------------------------------------')
        print('################################################################')
        print('P(r) calculations beginning')

        atsas_dir=self.atsas_dir
        if atsas_dir == '':
            print('Cannot run "runGNOM" function without knowing the ATSAS directory path')


        def concatenate_list_data(list):
            result= ''
            for element in list:
                result += str(element)
            return result

        if rmax=='Auto':
            print('Attempting to determine Rg using autoRg function on input file...')
            data = np.loadtxt(file_path,dtype={'names':('Q','I(Q)','ERROR'),'formats':(np.float,np.float,np.float)},skiprows=4)
            auto_rg_data = self.autoRg(q=data['Q'],i=data['I(Q)'],err=data['ERROR'],output_suppress=True)
            rmax= str(auto_rg_data[0]) + ' ' # fetches Rg from autoRg
            rmax_print=auto_rg_data[0]
            print('It worked! The input Rg is will be %.2f'%rmax_print)
        else:
            rmax = str(rmax) + ' '


        if radius56==None:
            radius56=radius56
        else:
            radius56=radius56 + ' '

        if force_zero_rmin==None:
            force_zero_rmin=force_zero_rmin
        else:
            force_zero_rmin=str(force_zero_rmin) + ' '

        if force_zero_rmax==None:
            force_zero_rmax=force_zero_rmax
        else:
            force_zero_rmax=str(force_zero_rmax) + ' '

        if alpha==None:
            alpha=alpha
        else:
            alpha=alpha + ' '


        '''
        Dealing with specific exceptions

        I do need to be a bit careful this may suppress all traceback errors once it is called.
        '''

        def exception_handler(exception_type, exception, traceback):
            # All your trace are belong to us!
            # your format
            print ("%s: %s" % (exception_type.__name__, exception))

        sys.excepthook = exception_handler

        if system==2:
            raise Exception('You cannot run --system=2 in command line mode')

        system=str(system) + ' '

        '''
        Setting up the output. Puts the file in the current working directory.
        '''

        output=str(os.getcwd()) + '/' + output_name # dump output file in current working directory.
        output=concatenate_list_data(output)
        output=(f'"{output}"') # ensures double quotes around the outfile name which is a requirement for the GNOM input
        

        '''
        Build GNOM command inputs
        '''


        GNOM = {
        '--force-zero-rmin=':force_zero_rmin,
        '--force-zero-rmax=':force_zero_rmax,
        '--rmax=':rmax,
        '--system=':system,
        '--rad56=':radius56,
        '--alpha=':alpha,
        '--output ':output,
        }

        '''
        If no argument is provided for a GNOM input by default it is set to None.
        If a value in the dictionary if None, we will remove it.
        That way it won't be passed into the final command.
        '''

        filtered = {k: v for k, v in GNOM.items() if v is not None}
        GNOM.clear()
        GNOM.update(filtered)


        GNOM_params=[]
        for key,value in GNOM.items():
            GNOM_params.append(key + value)

        GNOM_params=concatenate_list_data(GNOM_params)


        '''
        We need to build a way to cleave points at the beginning of the input file.
        Currently, I'll do it in a hacky way. We'll import the data file.
        Chop off what we want, create a new data file, and read that in!
        Remember, this file that it temporarily creates will get over-written each time you run this command!
        '''

        scatteringData=np.loadtxt(file_path,
            dtype={'names': ('Q', 'I(Q)','ERROR'), 'formats': (np.float,np.float,np.float)},
            skiprows=4)



        print('The length of the dataframe before cleaving: %s'%(str(len(scatteringData['Q']))))
        scatteringData=scatteringData[gnom_nmin:]
        entryCount=0
        print('The length of the dataframe after cleaving: %s'%(str(len(scatteringData['Q']))))

        with open('gnom_import.dat','w',newline='') as csvfile:
            fieldnames=['Q', 'I(Q)','ERROR']
            thewriter=csv.DictWriter(csvfile,fieldnames=fieldnames)
            thewriter.writeheader()
            n=len(scatteringData['Q'])
            for i in range(n):
                entryCount += 1
                thewriter.writerow({'Q':scatteringData['Q'][i],'I(Q)':scatteringData['I(Q)'][i],'ERROR':scatteringData['ERROR'][i]})


        file_path = str(os.getcwd()) + '/' + 'gnom_import.dat'
        # print(file_path)

        '''
        Finally, generate the command we will pass to GNOM
        '''

        cmd =  'cd ;'
        cmd = cmd + atsas_dir + ' ' + file_path + ' ' + GNOM_params

        print('################################################################')
        print('The final GNOM command that was passed through your terminal is:')
        n=int(len(cmd)/2)
        print(cmd[:n])
        print(cmd[n:])
        print('################################################################')
        # execute the GNOM command

        '''
        We will pass in the UNIX command with try and except python commands.
        This will allow us to better handle potential issues that will arise
        when attempting to run the calculation.
        '''

        try:
            os.system(cmd)
        except FileNotFoundError: # specific exceptions first!
            print('Could not locate the input file for runGNOM function')
        except Exception as e:
            print('Sorry, something went wrong that I had not anticipated', e)
        else:
            print('GNOM calculation execution was successful')
        finally:
            print('Now moving forward with downstream processing!')


        '''
        Now need to read the file that it just generated and report the important parameters:
        Dmax
        Beautiful Plot

        use filepath to find output file, with name = output
        but must be able to parse GNOM output file..
        '''
        IFT_input=str(os.getcwd()) + '/' + output_name
        IFT=self.FileParser.loadOutFile(str(IFT_input))
        Pr, R, Pr_err, Jexp, qshort, Jerr, Jreg, results, Ireg, qfull =IFT[0],IFT[1],IFT[2],IFT[3],IFT[4],IFT[5],IFT[6],IFT[7],IFT[8],IFT[9]

        Pr_norm=self.calcs.quickNormalize(Pr)

        output_dict = {'Pr':Pr,'Pr_norm':Pr_norm,'R':R,'Pr_err':Pr_err,'Jexp':Jexp,'qshort':qshort,'Jerr':Jerr,'Jreg':Jreg,
                       'results':results,'Ireg':Ireg,'qfull':qfull}
        print('GNOM reports the quality of the regularized fit to be: %s'%results['quality'])
        print('From GNOM calculated P(r) the Rg is reported as: %.2f +/- %.2f'%(results['rg'],results['rger']))
        print('From GNOM calculated P(r) the I0 is reported as: %.2f +/- %.2f' % (results['i0'],results['i0er']))
        print('From GNOM calculated P(r) the Dmax is reported as: %.2f'%results['dmax'])

        if plot==True:
            self.plots.twoPlot(X=R,Y1=Pr,Y2=[0]*len(Pr),savelabel=output_name+'_PDDF',
                plotlabel1='Pair Distance Distribution',plotlabel2='Baseline',
                    xlabel='r($\\AA$)',ylabel='P(r)',linewidth=4,darkmode=darkmode)
            self.plots.twoPlot_variX(X1=qshort,Y1=Jexp, X2=qshort,Y2=Jreg,plotlabel1='Expt',plotlabel2='Regularized Fit',
                               savelabel=output_name+'RegularizedFit_GNOM',xlabel='q $\\AA^{-1}$',ylabel='I(q)',LogLin=True,darkmode=darkmode)
        else:
            print('We did not plot the PDDF. If you want to see the PDDF, set plot=True in the runGNOM() arguments.')


        print('--------------------------------------------------------------------------')
        return output_dict


    def runDatgnom(self, file_path, save_path, outname,rg='Auto', first_pt=None, last_pt=None,plot=True,
                   lowclip=0):
        '''
        Adopted from RAW2.0.2
        Parameters
        ----------
        rg: input hydrodynamic radius calculated from GuinerError() function.
        file_path
        save_path
        outname
        first_pt
        last_pt

        Returns
        -------

        '''
        # This runs the ATSAS package DATGNOM program, to automatically find the Dmax and P(r) function
        # of a scattering profile.

        print('--------------------------------------------------------------------------')
        print('###################################################################')
        print('DATGNOM P(r) calculations beginning')



        '''
        If no save path is provided, dump the file in the current working directory.
        '''

        if save_path=='':
            save_path=os.getcwd()

        '''
        Building some extra compatibility if running on a windows computer.
        '''

        opsys = platform.system()

        if opsys == 'Windows':
            datgnomDir = os.path.join(self.atsas_dir.replace('gnom','datgnom.exe'))
            shell = False
        else:
            datgnomDir = os.path.join(self.atsas_dir.replace('gnom','datgnom'))
            shell = True

        '''
        Building the command to pass to DATGNOM
        And adding AutoRG capabilities! I've tried to build it to catch most errors.. 
        again, slightly modified from RAW code.
        '''
        if rg=='Auto':
            print('Attempting to determine Rg using autoRg function on input file...')
            try:
                data = np.loadtxt(file_path,dtype={'names':('Q','I(Q)','ERROR'),'formats':(np.float,np.float,np.float)},skiprows=4)
                auto_rg_data = self.autoRg(q=data['Q'][lowclip:],i=data['I(Q)'][lowclip:],err=data['ERROR'][lowclip:],output_suppress=False)
                rg= auto_rg_data[0] # fetches Rg from autoRg
                if rg == -1.00:
                    print('Terminating analysis..')
                    sys.exit()
            except Exception as e:
                print('Terminating analysis..')
                print(e)
                sys.exit()
        else:
            rg = rg

        if os.path.exists(datgnomDir):
            cmd = 'cd; ' + '"{}" -o "{}" -r {} '.format(datgnomDir,save_path+outname,rg)

            if first_pt is not None:
                cmd = cmd + '--first={} '.format(first_pt + 1)

            if last_pt is not None:
                cmd = cmd + ' --last={} '.format(last_pt + 1)

            cmd = cmd + '"{}"'.format(file_path)

            print('###################################################################')
            print('The final DATGNOM command that was passed through your terminal is:')
            n = int(len(cmd) / 2)
            print(cmd[:n])
            print(cmd[n:])
            print('###################################################################')

            process = subprocess.Popen(cmd,stdout=subprocess.PIPE,
                                       stderr=subprocess.PIPE,shell=shell,cwd=os.getcwd())

            output,error = process.communicate()


            if not isinstance(output,str):
                output = str(output,encoding='UTF-8')

            if not isinstance(error,str):
                error = str(error,encoding='UTF-8')

            error = error.strip()

            # print(process.communicate())

            if (error == 'Cannot define Dmax' or error == 'Could not find Rg'
                    or error == 'No intensity values (positive) found'
                    or error == 'LOADATF --E- No data lines recognized.'
                    or error == 'error: rg not specified'
                    or 'error' in error):
                datgnom_success = False
            else:
                datgnom_success = True

            print('The final output file/filepath is:\n%s'%os.path.join(os.getcwd(),outname))
            # print(os.path.join(os.getcwd(),outname))

            if datgnom_success:
                try:
                    iftm = self.FileParser.loadOutFile(os.path.join(save_path,outname))
                except Exception:
                    iftm = None
            else:
                iftm = None

            Pr,R,Pr_err,Jexp,qshort,Jerr,Jreg,results,Ireg,qfull = iftm[0],iftm[1],iftm[2],iftm[3],iftm[4],iftm[5],iftm[6],iftm[7],iftm[8],iftm[9]
            Pr_norm = self.calcs.quickNormalize(Pr)

            output_dict = {'Pr':Pr,'Pr_norm':Pr_norm,'R':R,'Pr_err':Pr_err,'Jexp':Jexp,'qshort':qshort,'Jerr':Jerr,
                           'Jreg':Jreg,'results':results,'Ireg':Ireg,'qfull':qfull}


            print('DATGNOM reports the quality of the regularized fit to be: %s' % results['quality'])
            print('From DATGNOM calculated P(r) the Rg is reported as: %.2f +/- %.2f' % (results['rg'],results['rger']))
            print('From DATGNOM calculated P(r) the I0 is reported as: %.2f +/- %.2f' % (results['i0'],results['i0er']))
            print('From DATGNOM calculated P(r) the Dmax is reported as: %.2f' % results['dmax'])

            if plot == True:
                self.plots.twoPlot(X=R,Y1=Pr,Y2=[0] * len(Pr),savelabel=outname+'_DatGNOM_PDDF',
                                   plotlabel1='Pair Distance Distribution',plotlabel2='Baseline',
                                   xlabel='r($\\AA$)',ylabel='P(r)',linewidth=4)
                self.plots.twoPlot_variX(X1=qshort,Y1=Jexp,X2=qshort,Y2=Jreg,plotlabel1='Expt',
                                         plotlabel2='Regularized Fit',
                                         savelabel=outname+'RegularizedFit_DATGNOM',xlabel='q $\\AA^{-1}$',ylabel='I(q)',
                                         LogLin=True)
            else:
                print('We did not plot the PDDF. If you want to see the PDDF, set plot=True in the runDatgnom() arguments.')

            print('###################################################################')
            print('Given these results, you should attempt to manually refine the P(r) using the runGNOM() function available in this class.')
            print('--------------------------------------------------------------------------')
            return iftm,output_dict

        else:
            print('Cannot find ATSAS')
            raise Exception('Cannot find datgnom.')
            
            

    # def peakScalar(self,peakList):
    #     '''
    #     Input: list of N column vectors containing spectrumy like data to scale to the maximum peak
        
    #     output: a matrix containing N scaled column vectors
    #     '''
    #     # How do I want to normalize??
    #     # First need to align peaks? or just by total area?


    def PDDF(self,shape,Dmax,I,q):
        """plot pair distance distribution
        This function is now essentially useless now that I have GNOM running in runGNOM()..
        """
        # reference:"Svergen&Koch,Rep.Phys.Prog 66(2003) 1735-82
        r_range = np.arange(0,Dmax * 1.4,Dmax / 100)
        if shape == "FoxS":
            if r_range[0] == 0: # RM! first point in r_range is 0, this will create issues.
                q = q[1:]
                I = I[1:]
                r_range=r_range[1:]
                # Taking the first point in exp_q out if it's 0, avoiding dividing by 0 problem
            else:
                q = q

        '''
        Initialize an array to stick the data in
        Note: this MUST be outside of the loop, otherwise we will clear the array for each iteration of the loop.
        '''
        P_r = np.array([],dtype=float)


        for r in r_range:
            p_r = np.sum(q ** 2 * I * np.sin(q * r) / (q * r) * 0.02) * (r ** 2) / (2.0 * np.pi ** 2)
            P_r = np.append(P_r,p_r)

        # print(q[0],P_r[0],r_range[0])

        return P_r,r_range


    def nanCheck(self,data):
        for entry in data:
            c=0
            if np.isnan(entry) == True:
                print("Uh oh, we've generated an nan value")
                nanIndex = np.where(entry)
                print("The index of the nan value is: %s" % str(nanIndex[0]))
                c+=1
        if c == 0:
            print('No nan values generated')
