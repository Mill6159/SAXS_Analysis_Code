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
import math
from PlotClass import *
from FileParser import *

#######################################################
# Define Class
#######################################################

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
                print('The GNOM library was found (!!!) in the directory: %s. If this is the incorrect library, the calculation may fail'%gnom_dir)
            print('##############################################################')
        else:
            print('#######################')
            print('GNOM directory provided')
            print('#######################')
            print('--------------------------------------------------------------')

        self.plots=PlotClass(notify=False)
        self.FileParser=FileParser(notify=False)

    def lineModel(self,x,m,b):
        return m*x+b

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

    def Guiner_Error(self,q,I,I_Err,nmin,nmax,file='No File Description Provided'):
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



        return hbI0,hbRg,hbRg_Err,hb_qminRg,hb_qmaxRg,model


    def runGNOM(self,file_path,output_name, gnom_nmin=50,
        rmax=100,force_zero_rmin=None,force_zero_rmax=None,system=0,radius56=None,alpha=None,
        plot=True):
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



        # rmax=str(120) + ' ' # comment out once finished building
        rmax=str(rmax) + ' '
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

        print('GNOM reports the quality of the regularized fit to be: %s'%results['quality'])
        print('From GNOM calculated P(r) the Rg is reported as: %.2f +/- %.2f'%(results['rg'],results['rger']))
        print('From GNOM calculated P(r) the I0 is reported as: %.2f +/- %.2f' % (results['i0'],results['i0er']))
        print('From GNOM calculated P(r) the Dmax is reported as: %.2f'%results['dmax'])

        if plot==True:
            self.plots.twoPlot(X=R,Y1=Pr,Y2=[0]*len(Pr),savelabel='tkRubisCO_0MPa_GNOM_PDDF',
                plotlabel1='Pair Distance Distribution',plotlabel2='Baseline',
                    xlabel='r($\\AA$)',ylabel='P(r)',linewidth=4)
            self.plots.twoPlot_variX(X1=qshort,Y1=Jexp, X2=qshort,Y2=Jreg,plotlabel1='Expt',plotlabel2='Regularized Fit',
                               savelabel='RegularizedFit_GNOM',xlabel='q $\\AA^{-1}$',ylabel='I(q)',LogLin=True)
        else:
            print('We did not plot the PDDF. If you want to see the PDDF, set plot=True in the runGNOM() arguments.')


        print('--------------------------------------------------------------------------')

    def runDatgnom(self,rg, file_path, save_path, outname, first_pt=None, last_pt=None,plot=True):
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
        Add an AUTORG functionaility
        '''

        # if rg=='':
        #     hbI0,Rg,hbRg_Err,hb_qminRg,hb_qmaxRg,model=Guiner_Error(q,I,I_Err,nmin,nmax,file='No File Description Provided')

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
        '''

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

            print('DATGNOM reports the quality of the regularized fit to be: %s' % results['quality'])
            print('From DATGNOM calculated P(r) the Rg is reported as: %.2f +/- %.2f' % (results['rg'],results['rger']))
            print('From DATGNOM calculated P(r) the I0 is reported as: %.2f +/- %.2f' % (results['i0'],results['i0er']))
            print('From DATGNOM calculated P(r) the Dmax is reported as: %.2f' % results['dmax'])

            if plot == True:
                self.plots.twoPlot(X=R,Y1=Pr,Y2=[0] * len(Pr),savelabel=outname,
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
            return iftm

        else:
            print('Cannot find ATSAS')
            raise Exception('Cannot find datgnom.')



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
