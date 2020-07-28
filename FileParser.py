import sys
import os 
import re
import numpy as np  

class FileParser():


    def __init__(self,notify=True):
        '''
        Describe class here.
        '''
        if notify==True:
            print('---------------------------')
            print('FileParser Class was called')
            print('---------------------------')



    def loadOutFile(self,filename):
        '''
        This function reads the output file from the GNOM P(R) calculations and was adopted from
        bioxtsas RAW2.0.2 (!REF)
        '''

        five_col_fit = re.compile('\s*-?\d*[.]\d*[+eE-]*\d+\s+-?\d*[.]\d*[+eE-]*\d+\s+-?\d*[.]\d*[+eE-]*\d+\s+-?\d*[.]\d*[+eE-]*\d+\s+-?\d*[.]\d*[+eE-]*\d+\s*$')
        three_col_fit = re.compile('\s*-?\d*[.]\d*[+eE-]*\d+\s+-?\d*[.]\d*[+eE-]*\d+\s+-?\d*[.]\d*[+eE-]*\d+\s*$')
        two_col_fit = re.compile('\s*-?\d*[.]\d*[+eE-]*\d+\s+-?\d*[.]\d*[+eE-]*\d+\s*$')

        results_fit = re.compile('\s*Current\s+\d*[.]\d*[+eE-]*\d*\s+\d*[.]\d*[+eE-]*\d*\s+\d*[.]\d*[+eE-]*\d*\s+\d*[.]\d*[+eE-]*\d*\s+\d*[.]\d*[+eE-]*\d*\s+\d*[.]\d*[+eE-]*\d*\s*\d*[.]?\d*[+eE-]*\d*\s*$')

        te_fit = re.compile('\s*Total\s+[Ee]stimate\s*:\s+\d*[.]\d+\s*\(?[A-Za-z\s]+\)?\s*$')
        te_num_fit = re.compile('\d*[.]\d+')
        te_quality_fit = re.compile('[Aa][A-Za-z\s]+\)?\s*$')

        p_rg_fit = re.compile('\s*Real\s+space\s*\:?\s*Rg\:?\s*\=?\s*\d*[.]\d+[+eE-]*\d*\s*\+-\s*\d*[.]\d+[+eE-]*\d*')
        q_rg_fit = re.compile('\s*Reciprocal\s+space\s*\:?\s*Rg\:?\s*\=?\s*\d*[.]\d+[+eE-]*\d*\s*')

        p_i0_fit = re.compile('\s*Real\s+space\s*\:?[A-Za-z0-9\s\.,+-=]*\(0\)\:?\s*\=?\s*\d*[.]\d+[+eE-]*\d*\s*\+-\s*\d*[.]\d+[+eE-]*\d*')
        q_i0_fit = re.compile('\s*Reciprocal\s+space\s*\:?[A-Za-z0-9\s\.,+-=]*\(0\)\:?\s*\=?\s*\d*[.]\d+[+eE-]*\d*\s*')

        alpha_fit = re.compile('\s*Current\s+ALPHA\s*\:?\s*\=?\s*\d*[.]\d+[+eE-]*\d*\s*')

        qfull = []
        qshort = []
        Jexp = []
        Jerr  = []
        Jreg = []
        Ireg = []

        R = []
        P = []
        Perr = []

        outfile = []

        #In case it returns NaN for either value, and they don't get picked up in the regular expression
        q_rg=None         #Reciprocal space Rg
        q_i0=None         #Reciprocal space I0

        #Set some defaults in case the .out file isn't perfect. I've encountered
        #at least one case where no DISCRIP is returned, which messes up loading in
        #the .out file.
        Actual_DISCRP = -1
        Actual_OSCILL = -1
        Actual_STABIL = -1
        Actual_SYSDEV = -1
        Actual_POSITV = -1
        Actual_VALCEN = -1
        Actual_SMOOTH = -1


        with open(filename, 'rU') as f:
            for line in f:
                twocol_match = two_col_fit.match(line)
                threecol_match = three_col_fit.match(line)
                fivecol_match = five_col_fit.match(line)
                results_match = results_fit.match(line)
                te_match = te_fit.match(line)
                p_rg_match = p_rg_fit.match(line)
                q_rg_match = q_rg_fit.match(line)
                p_i0_match = p_i0_fit.match(line)
                q_i0_match = q_i0_fit.match(line)
                alpha_match = alpha_fit.match(line)

                outfile.append(line)

                if twocol_match:
                    # print line
                    found = twocol_match.group().split()

                    qfull.append(float(found[0]))
                    Ireg.append(float(found[1]))

                elif threecol_match:
                    #print line
                    found = threecol_match.group().split()

                    R.append(float(found[0]))
                    P.append(float(found[1]))
                    Perr.append(float(found[2]))

                elif fivecol_match:
                    #print line
                    found = fivecol_match.group().split()

                    qfull.append(float(found[0]))
                    qshort.append(float(found[0]))
                    Jexp.append(float(found[1]))
                    Jerr.append(float(found[2]))
                    Jreg.append(float(found[3]))
                    Ireg.append(float(found[4]))

                elif results_match:
                    found = results_match.group().split()
                    Actual_DISCRP = float(found[1])
                    Actual_OSCILL = float(found[2])
                    Actual_STABIL = float(found[3])
                    Actual_SYSDEV = float(found[4])
                    Actual_POSITV = float(found[5])
                    Actual_VALCEN = float(found[6])

                    if len(found) == 8:
                        Actual_SMOOTH = float(found[7])

                elif te_match:
                    te_num_search = te_num_fit.search(line)
                    te_quality_search = te_quality_fit.search(line)

                    TE_out = float(te_num_search.group().strip())
                    quality = te_quality_search.group().strip().rstrip(')').strip()


                if p_rg_match:
                    found = p_rg_match.group().split()
                    try:
                        rg = float(found[-3])
                    except:
                        rg = float(found[-2])
                    try:
                        rger = float(found[-1])
                    except:
                        rger = float(found[-1].strip('+-'))

                elif q_rg_match:
                    found = q_rg_match.group().split()
                    q_rg = float(found[-1])

                if p_i0_match:
                    found = p_i0_match.group().split()
                    i0 = float(found[-3])
                    i0er = float(found[-1])

                elif q_i0_match:
                    found = q_i0_match.group().split()
                    q_i0 = float(found[-1])

                if alpha_match:
                    found = alpha_match.group().split()
                    alpha = float(found[-1])

        # Output variables not in the results file:
            # 'r'         : R,            #R, note R[-1] == Dmax
            # 'p'         : P,            #P(r)
            # 'perr'      : Perr,         #P(r) error
            # 'qlong'     : qfull,        #q down to q=0
            # 'qexp'      : qshort,       #experimental q range
            # 'jexp'      : Jexp,         #Experimental intensities
            # 'jerr'      : Jerr,         #Experimental errors
            # 'jreg'      : Jreg,         #Experimental intensities from P(r)
            # 'ireg'      : Ireg,         #Experimental intensities extrapolated to q=0

        name = os.path.basename(filename)

        chisq = np.sum(np.square(np.array(Jexp)-np.array(Jreg))/np.square(Jerr))/(len(Jexp)-1) #DOF normalied chi squared

        results = { 'dmax'      : R[-1],        #Dmax
                    'TE'        : TE_out,       #Total estimate
                    'rg'        : rg,           #Real space Rg
                    'rger'      : rger,         #Real space rg error
                    'i0'        : i0,           #Real space I0
                    'i0er'      : i0er,         #Real space I0 error
                    'q_rg'      : q_rg,         #Reciprocal space Rg
                    'q_i0'      : q_i0,         #Reciprocal space I0
                    'out'       : outfile,      #Full contents of the outfile, for writing later
                    'quality'   : quality,      #Quality of GNOM out file
                    'discrp'    : Actual_DISCRP,#DISCRIP, kind of chi squared (normalized by number of points, with a regularization parameter thing thrown in)
                    'oscil'     : Actual_OSCILL,#Oscillation of solution
                    'stabil'    : Actual_STABIL,#Stability of solution
                    'sysdev'    : Actual_SYSDEV,#Systematic deviation of solution
                    'positv'    : Actual_POSITV,#Relative norm of the positive part of P(r)
                    'valcen'    : Actual_VALCEN,#Validity of the chosen interval in real space
                    'smooth'    : Actual_SMOOTH,#Smoothness of the chosen interval? -1 indicates no real value, for versions of GNOM < 5.0 (ATSAS <2.8)
                    'filename'  : name,         #GNOM filename
                    'algorithm' : 'GNOM',       #Lets us know what algorithm was used to find the IFT
                    'chisq'     : chisq,        #Actual chi squared value
                    'alpha'     : alpha,        #Alpha used for the IFT
                    'qmin'      : qshort[0],    #Minimum q
                    'qmax'      : qshort[-1],   #Maximum q
                        }

        # iftm = SASM.IFTM(P, R, Perr, Jexp, qshort, Jerr, Jreg, results, Ireg, qfull)

        # return [iftm]
        iftm = P, R, Perr, Jexp, qshort, Jerr, Jreg, results, Ireg, qfull

        return iftm

