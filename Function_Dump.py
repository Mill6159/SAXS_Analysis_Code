### REF line 2223 from RAWAPI.py from RAW for how to build this yourself.
### It requires quite of bit of backend to get this bad boy to run...

from numba import jit
import math


### 07/23/20 We are far far far away from complete..

class GNOMError(Exception):
    def __init__(self, value):
       self.parameter = value
    def __str__(self):
       return repr(self.parameter)


class NoATSASError(Exception):
    def __init__(self, value):
       self.parameter = value
    def __str__(self):
       return repr(self.parameter)


def loadOutFile(filename):

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

    iftm = SASM.IFTM(P, R, Perr, Jexp, qshort, Jerr, Jreg, results, Ireg, qfull)

    return [iftm]

def writeGnomCFG(fname, outname, dmax, args):
    #This writes the GNOM CFG file, using the arguments passed into the function.
    datadir = os.path.dirname(fname)

    f = open(os.path.join(datadir, 'gnom.cfg'),'w')

    f.write('This line intentionally left blank\n')
    f.write('PRINTER C [      postscr     ]  Printer type\n')
    if 'form' in args and args['form'] != '':
        f.write('FORFAC  C [         %s         ]  Form factor file (valid for JOB=2)\n' %(args['form']))
    else:
        f.write('FORFAC  C [                  ]  Form factor file (valid for JOB=2)\n')
    if 'expert' in args and args['expert'] != '':
        f.write('EXPERT  C [     %s         ]  File containing expert parameters\n' %(args['expert']))
    else:
        f.write('EXPERT  C [     none         ]  File containing expert parameters\n')
    f.write('INPUT1  C [        %s        ]  Input file name (first file)\n' %(fname))
    f.write('INPUT2  C [       none       ]  Input file name (second file)\n')
    f.write('NSKIP1  I [       0         ]  Number of initial points to skip\n')
    f.write('NSKIP2  I [        0         ]  Number of final points to skip\n')
    f.write('OUTPUT  C [       %s         ]  Output file\n' %(outname))
    if 'angular' in args and args['angular'] != 1:
        f.write('ISCALE  I [       %s         ]  Angular scale of input data\n' %(args['angular']))
    else:
        f.write('ISCALE  I [        1         ]  Angular scale of input data\n')
    f.write('PLOINP  C [       n          ]  Plotting flag: input data (Y/N)\n')
    f.write('PLORES  C [       n          ]  Plotting flag: results    (Y/N)\n')
    f.write('EVAERR  C [       y          ]  Error flags: calculate errors   (Y/N)\n')
    f.write('PLOERR  C [       n          ]  Plotting flag: p(r) with errors (Y/N)\n')
    f.write('PLOALPH C [       n          ]  Plotting flag: alpha distribution (Y/N)\n')
    f.write('LKERN   C [       n          ]  Kernel file status (Y/N)\n')
    if 'system' in args and args['system'] != 0:
        f.write('JOBTYP  I [       %s         ]  Type of system (0/1/2/3/4/5/6)\n' %(args['system']))
    else:
        f.write('JOBTYP  I [       0          ]  Type of system (0/1/2/3/4/5/6)\n')
    if 'rmin' in args and args['rmin'] != -1:
        f.write('RMIN    R [        %s         ]  Rmin for evaluating p(r)\n' %(args['rmin']))
    else:
        f.write('RMIN    R [                 ]  Rmin for evaluating p(r)\n')
    f.write('RMAX    R [        %s        ]  Rmax for evaluating p(r)\n' %(str(dmax)))
    if 'rmin_zero' in args and args['rmin_zero'] != '':
        f.write('LZRMIN  C [      %s          ]  Zero condition at r=RMIN (Y/N)\n' %(args['rmin_zero']))
    else:
        f.write('LZRMIN  C [       Y          ]  Zero condition at r=RMIN (Y/N)\n')
    if 'rmax_zero' in args and args['rmax_zero'] != '':
        f.write('LZRMAX  C [      %s          ]  Zero condition at r=RMAX (Y/N)\n' %(args['rmax_zero']))
    else:
        f.write('LZRMAX  C [       Y          ]  Zero condition at r=RMAX (Y/N)\n')
    f.write('KERNEL  C [       kern.bin   ]  Kernel-storage file\n')
    f.write('DEVIAT  R [      0.0         ]  Default input errors level\n')
    f.write('IDET    I [       0          ]  Experimental set up (0/1/2)\n')
    f.write('FWHM1   R [       0.0        ]  FWHM for 1st run\n')
    f.write('FWHM2   R [                  ]  FWHM for 2nd run\n')
    f.write('AH1     R [                  ]  Slit-height parameter AH (first  run)\n')
    f.write('LH1     R [                  ]  Slit-height parameter LH (first  run)\n')
    f.write('AW1     R [                  ]  Slit-width  parameter AW (first  run)\n')
    f.write('LW1     R [                  ]  Slit-width  parameter LW (first  run)\n')
    f.write('AH2     R [                  ]  Slit-height parameter AH (second run)\n')
    f.write('LH2     R [                  ]  Slit-height parameter LH (second run)\n')
    f.write('AW2     R [                  ]  Slit-width  parameter AW (second run)\n')
    f.write('LW2     R [                  ]  Slit-width  parameter LW (second run)\n')
    f.write('SPOT1   C [                  ]  Beam profile file (first run)\n')
    f.write('SPOT2   C [                  ]  Beam profile file (second run)\n')
    if 'alpha' in args and args['alpha'] !=0.0:
        f.write('ALPHA   R [      %s         ]  Initial ALPHA\n' %(str(args['alpha'])))
    else:
        f.write('ALPHA   R [      0.0         ]  Initial ALPHA\n')
    if 'npts' in args and args['npts'] !=101:
        f.write('NREAL   R [       %s        ]  Number of points in real space\n' %(str(args['npts'])))
    else:
        f.write('NREAL   R [       101        ]  Number of points in real space\n')
    f.write('COEF    R [                  ]\n')
    if 'radius56' in args and args['radius56'] != -1:
        f.write('RAD56   R [         %s         ]  Radius/thickness (valid for JOB=5,6)\n' %(args['radius56']))
    else:
        f.write('RAD56   R [                  ]  Radius/thickness (valid for JOB=5,6)\n')
    f.write('NEXTJOB C [        n         ]\n')


    f.close()




def setATSASEnv(atsasDir):
    my_env = os.environ.copy()
    my_env["PATH"] = my_env["PATH"] + '{}{}'.format(os.pathsep, atsasDir) #Not ideal, what if there's another ATSAS path?
    my_env["ATSAS"] = os.path.split(atsasDir.rstrip(os.sep))[0] #Can only have one thing in ATSAS env variable!

    return my_env

def runGnom(fname, outname, dmax, args, path, atsasDir, new_gnom = False, ):
    #This function runs GNOM from the atsas package. It can do so without writing a GNOM cfg file.
    #It takes as input the filename to run GNOM on, the output name from the GNOM file, the dmax to evaluate
    #at, and a dictionary of arguments, which can be used to set the optional GNOM arguments.
    #Using the GNOM cfg file is significantly faster than catching the interactive input and pass arguments there.

    #Solution for non-blocking reads adapted from stack overflow
    #http://stackoverflow.com/questions/375427/non-blocking-read-on-a-subprocess-pipe-in-python
    def enqueue_output(out, queue):
        line = 'test'
        line2=''
        while line != '':
            line = out.read(1)

            if not isinstance(line, str):
                line = str(line, encoding='UTF-8')

            line2+=line
            if line == ':':
                queue.put_nowait([line2])
                line2=''

        out.close()

    if new_gnom:
        cfg = False
    else:
        cfg = True

    if new_gnom:
        #Check whether everything can be set at the command line:
        fresh_settings = RAWSettings.RawGuiSettings().getAllParams()

        key_ref = { 'gnomExpertFile'    : 'expert',
                    'gnomForceRminZero' : 'rmin_zero',
                    'gnomForceRmaxZero' : 'rmax_zero',
                    'gnomNPoints'       : 'npts',
                    'gnomInitialAlpha'  : 'alpha',
                    'gnomAngularScale'  : 'angular',
                    'gnomSystem'        : 'system',
                    'gnomFormFactor'    : 'form',
                    'gnomRadius56'      : 'radius56',
                    'gnomRmin'          : 'rmin',
                    'gnomFWHM'          : 'fwhm',
                    'gnomAH'            : 'ah',
                    'gnomLH'            : 'lh',
                    'gnomAW'            : 'aw',
                    'gnomLW'            : 'lw',
                    'gnomSpot'          : 'spot',
                    'gnomExpt'          : 'expt'
                    }

        cmd_line_keys = {'rmin_zero', 'rmax_zero', 'system', 'rmin',
            'radiu56', 'npts', 'alpha'}

        changed = []

        for key in fresh_settings:
            if key in key_ref:
                if fresh_settings[key][0] != args[key_ref[key]]:
                    changed.append((key_ref[key]))

        if set(changed) <= cmd_line_keys:
            use_cmd_line = True
        else:
            use_cmd_line = False

    opsys = platform.system()
    if opsys == 'Windows':
        gnomDir = os.path.join(atsasDir, 'gnom.exe')
        shell=False
    else:
        gnomDir = os.path.join(atsasDir, 'gnom')
        shell=True

    if os.path.exists(gnomDir):

        my_env = setATSASEnv(atsasDir)

        if cfg:
            writeGnomCFG(fname, outname, dmax, args)

            proc = subprocess.Popen('"%s"' %(gnomDir), shell=shell,
                stdin=subprocess.PIPE, stdout=subprocess.PIPE,
                stderr=subprocess.STDOUT, cwd=path, env=my_env)
            proc.communicate('\r\n')

        else:
            if os.path.isfile(os.path.join(path, 'gnom.cfg')):
                os.remove(os.path.join(path, 'gnom.cfg'))

            if new_gnom and use_cmd_line:
                cmd = '"%s" --rmax=%s --output="%s"' %(gnomDir, str(dmax), outname)

                if args['npts'] > 0:
                    cmd = cmd + ' --nr=%s'%(str(args['npts']))

                if 'system' in changed:
                    cmd = cmd+' --system=%s' %(str(args['system']))

                if 'rmin' in changed:
                    cmd = cmd+' --rmin=%s' %(str(args['rmin']))

                if 'radius56' in changed:
                    cmd = cmd + ' --rad56=%s' %(str(args['radius56']))

                if 'rmin_zero' in changed:
                    cmd = cmd + ' --force-zero-rmin=%s' %(args['rmin_zero'])

                if 'rmax_zero' in changed:
                    cmd = cmd + ' --force-zero-rmax=%s' %(args['rmax_zero'])

                if 'alpha' in changed:
                    cmd = cmd + ' --alpha=%s' %(str(args['alpha']))

                cmd = cmd + ' "%s"' %(fname)

                proc = subprocess.Popen(cmd, shell=shell, cwd=path, env=my_env)

                proc.wait()

            else:

                gnom_q = queue.Queue()

                proc = subprocess.Popen('"%s"' %(gnomDir), shell=shell,
                    stdin=subprocess.PIPE, stdout=subprocess.PIPE,
                    stderr=subprocess.STDOUT, cwd=path, universal_newlines=True,
                    bufsize=1, env=my_env)
                gnom_t = threading.Thread(target=enqueue_output, args=(proc.stdout, gnom_q))
                gnom_t.daemon = True
                gnom_t.start()

                previous_line = ''
                previous_line2 = ''

                while proc.poll() is None:
                    data = None
                    try:
                        data = gnom_q.get_nowait()
                        data = data[0]
                        gnom_q.task_done()
                    except queue.Empty:
                        pass

                    if data is not None:
                        current_line = data

                        if data.find('[ postscr     ] :') > -1:
                            proc.stdin.write('\r\n') #Printer type, default is postscr

                        elif data.find('Input data') > -1:
                            proc.stdin.write('%s\r\n' %(fname)) #Input data, first file. No default.

                        elif data.find('Output file') > -1:
                            proc.stdin.write('%s\r\n' %(outname)) #Output file, default is gnom.out

                        elif data.find('No of start points to skip') > -1:
                            if 's_skip' in args and args['s_skip'] != '':
                                proc.stdin.write('%s\r\n' %(args['s_skip']))
                            else:
                                proc.stdin.write('\r\n') #Number of start points to skip, default is 0

                        elif data.find('Input data, second file') > -1:
                            proc.stdin.write('\r\n') #Input data, second file, default is none

                        elif data.find('No of end points to omit') > -1:
                            if 'e_skip' in args and args['e_skip'] != '':
                                proc.stdin.write('%s\r\n' %(args['e_skip']))
                            else:
                                proc.stdin.write('\r\n') #Number of end poitns to omit, default is 0

                        elif data.find('Default input errors level') > -1:
                            proc.stdin.write('\r\n') #Default input errors level, default 0.0

                        elif data.find('Angular scale') > -1:
                            if 'angular' in args and args['angular'] != '':
                                proc.stdin.write('%s\r\n' %(args['angular']))
                            else:
                                proc.stdin.write('\r\n') #Angular scale, default 1

                        elif data.find('Plot input data') > -1:
                            proc.stdin.write('n\r\n') #Plot input data, default yes

                        elif data.find('File containing expert parameters') > -1:
                            if 'expert' in args and args['expert'] != '':
                                proc.stdin.write('%s\r\n' %(args['expert']))
                            else:
                                proc.stdin.write('\r\n') #File containing expert parameters, default none

                        elif data.find('Kernel already calculated') > -1:
                            proc.stdin.write('\r\n') #Kernel already calculated, default no

                        elif data.find('Type of system') > -1 or data.find('arbitrary monodisperse)') > -1:
                            if 'system' in args and args['system'] != '':
                                proc.stdin.write('%s\r\n' %(args['system']))
                            else:
                                proc.stdin.write('\r\n') #Type of system, default 0 (P(r) function)

                        elif (data.find('Zero condition at r=rmin') > -1 and data.find('[') > -1) or (previous_line.find('Zero condition at r=rmin') > -1 and previous_line.find('(') > -1):
                            if 'rmin_zero' in args and args['rmin_zero'] != '':
                                proc.stdin.write('%s\r\n' %(args['rmin_zero']))
                            else:
                                proc.stdin.write('\r\n') #Zero condition at r=rmin, default is yes

                        elif (data.find('Zero condition at r=rmax') > -1 and data.find('[') > -1) or (previous_line.find('Zero condition at r=rmax') > -1 and previous_line.find('(') > -1):
                            if 'rmax_zero' in args and args['rmax_zero'] != '':
                                proc.stdin.write('%s\r\n' %(args['rmax_zero']))
                            else:
                                proc.stdin.write('\r\n') #Zero condition at r=rmax, default is yes

                        elif data.find('Rmax for evaluating p(r)') > -1 or data.find('Maximum particle diameter') > -1 or data.find('Maximum characteristic size') > -1 or data.find('Maximum particle thickness') > -1 or data.find('Maximum diameter of particle') > -1 or data.find('Maximum height of cylinder') > -1 or data.find('Maximum outer shell radius') > -1:
                            proc.stdin.write('%s\r\n' %(str(dmax))) #Rmax for evaluating p(r), no default (DMAX!)

                        elif (data.find('Number of points in real space') > -1 and data.find('[') > -1) or previous_line.find('Number of points in real space?') > -1:
                            if 'npts' in args and args['npts'] != -1:
                                proc.stdin.write('%s\r\n' %(str(args['npts'])))
                            else:
                                proc.stdin.write('\r\n') #Number of points in real space, default is 171

                        elif data.find('Kernel-storage file name') > -1:
                            proc.stdin.write('\r\n') #Kernal-storage file name, default is kern.bin

                        elif (data.find('Experimental setup') > -1 and data.find('[') > -1) or data.find('point collimation)') > -1:
                            if 'gnomExp' in args:
                                proc.stdin.write('%s\r\n' %(str(args['gnomExp'])))
                            else:
                                proc.stdin.write('\r\n') #Experimental setup, default is 0 (no smearing)

                        elif data.find('Initial ALPHA') > -1 or previous_line.find('Initial alpha') > -1:
                            if 'alpha' in args and args['alpha'] != 0.0:
                                proc.stdin.write('%s\r\n' %(str(args['alpha'])))
                            else:
                                proc.stdin.write('\r\n') #Initial ALPHA, default is 0.0

                        elif data.find('Plot alpha distribution') > -1:
                            proc.stdin.write('n\r\n') #Plot alpha distribution, default is yes

                        elif data.find('Plot results') > -1:
                            proc.stdin.write('n\r\n') #Plot results, default is no

                        elif data.find('Your choice') > -1:
                            proc.stdin.write('\r\n') #Choice when selecting one of the following options, CR for EXIT

                        elif data.find('Evaluate errors') > -1:
                            proc.stdin.write('\r\n') #Evaluate errors, default yes

                        elif data.find('Plot p(r) with errors') > -1:
                            proc.stdin.write('n\r\n') #Plot p(r) with errors, default yes

                        elif data.find('Next data set') > -1:
                            proc.stdin.write('\r\n') #Next data set, default no

                        elif data.find('Rmin for evaluating p(r)') > -1 or data.find('Minimum characteristic size') > -1 or previous_line.find('Minimum height of cylinder') > -1 or previous_line.find('Minimum outer shell radius') > -1:
                            if 'rmin' in args and args['rmin'] != -1 and args['rmin'] >= 0:
                                proc.stdin.write('%s\r\n' %(str(args['rmin']))) #Rmin, required for some job types
                            else:
                                proc.stdin.write('\r\n' %(str(args['rmin']))) #Default is 0

                        elif data.find('Form factor file for JOB=2') > -1 or data.find('Form Factor file') > -1:
                            proc.stdin.write('%s\r\n' %(str(args['form'])))

                        elif data.find('Cylinder radius') > -1 or data.find('Relative shell thickness') > -1:
                            if 'radius56' in args and args['radius56'] != -1:
                                proc.stdin.write('%s\r\n') %(str(args['radius56'])) #Cylinder radius / relative thickness
                            else:
                                proc.stdin.write('\r\n') #Default is 0

                        elif data.find('FWHM for the first run') > 1:
                            #Need something here
                            if 'fwhm' in args and args['fwhm'] != -1:
                                proc.stdin.write('%s\r\n') %(str(args['fwhm'])) #Beam FWHM
                            else:
                                proc.stdin.write('\r\n') #Default is 0

                        elif data.find('Slit-height parameter AH') > -1 or previous_line.find('Slit height parameter A') > -1:
                            if 'ah' in args and args['ah'] != -1:
                                proc.stdin.write('%s\r\n') %(str(args['ah'])) #Beam height in the detector plane
                            else:
                                proc.stdin.write('\r\n') #Default is 0

                        elif data.find('Slit-height parameter LH') > -1 or previous_line.find('Slit height parameter L') > -1:
                            if 'lh' in args and args['lh'] != -1:
                                proc.stdin.write('%s\r\n') %(str(args['lh'])) #Half the height  difference between top and bottom edge of beam in detector plane
                            else:
                                proc.stdin.write('\r\n') #Default is 0

                        elif data.find('parameter AW') > -1 or previous_line.find('Slit width parameter A') > -1:
                            if 'aw' in args and args['aw'] != -1:
                                proc.stdin.write('%s\r\n') %(str(args['aw'])) #Projection of beam width in detectgor plane
                            else:
                                proc.stdin.write('\r\n') #Default is 0

                        elif data.find('parameter LW') > -1 or previous_line.find('Slit width parameter L') > -1:
                            if 'lw' in args and args['lw'] != -1:
                                proc.stdin.write('%s\r\n') %(str(args['lw'])) #Half of the width difference bewteen top and bottom edge of beam projection in detector plane
                            else:
                                proc.stdin.write('\r\n') #Default is 0

                        elif data.find('Beam profile file') > -1:
                            if 'spot' in args and args['spot'] != -1:
                                proc.stdin.write('%s\r\n') %(str(args['spot'])) #Beam profile file
                            else:
                                proc.stdin.write('\r\n') #Default is none, doesn't use a profile

                        #Prompts from GNOM5
                        elif previous_line.find('(e) expert') > -1:
                            proc.stdin.write('\r\n') #Default is user, good for now. Looks like setting weights is now done in expert mode rather than with a file, so eventually incorporate that.

                        elif previous_line.find('First point to use') > -1:
                            if 's_skip' in args and args['s_skip'] != '':
                                proc.stdin.write('%i\r\n' %(int(args['s_skip'])+1))
                            else:
                                proc.stdin.write('\r\n') #Number of start points to skip, plus one, default is 1

                        elif previous_line.find('Last point to use') > -1:
                            tot_pts = int(current_line.split()[0].strip().rstrip(')'))
                            if 'e_skip' in args and args['e_skip'] != '':
                                proc.stdin.write('%i\r\n' %(tot_pts-int(args['e_skip'])))
                            else:
                                proc.stdin.write('\r\n') #Number of start points to skip, plus one, default is 1

                        #Not implimented yet in RAW.
                        elif previous_line2.find('Slit height setup') > -1:
                            pass

                        elif previous_line2.find('Slight width setup') > -1:
                            pass

                        elif previous_line2.find('Wavelength distribution setup') > -1:
                            pass

                        elif previous_line.find('FWHM of wavelength') > -1:
                            pass

                        elif data.find('Slit height experimental profile file') > -1:
                            pass

                        previous_line2 = previous_line
                        previous_line = current_line

                gnom_t.join()
        try:
            iftm=SASFileIO.loadOutFile(os.path.join(path, outname))[0]
        except IOError:
            raise SASExceptions.GNOMError('No GNOM output file present. GNOM failed to run correctly')

        if cfg:
            try:
                os.remove(os.path.join(path, 'gnom.cfg'))
            except Exception as e:
                print(e)
                print('GNOM cleanup failed to delete gnom.cfg!')

        if not new_gnom:
            try:
                os.remove(os.path.join(path, 'kern.bin'))
            except Exception as e:
                print(e)
                print('GNOM cleanup failed to delete kern.bin!')

        return iftm

    else:
        print('Cannot find ATSAS')
        raise SASExceptions.NoATSASError('Cannot find gnom.')
        return Non



def gnom(profile, dmax, rg=None, idx_min=None, idx_max=None, dmax_zero=True, alpha=0,
    atsas_dir=None, use_rg_from='guinier', use_guinier_start=True,
    cut_8rg=False, write_profile=True, datadir=None, filename=None,
    save_ift=False, savename=None, settings=None, dmin_zero=True, npts=0,
    angular_scale=1, system=0, form_factor='', radius56=-1, rmin=-1, fwhm=-1,
    ah=-1, lh=-1, aw=-1, lw=-1, spot=''):
    """
    Calculates the IFT and resulting P(r) function using gnom from the
    ATSAS package. This requires a separate installation of the ATSAS package
    to use. If gnom fails, values of ``None``, -1, or ``''`` are returned.

    Parameters
    ----------
    profile: :class:`bioxtasraw.SASM.SASM`
        The profile to calculate the IFT for. If using write_file false, you
        can pass None here.
    dmax: float
        The Dmax to be used in calculating the IFT.
    rg: float, optional
        The Rg to be used in calculating the 8/rg cutoff, if cut_8/rg is
        True. If not provided, then the Rg is taken from the analysis
        dictionary of the profile, in conjunction with the use_rg_from setting.
    idx_min: int, optional
        The index of the q vector that corresponds to the minimum q point
        to be used in the IFT. Defaults to the first point in the q vector.
        Overrides use_guinier_start.
    idx_max: int, optional
        The index of the q vector that corresponds to the maximum q point
        to be used in the IFT. Defaults to the last point in the q vector.
        If write_profile is false and no profile is provided then this cannot
        be set, and datgnom will truncate to 8/Rg automatically. Overrides
        cut_8rg.
    dmax_zero: bool, optional
        If True, force P(r) function to zero at Dmax.
    alpha: bool, optional
        If not zero, force alpha value to the input value. If zero (default),
        then alpha is automatically determined by GNOM.
    atsas_dir: str, optional
        The directory of the atsas programs (the bin directory). If not provided,
        the API uses the auto-detected directory.
    use_rg_from: {'guinier', 'gnom', 'bift'} str, optional
        Determines whether the Rg value used for the 8/rg cutoff calculation
        is from the Guinier fit, or the GNOM or BIFT P(r) function. Ignored if
        the rg parameter is provided. Only used if cut_8/rg is True.
    use_guiner_start: bool, optional
        If set to True, and no idx_min is provided, if a Guinier fit has
        been done for the input profile, the start point of the Guinier fit is
        used as the start point for the IFT. Ignored if there is no input profile.
    cut_8rg: bool, optional
        If set to True and no idx_max is provided, then the profile is
        automatically truncated at q=8/Rg.
    write_profile: bool, optional
        If True, the input profile is written to file. If False, then the
        input profile is ignored, and the profile specified by datadir and
        filename is used. This is convenient if you are trying to process
        a lot of files that are already on disk, as it saves having to read
        in each file and then save them again. Defaults to True. If False,
        you must provide a value for the rg parameter.
    datadir: str, optional
        If write_file is False, this is used as the path to the scattering
        profile on disk.
    filename: str, optional.
        If write_file is False, this is used as the filename of the scattering
        profile on disk.
    save_ift: bool, optional
        If True, the IFT from datgnom (.out file) is saved on disk. Requires
        specification of datadir and savename parameters.
    savename: str, optional
        If save_ift is True, this is used as the filename of the .out file on
        disk. This should just be the filename, no path. The datadir parameter
        is used as the parth.
    settings: :class:`bioxtasraw.RAWSettings.RAWSettings`, optional
        RAW settings containing relevant parameters. If provided, the
        dmin_zero, npts, angular_scale, system, form_factor, radius56, rmin,
        fwhm, ah, lh, aw, lw, and spot parameters will be overridden with the
        values in the settings. Default is None.
    dmin_zero: bool, optional
        If True, force P(r) function to zero at Dmin.
    npts: int, optional
        If provided, fixes the number of points in the P(r) function. If 0
        (default), number of points in th P(r) function is automatically
        determined.
    angular_scale: int, optional
        Defines the angular scale of the data as given in the GNOM manual.
        Default is 1/Angstrom.
    system: int, optional
        Defines the job type as in the GNOM manual. Default is 0, a
        monodisperse system.
    form_factor: str, optional
        Path to the form factor file for system type 2. Default is not used.
    radius56: float, optional
        The radius/thickness for system type 5/6. Default is not used.
    rmin: float, optional
        Minimum size for system types 1-6. Default is not used.
    fwhm: float, optional
        Beam FWHM. Default is not used.
    ah: float, optional
        Slit height parameter A as defined in the GNOM manual. Default is not
        used.
    lh: float, optional
        Slit height parameter L as defined in the GNOM manual. Default is not
        used.
    aw: float, optional
        Slit width parameter A as defined in the GNOM manual. Default is not
        used.
    lw: float, optional
        Slit width parameter L as defined in the GNOM manual. Default is not
        used.
    spot: str, optional
        Beam profile file. Default is not used.


    Returns
    -------
    ift: :class:`bioxtasraw.SASM.IFTM`
        The IFT calculated by GNOM from the input profile.
    dmax: float
        The maximum dimension of the P(r) function.
    rg: float
        The real space radius of gyration (Rg) from the P(r) function.
    i0: float
        The real space scattering at zero angle (I(0)) from the P(r) function.
    rg_err: float
        The uncertainty in the real space radius of gyration (Rg) from the P(r)
        function.
    i0_err: float
        The uncertainty in the real space scattering at zero angle (I(0)) from
        the P(r) function.
    total_est: float
        The GNOM total estimate.
    chi_sq: float
        The chi squared value of the fit of the scattering profile calculated
        from the P(r) function to the input scattering profile.
    alpha: float
        The alpha value determined by datgnom.
    quality: str
        The GNOM qualitative interpretation of the total estimate.
    """

    if atsas_dir is None:
        atsas_dir = __default_settings.get('ATSASDir')

    # Set input and output filenames and directory
    if not save_ift and write_profile:
        datadir = os.path.abspath(os.path.expanduser(tempfile.gettempdir()))

        filename = tempfile.NamedTemporaryFile(dir=datadir).name
        filename = os.path.split(filename)[-1] + '.dat'

    elif save_ift and write_profile:
        filename = profile.getParameter('filename')

    datadir = os.path.abspath(os.path.expanduser(datadir))

    if not save_ift:
        savename = tempfile.NamedTemporaryFile(dir=datadir).name
        while os.path.isfile(savename):
            savename = tempfile.NamedTemporaryFile(dir=datadir).name

        savename = os.path.split(savename)[1]
        savename = savename+'.out'



    # Save profile if necessary, truncating q range as appropriate
    if write_profile:
        analysis_dict = profile.getParameter('analysis')

        if rg is None:
            if use_rg_from == 'guinier':
                guinier_dict = analysis_dict['guinier']
                rg = float(guinier_dict['Rg'])

            elif use_rg_from == 'gnom':
                gnom_dict = analysis_dict['GNOM']
                rg = float(gnom_dict['Real_Space_Rg'])

            elif use_rg_from == 'bift':
                bift_dict = analysis_dict['BIFT']
                rg = float(bift_dict['Real_Space_Rg'])

        save_profile = copy.deepcopy(profile)

        if idx_min is None and use_guinier_start:
            analysis_dict = profile.getParameter('analysis')
            if 'guinier' in analysis_dict:
                guinier_dict = analysis_dict['guinier']
                idx_min = int(guinier_dict['nStart']) - profile.getQrange()[0]
            else:
                idx_min = 0

        elif idx_min is None:
            idx_min = 0

        if idx_max is not None:
            save_profile.setQrange((idx_min, idx_max+1))
        else:
            if cut_8rg:
                q = save_profile.getQ()
                idx_max = np.argmin(np.abs(q-(8/rg)))
            else:
                _, idx_max = save_profile.getQrange()

            save_profile.setQrange((idx_min, idx_max))

        SASFileIO.writeRadFile(save_profile, os.path.join(datadir, filename),
            False)

    else:
        if idx_min is None and use_guinier_start and profile is not None:
            analysis_dict = profile.getParameter('analysis')
            if 'guinier' in analysis_dict:
                guinier_dict = analysis_dict['guinier']
                idx_min = int(guinier_dict['nStart']) - profile.getQrange()[0]
            else:
                idx_min = 0

        elif idx_min is None:
            idx_min=0

        if idx_max is None and profile is not None:
            _, idx_max = save_profile.getQrange()


    #Initialize settings
    if settings is not None:
        gnom_settings = {
            'expert'        : settings.get('gnomExpertFile'),
            'rmin_zero'     : settings.get('gnomForceRminZero'),
            'rmax_zero'     : dmax_zero,
            'npts'          : settings.get('gnomNPoints'),
            'alpha'         : alpha,
            'angular'       : settings.get('gnomAngularScale'),
            'system'        : settings.get('gnomSystem'),
            'form'          : settings.get('gnomFormFactor'),
            'radius56'      : settings.get('gnomRadius56'),
            'rmin'          : settings.get('gnomRmin'),
            'fwhm'          : settings.get('gnomFWHM'),
            'ah'            : settings.get('gnomAH'),
            'lh'            : settings.get('gnomLH'),
            'aw'            : settings.get('gnomAW'),
            'lw'            : settings.get('gnomLW'),
            'spot'          : settings.get('gnomSpot'),
            'expt'          : settings.get('gnomExpt')
            }

    else:
        settings = __default_settings

        gnom_settings = {
            'expert'        : settings.get('gnomExpertFile'),
            'rmin_zero'     : dmin_zero,
            'rmax_zero'     : dmax_zero,
            'npts'          : npts,
            'alpha'         : alpha,
            'angular'       : angular_scale,
            'system'        : system,
            'form'          : form_factor,
            'radius56'      : radius56,
            'rmin'          : rmin,
            'fwhm'          : fwhm,
            'ah'            : ah,
            'lh'            : lh,
            'aw'            : aw,
            'lw'            : lw,
            'spot'          : spot,
            'expt'          : settings.get('gnomExpt')
            }


    # Run the IFT
    ift = SASCalc.runGnom(filename, savename, dmax, gnom_settings, datadir,
        atsas_dir, True)

    # Clean up
    if write_profile and os.path.isfile(os.path.join(datadir, filename)):
        try:
            os.remove(os.path.join(datadir, filename))
        except Exception:
            pass

    if not save_ift and os.path.isfile(os.path.join(datadir, savename)):
        try:
            os.remove(os.path.join(datadir, savename))
        except Exception:
            pass

        if write_profile:
            ift_name = profile.getParameter('filename')
        else:
            ift_name = filename

        ift_name = os.path.splitext(ift_name)[0] + '.out'
        ift.setParameter('filename', ift_name)

    # Save results
    if ift is not None:
        dmax = float(ift.getParameter('dmax'))
        rg = float(ift.getParameter('rg'))
        rg_err = float(ift.getParameter('rger'))
        i0 = float(ift.getParameter('i0'))
        i0_err = float(ift.getParameter('i0er'))
        chi_sq = float(ift.getParameter('chisq'))
        alpha = float(ift.getParameter('alpha'))
        total_est = float(ift.getParameter('TE'))
        quality = ift.getParameter('quality')

        if profile is not None:
            results_dict = {}
            results_dict['Dmax'] = str(dmax)
            results_dict['Total_Estimate'] = str(total_est)
            results_dict['Real_Space_Rg'] = str(rg)
            results_dict['Real_Space_Rg_Err'] = str(rg_err)
            results_dict['Real_Space_I0'] = str(i0)
            results_dict['Real_Space_I0_Err'] = str(i0_err)
            results_dict['GNOM_ChiSquared'] = str(chi_sq)
            results_dict['Alpha'] = str(alpha)
            results_dict['qStart'] = save_profile.getQ()[0]
            results_dict['qEnd'] = save_profile.getQ()[0]
            results_dict['GNOM_Quality_Assessment'] = quality
            analysis_dict['GNOM'] = results_dict
            profile.setParameter('analysis', analysis_dict)

    else:
        dmax = -1
        rg = -1
        rg_err = -1
        i0 = -1
        i0_err = -1
        chi_sq = -1
        alpha = -1
        total_est = -1
        quality = ''

    return ift, dmax, rg, i0, rg_err, i0_err, total_est, chi_sq, alpha, quality








#### AutoRG workup ####


@jit(nopython=True, cache=True, parallel=False)
def linear_func(x, a, b):
    return a+b*x


@jit(nopython=True, cache=True, parallel=False)
def calcRg(q, i, err, transform=True, error_weight=True):
    if transform:
        #Start out by transforming as usual.
        x = np.square(q)
        y = np.log(i)
        yerr = np.absolute(err/i) #I know it looks odd, but it's correct for a natural log
    else:
        x = q
        y = i
        yerr = err

    if error_weight:
        a, b, cov_a, cov_b = weighted_lin_reg(x, y, yerr)
    else:
        a, b, cov_a, cov_b = lin_reg(x, y)

    if b < 0:
        RG=np.sqrt(-3.*b)
        I0=np.exp(a)

        #error in rg and i0 is calculated by noting that q(x)+/-Dq has Dq=abs(dq/dx)Dx, where q(x) is your function you're using
        #on the quantity x+/-Dx, with Dq and Dx as the uncertainties and dq/dx the derviative of q with respect to x.
        RGer=np.absolute(0.5*(np.sqrt(-3./b)))*np.sqrt(np.absolute(cov_b))
        I0er=I0*np.sqrt(np.absolute(cov_a))

    else:
        RG = -1
        I0 = -1
        RGer = -1
        I0er = -1

    return RG, I0, RGer, I0er, a, b



def autoRg(sasm, single_fit=False, error_weight=True):
    #This function automatically calculates the radius of gyration and scattering intensity at zero angle
    #from a given scattering profile. It roughly follows the method used by the autorg function in the atsas package

    q = sasm.q
    i = sasm.i
    err = sasm.err
    qmin, qmax = sasm.getQrange()

    # RM!
    q=q
    i=i
    err=err
    qmin, qmax = (0, len(q))

    q = q[qmin:qmax]
    i = i[qmin:qmax]
    err = err[qmin:qmax]

    try:
        rg, rger, i0, i0er, idx_min, idx_max = autoRg_inner(q, i, err, qmin, single_fit, error_weight)
    except Exception: #Catches unexpected numba errors, I hope
        traceback.print_exc()
        rg = -1
        rger = -1
        i0 = -1
        i0er = -1
        idx_min = -1
        idx_max = -1

    return rg, rger, i0, i0er, idx_min, idx_max

@jit(nopython=True, cache=True, parallel=False)
def autoRg_inner(q, i, err, qmin, single_fit, error_weight):
    #Pick the start of the RG fitting range. Note that in autorg, this is done
    #by looking for strong deviations at low q from aggregation or structure factor
    #or instrumental scattering, and ignoring those. This function isn't that advanced
    #so we start at 0.

    # Note, in order to speed this up using numba, I had to do some unpythonic things
    # with declaring lists ahead of time, and making sure lists didn't have multiple
    # object types in them. It makes the code a bit more messy than the original
    # version, but numba provides a significant speedup.
    data_start = 0

    #Following the atsas package, the end point of our search space is the q value
    #where the intensity has droped by an order of magnitude from the initial value.
    data_end = np.abs(i-i[data_start]/10.).argmin()

    #This makes sure we're not getting some weird fluke at the end of the scattering profile.
    if data_end > len(q)/2.:
        found = False
        idx = 0
        while not found:
            idx = idx +1
            if i[idx]<i[0]/10.:
                found = True
            elif idx == len(q) -1:
                found = True
        data_end = idx

    #Start out by transforming as usual.
    qs = np.square(q)
    il = np.log(i)
    iler = np.absolute(err/i)

    #Pick a minimum fitting window size. 10 is consistent with atsas autorg.
    min_window = 10

    max_window = data_end-data_start

    if max_window<min_window:
        max_window = min_window

    #It is very time consuming to search every possible window size and every possible starting point.
    #Here we define a subset to search.
    tot_points = max_window
    window_step = tot_points//10
    data_step = tot_points//50

    if window_step == 0:
        window_step =1
    if data_step ==0:
        data_step =1

    window_list = [0 for k in range(int(math.ceil((max_window-min_window)/float(window_step)))+1)]

    for k in range(int(math.ceil((max_window-min_window)/float(window_step)))):
        window_list[k] = min_window+k*window_step

    window_list[-1] = max_window

    num_fits = 0

    for w in window_list:
        num_fits = num_fits + int(math.ceil((data_end-w-data_start)/float(data_step)))

    if num_fits < 0:
        num_fits = 1

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
    #This function takes every window size in the window list, stepts it through the data range, and
    #fits it to get the RG and I0. If basic conditions are met, qmin*RG<1 and qmax*RG<1.35, and RG>0.1,
    #We keep the fit.
    for w in window_list:
        for start in range(data_start,data_end-w, data_step):
            x = qs[start:start+w]
            y = il[start:start+w]
            yerr = iler[start:start+w]

            #Remove NaN and Inf values:
            x = x[np.where(np.isfinite(y))]
            yerr = yerr[np.where(np.isfinite(y))]
            y = y[np.where(np.isfinite(y))]


            RG, I0, RGer, I0er, a, b = calcRg(x, y, yerr, transform=False, error_weight=error_weight)

            if RG>0.1 and q[start]*RG<1 and q[start+w-1]*RG<1.35 and RGer/RG <= 1:

                r_sqr = 1 - np.square(il[start:start+w]-linear_func(qs[start:start+w], a, b)).sum()/np.square(il[start:start+w]-il[start:start+w].mean()).sum()

                if r_sqr > .15:
                    chi_sqr = np.square((il[start:start+w]-linear_func(qs[start:start+w], a, b))/iler[start:start+w]).sum()

                    #All of my reduced chi_squared values are too small, so I suspect something isn't right with that.
                    #Values less than one tend to indicate either a wrong degree of freedom, or a serious overestimate
                    #of the error bars for the system.
                    dof = w - 2.
                    reduced_chi_sqr = chi_sqr/dof

                    start_list[current_fit] = start
                    w_list[current_fit] = w
                    q_start_list[current_fit] = q[start]
                    q_end_list[current_fit] = q[start+w-1]
                    rg_list[current_fit] = RG
                    rger_list[current_fit] = RGer
                    i0_list[current_fit] = I0
                    i0er_list[current_fit] = I0er
                    qrg_start_list[current_fit] = q[start]*RG
                    qrg_end_list[current_fit] = q[start+w-1]*RG
                    rsqr_list[current_fit] = r_sqr
                    chi_sqr_list[current_fit] = chi_sqr
                    reduced_chi_sqr_list[current_fit] = reduced_chi_sqr

                    # fit_list[current_fit] = [start, w, q[start], q[start+w-1], RG, RGer, I0, I0er, q[start]*RG, q[start+w-1]*RG, r_sqr, chi_sqr, reduced_chi_sqr]
                    success[current_fit] = 1

            current_fit = current_fit + 1
    #Extreme cases: may need to relax the parameters.
    if np.sum(success) == 0:
        #Stuff goes here
        pass

    if np.sum(success) > 0:

        fit_array = np.array([[start_list[k], w_list[k], q_start_list[k],
            q_end_list[k], rg_list[k], rger_list[k], i0_list[k], i0er_list[k],
            qrg_start_list[k], qrg_end_list[k], rsqr_list[k], chi_sqr_list[k],
            reduced_chi_sqr_list[k]] for k in range(num_fits) if success[k]==1])

        #Now we evaluate the quality of the fits based both on fitting data and on other criteria.

        #Choice of weights is pretty arbitrary. This set seems to yield results similar to the atsas autorg
        #for the few things I've tested.
        qmaxrg_weight = 1
        qminrg_weight = 1
        rg_frac_err_weight = 1
        i0_frac_err_weight = 1
        r_sqr_weight = 4
        reduced_chi_sqr_weight = 0
        window_size_weight = 4

        weights = np.array([qmaxrg_weight, qminrg_weight, rg_frac_err_weight, i0_frac_err_weight, r_sqr_weight,
                            reduced_chi_sqr_weight, window_size_weight])

        quality = np.zeros(len(fit_array))

        max_window_real = float(window_list[-1])

        # all_scores = [np.array([]) for k in range(len(fit_array))]

        #This iterates through all the fits, and calculates a score. The score is out of 1, 1 being the best, 0 being the worst.
        indices =list(range(len(fit_array)))
        for a in indices:
            k=int(a) #This is stupid and should not be necessary. Numba bug?

            #Scores all should be 1 based. Reduced chi_square score is not, hence it not being weighted.
            qmaxrg_score = 1-abs((fit_array[k,9]-1.3)/1.3)
            qminrg_score = 1-fit_array[k,8]
            rg_frac_err_score = 1-fit_array[k,5]/fit_array[k,4]
            i0_frac_err_score = 1 - fit_array[k,7]/fit_array[k,6]
            r_sqr_score = fit_array[k,10]
            reduced_chi_sqr_score = 1/fit_array[k,12] #Not right
            window_size_score = fit_array[k,1]/max_window_real

            scores = np.array([qmaxrg_score, qminrg_score, rg_frac_err_score, i0_frac_err_score, r_sqr_score,
                               reduced_chi_sqr_score, window_size_score])

            total_score = (weights*scores).sum()/weights.sum()

            quality[k] = total_score

            # all_scores[k] = scores


        #I have picked an aribtrary threshold here. Not sure if 0.6 is a good quality cutoff or not.
        if quality.max() > 0.6:
            if not single_fit:
                idx = quality.argmax()
                rg = fit_array[:,4][quality>quality[idx]-.1].mean()
                rger = fit_array[:,5][quality>quality[idx]-.1].std()
                i0 = fit_array[:,6][quality>quality[idx]-.1].mean()
                i0er = fit_array[:,7][quality>quality[idx]-.1].std()
                idx_min = int(fit_array[idx,0])
                idx_max = int(fit_array[idx,0]+fit_array[idx,1]-1)
            else:
                idx = quality.argmax()
                rg = fit_array[idx,4]
                rger = fit_array[idx,5]
                i0 = fit_array[idx,6]
                i0er = fit_array[idx,7]
                idx_min = int(fit_array[idx,0])
                idx_max = int(fit_array[idx,0]+fit_array[idx,1]-1)

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

    #We could add another function here, if not good quality fits are found, either reiterate through the
    #the data and refit with looser criteria, or accept lower scores, possibly with larger error bars.

    #returns Rg, Rg error, I0, I0 error, the index of the first q point of the fit and the index of the last q point of the fit
    return rg, rger, i0, i0er, idx_min, idx_max



