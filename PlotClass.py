#!/Users/robmiller/opt/anaconda3/bin/python3
'''
# !/usr/bin/env python3
# -*- coding: utf-8 -*-
'''
"""
Created on Thu Jun  4 14:47:22 2020

@author: Rob
"""

from matplotlib import pyplot as plt
from mpl_toolkits import mplot3d
from itertools import cycle
import matplotlib as mpl
# from Basic_SAXS_Calcs import *
# from SAXS_Calcs import *
import numpy as np
import itertools
from matplotlib.legend_handler import HandlerLine2D, HandlerTuple


class PlotClass:
    
    def __init__(self,notify=True):
        '''
        Describe the class here.
        The initialization of the class sets a few global parameters
        (1) Axes width
        (2) Line width
        (3) Set font to Avenir. Note, if when running locally you receive a traceback error
            along the lines of "falling back to Arial" you need to enable the base Avenir font
            provided by MacOS. If this is important to you, let me know, it's a pretty quick fix!
            Rob: rcm347@cornell.edu
        '''
        self.notify=notify
        if self.notify==True:
            print('--------------------------------------------------------------')
            print('Plot class was called')
            print('--------------------------------------------------------------')
        
        self.axes=plt.rc('axes',linewidth=2)
        self.lines=plt.rc('lines',markeredgewidth=2)
        self.font=plt.rc('font',**{'sans-serif': ['Avenir']})
        # self.calcs=SAXSCalcs(notify=False)

        
    def basicPlot(self,X,Y,plotlabel='',savelabel='',xlabel='',ylabel='NOT PROVIDED'):
        '''
        As simple as it gets
        '''
        fig=plt.figure(figsize=(10,8)) # set figure dimensions
        ax1=fig.add_subplot(1,1,1) # allows us to build more complex plots
        for tick in ax1.xaxis.get_major_ticks():
            tick.label1.set_fontsize(20) # scale for publication needs
            tick.label1.set_fontname('Avenir')
        for tick in ax1.yaxis.get_major_ticks():
            tick.label1.set_fontsize(0.2) # scale for publication needs
            tick.label1.set_fontname('Avenir')
        plt.ylabel(ylabel,size=22)
        plt.xlabel(xlabel,size=22)
        plt.plot(X,Y,'-',label=plotlabel,
                 color='#E55334') # best to use HTML color codes: https://htmlcolorcodes.com
        plt.legend(numpoints=1,fontsize=18,loc='best')
        fig.tight_layout()
        plt.savefig(savelabel+'.png',format='png',
                    bbox_inches='tight',dpi=300)
        plt.show()

    def basicPointPlot(self,X,Y,plotlabel='',savelabel='',xlabel='',ylabel='NOT PROVIDED'):
        '''
        As simple as it gets
        '''
        fig=plt.figure(figsize=(10,8)) # set figure dimensions
        ax1=fig.add_subplot(1,1,1) # allows us to build more complex plots
        for tick in ax1.xaxis.get_major_ticks():
            tick.label1.set_fontsize(20) # scale for publication needs
            tick.label1.set_fontname('Avenir')
        for tick in ax1.yaxis.get_major_ticks():
            tick.label1.set_fontsize(20) # scale for publication needs
            tick.label1.set_fontname('Avenir')
        plt.ylabel(ylabel,size=22)
        plt.xlabel(xlabel,size=22)
        plt.plot(X,Y,marker='o',label=plotlabel,
                 color='#E55334') # best to use HTML color codes: https://htmlcolorcodes.com
        plt.legend(numpoints=1,fontsize=18,loc='best')
        fig.tight_layout()
        plt.savefig(savelabel+'.png',format='png',
                    bbox_inches='tight',dpi=300)
        plt.show()
        
    def singlePlot(self,X,Y,plotlabel='NOT PROVIDED',savelabel='NOT PROVIDED',xlabel='NOT PROVIDED',
                   ylabel='NOT PROVIDED',color='k',linewidth=2,linestyle='-',LogLin=True,setylim=False,
                   ymin=10**(-4),ymax=1,darkmode=False):
        '''
        As simple as it gets
        '''
        if darkmode==True:
            plt.style.use('dark_background')
        else:
            mpl.rcParams.update(mpl.rcParamsDefault)

        fig=plt.figure(figsize=(10,8)) # set figure dimensions
        ax1=fig.add_subplot(1,1,1) # allows us to build more complex plots
        for tick in ax1.xaxis.get_major_ticks():
            tick.label1.set_fontsize(20) # scale for publication needs
            tick.label1.set_fontname('Avenir')
        for tick in ax1.yaxis.get_major_ticks():
            tick.label1.set_fontsize(20) # scale for publication needs
            tick.label1.set_fontname('Avenir')
        plt.ylabel(ylabel,size=22)
        plt.xlabel(xlabel,size=22)
        if LogLin==True:
            plt.semilogy(X,Y,linestyle=linestyle,
                 linewidth=linewidth,
                 label=plotlabel,
                 color=color)
        else:
            plt.plot(X,Y,linestyle=linestyle,
                 linewidth=linewidth,
                 label=plotlabel,
                 color=color) # best to use HTML color codes: https://htmlcolorcodes.com
        if setylim==True:
            ax1.set_ylim(ymin,ymax)
        plt.legend(numpoints=1,fontsize=18,loc='best')
        fig.tight_layout()
        plt.savefig(savelabel+'.png',format='png',
                    bbox_inches='tight',dpi=300)
        plt.show()
        
    def semilogyPlot(self,X,Y,plotlabel='',savelabel='',xlabel='',ylabel='',linewidth=4,
                     darkmode=False):
        '''
        As simple as it gets
        '''
        if darkmode==True:
            plt.style.use('dark_background')
        else:
            mpl.rcParams.update(mpl.rcParamsDefault)

        fig=plt.figure(figsize=(10,8)) # set figure dimensions
        ax1=fig.add_subplot(1,1,1) # allows us to build more complex plots
        for tick in ax1.xaxis.get_major_ticks():
            tick.label1.set_fontsize(20) # scale for publication needs
            tick.label1.set_fontname('Avenir')
        for tick in ax1.yaxis.get_major_ticks():
            tick.label1.set_fontsize(20) # scale for publication needs
            tick.label1.set_fontname('Avenir')
        plt.ylabel(ylabel,size=22)
        plt.xlabel(xlabel,size=22)
        plt.semilogy(X,Y,'-',label=plotlabel,
                     linewidth=linewidth,
                     color='#E55334')
        plt.legend(numpoints=1,fontsize=18,loc='best')
        fig.tight_layout()
        plt.savefig(savelabel+'.png',format='png',
                    bbox_inches='tight',dpi=300)
        plt.show()
        
        
    def twoPlot(self,X,Y1,Y2,plotlabel1='',plotlabel2='',savelabel='',xlabel='',ylabel='',linewidth=2,
                darkmode=False,plot=True):
        '''
        As simple as it gets P2
        '''
        if darkmode==True:
            plt.style.use('dark_background')
            c2='#22A2AA'
        else:
            mpl.rcParams.update(mpl.rcParamsDefault)
            c2='#616061'

        fig=plt.figure(figsize=(8,4)) # set figure dimensions
        ax1=fig.add_subplot(1,1,1) # allows us to build more complex plots
        for tick in ax1.xaxis.get_major_ticks():
            tick.label1.set_fontsize(20) # scale for publication needs
            tick.label1.set_fontname('Avenir')
        for tick in ax1.yaxis.get_major_ticks():
            tick.label1.set_fontsize(20) # scale for publication needs
            tick.label1.set_fontname('Avenir')
        plt.ylabel(ylabel,size=22)
        plt.xlabel(xlabel,size=22)
        plt.plot(X,Y1,'-',label=plotlabel1,
                 linewidth=linewidth, linestyle='dashed',
                 color='#C922BC')
        plt.plot(X,Y2,linestyle=(0, (1, 5)),label=plotlabel2,linewidth=linewidth,
                 color=c2)
        plt.legend(numpoints=1,fontsize=18,loc='best')
        fig.tight_layout()
        plt.savefig(savelabel+'.png',format='png',
                    bbox_inches='tight',dpi=300)
        if plot==True:
            plt.show()

    def twoPlot_variX(self,X1,Y1,X2,Y2,plotlabel1='',plotlabel2='',savelabel='',xlabel='',ylabel='',linewidth=2,
                      LogLin=False,set_ylim=False,ylow=0.0001,yhigh=1,darkmode=False,plot=True):
        '''
        As simple as it gets P2
        '''
        if darkmode==True:
            plt.style.use('dark_background')
        else:
            mpl.rcParams.update(mpl.rcParamsDefault)

        if LogLin==False:
            fig=plt.figure(figsize=(8,4)) # set figure dimensions
            ax1=fig.add_subplot(1,1,1) # allows us to build more complex plots
            for tick in ax1.xaxis.get_major_ticks():
                tick.label1.set_fontsize(20) # scale for publication needs
                tick.label1.set_fontname('Avenir')
            for tick in ax1.yaxis.get_major_ticks():
                tick.label1.set_fontsize(20) # scale for publication needs
                tick.label1.set_fontname('Avenir')
            plt.ylabel(ylabel,size=22)
            plt.xlabel(xlabel,size=22)
            plt.plot(X1,Y1,'-',label=plotlabel1,
                    linewidth=linewidth, linestyle='dashed',
                    color='#C922BC')
            plt.plot(X2,Y2,marker='x',markersize=linewidth,linestyle='None',label=plotlabel2,
                    color='#616061')
            plt.legend(numpoints=1,fontsize=18,loc='best')
            if set_ylim == True:
                ax1.set_ylim(ylow,yhigh)
            # else:
            #     print('Using default y-limit range')
            fig.tight_layout()
            plt.savefig(savelabel+'.png',format='png',
                    bbox_inches='tight',dpi=300)
            plt.show()
        else:
            fig=plt.figure(figsize=(8,4)) # set figure dimensions
            ax1=fig.add_subplot(1,1,1) # allows us to build more complex plots
            for tick in ax1.xaxis.get_major_ticks():
                tick.label1.set_fontsize(20) # scale for publication needs
                tick.label1.set_fontname('Avenir')
            for tick in ax1.yaxis.get_major_ticks():
                tick.label1.set_fontsize(20) # scale for publication needs
                tick.label1.set_fontname('Avenir')
            plt.ylabel(ylabel,size=22)
            plt.xlabel(xlabel,size=22)
            plt.semilogy(X1,Y1,'-',label=plotlabel1,
                    linewidth=linewidth, linestyle='dashed',
                    color='#C922BC')
            plt.semilogy(X2,Y2,marker='x',markersize=linewidth,linestyle='None',label=plotlabel2,
                    color='#616061')
            plt.legend(numpoints=1,fontsize=18,loc='best')
            if set_ylim == True:
                ax1.set_ylim(ylow,yhigh)
            # else:
            #     print('Using default y-limit range')
            fig.tight_layout()
            plt.savefig(savelabel+'.png',format='png',
                    bbox_inches='tight',dpi=300)
            if plot==True:
                plt.show()


    # def kratkyPlot(self,ref_q,ref_i,plotList=[],plotListLabels=[],xlabel='NoLabel',ylabel='NoLabel',scaled=False,
    #                    plotlabel='NoLabelProvided',linewidth=3,savelabel='NoSaveLabelProvided_Kratky',
    #                truncation_q=0.6,darkmode=False):
    #     '''
    #     Generates a basic kratky with options to scale or not scale to the max peak.
    #     input:
            
    #     ** need to deal with empty label list in a clever way
            
    #     output: goes to a max of q=0.6
    #     '''

    #     if scaled==False:
    #         initial_output='Unscaled'
    #     else:
    #         initial_output='Scaled'
    #     print('-------> Kratky plot - %s'%initial_output)

    #     if darkmode==True:
    #         plt.style.use('dark_background')
    #     else:
    #         mpl.rcParams.update(mpl.rcParamsDefault)

    #     fig=plt.figure(figsize=(10,8)) # set figure dimensions
    #     ax1=fig.add_subplot(1,1,1) # allows us to build more complex plots
    #     for tick in ax1.xaxis.get_major_ticks():
    #         tick.label1.set_fontsize(20) # scale for publication needs
    #         tick.label1.set_fontname('Avenir')
    #     for tick in ax1.yaxis.get_major_ticks():
    #         tick.label1.set_fontsize(20) # scale for publication needs
    #         tick.label1.set_fontname('Avenir')
    #     plt.ylabel(ylabel,size=22)
    #     plt.xlabel(xlabel,size=22)

    #     cycol = cycle(['#0E403E','#194E05','#57B036','#AA8B0B','#0F8985']) # set of colors to iterate through.

    #     # set truncation point at high-q for the plot
    #     try:
    #         lowclip = np.where(ref_q >= truncation_q)[0][0]
    #         print('Data was truncated at q=%.2f inverse angstroms'%truncation_q)
    #     except Exception:
    #         lowclip=len(ref_q)-1
    #         print('Data did not go out to q=0.6 inverse angstroms \nso the entire data frame was plotted')

    #     # clip at q=0.2 for determination of scaling constant
    #     scaleclip=np.where(ref_q >= 0.3)[0][0]
    #     # scaleclip=np.where(ref_q == np.max(ref_q))[0][0]
    #     if scaled==False:
    #         if plotList==[]:
    #             plt.plot(ref_q[:lowclip],ref_i[:lowclip]*ref_q[:lowclip]**2,'-',label=plotlabel,
    #                       linewidth=linewidth,color='#0F8985',linestyle='-')
    #             plt.legend(numpoints=1,fontsize=18,loc='best')
    #             fig.tight_layout()
    #             plt.savefig(savelabel + '.png',format='png',
    #                             bbox_inches='tight',dpi=500)
    #             plt.show()
    #             print('No list was provided to scale to. \nTherefore, only the reference profile was plotted.')
    #             # save/etc
    #         else:
    #             plt.plot(ref_q[:lowclip],ref_i[:lowclip]*ref_q[:lowclip]**2,'-',label=plotlabel,
    #                       linewidth=linewidth,color='#0F8985',linestyle='-')
    #             for i,j in zip(plotList,plotListLabels):
    #                 plt.plot(i[0][:lowclip],i[1][:lowclip]*i[0][:lowclip]**2,'-',label=j,
    #                       linewidth=linewidth,color=next(cycol),
    #                          linestyle='-')
    #             plt.legend(numpoints=1,fontsize=18,loc='best')
    #             fig.tight_layout()
    #             plt.savefig(savelabel + '.png',format='png',
    #                             bbox_inches='tight',dpi=500)
    #             plt.show()

    #             # save/etc

                
    #         # need to add in offset capability..
    #     else:
    #         if plotList==[]:
    #             plt.plot(ref_q[:lowclip],ref_i[:lowclip]*ref_q[:lowclip]**2,'-',label=plotlabel,
    #                       linewidth=linewidth,color='#E55334')
    #             plt.legend(numpoints=1,fontsize=18,loc='best')
    #             fig.tight_layout()
    #             plt.savefig(savelabel + '.png',format='png',
    #                             bbox_inches='tight',dpi=500)
    #             plt.show()
    #             print('No list was provided to scale to. \nTherefore, only the reference profile was plotted.')
    #             # save/etc
    #         else:
    #             plt.plot(ref_q[:lowclip],ref_i[:lowclip]*ref_q[:lowclip]**2,'-',label=plotlabel,linewidth=linewidth,color='#0F8985',linestyle='dashed')
    #             for i,j in zip(plotList,plotListLabels):
    #                 scale,offset=self.calcs.superimpose(ref_q[:scaleclip],ref_i[:scaleclip]*ref_q[:scaleclip]**2,0,
    #                                                     len(ref_q[:scaleclip]),[[i[0][:scaleclip],i[1][:scaleclip]*i[0][:scaleclip]**2]],choice='Scale')
    #                 # scale,offset=self.calcs.superimpose(ref_q[:scaleclip],ref_i[:scaleclip]*ref_q[:scaleclip]**2,0,
    #                 #                 0.3,[[i[0][:scaleclip],i[1][:scaleclip]*i[0][:scaleclip]**2]],choice='Scale')
    #                 i[1]=i[1]*scale # applying scaling factor and now we can plot
    #                 print('SCALE:',scale)
    #                 plt.plot(i[0][:lowclip],i[1][:lowclip]*i[0][:lowclip]**2,'-',label=j,
    #                          linewidth=linewidth,color=next(cycol),
    #                          linestyle='dashed')
    #             plt.legend(numpoints=1,fontsize=18,loc='best')
    #             fig.tight_layout()
    #             plt.savefig(savelabel + '.png',format='png',
    #                         bbox_inches='tight',dpi=500)
    #             plt.show()
    #             print('The scaled profiles were:')
    #             for j in plotListLabels:
    #                 print(j)
    #             print('& the reference profile was %s'%plotlabel)
            # apply scaling function..

    def dimless_kratkyPlot(self,plotList=[],plotListLabels=[],rgList=[],I0List=[],xlabel='qR$_{g}$',ylabel='(qR$_{g}$)$^{2}$ $\cdot$ $\\frac{I(q)}{I(0)}$',
                   linewidth=3,savelabel='NoSaveLabelProvided_dimlessKratky',
                   truncation_q=0.6,darkmode=False):
        '''
        Generates a basic kratky with options to scale or not scale to the max peak.
        input:

        ** need to deal with empty label list in a clever way

        output: goes to a max of q=0.6
        '''

        if darkmode == True:
            plt.style.use('dark_background')
        else:
            mpl.rcParams.update(mpl.rcParamsDefault)

        fig = plt.figure(figsize=(10,8))  # set figure dimensions
        ax1 = fig.add_subplot(1,1,1)  # allows us to build more complex plots
        ax1.axvline(x=1.73,linestyle='dashed',color='k')
        ax1.axhline(1.1,linestyle='dashed',color='k')
        ax1.axhline(0.0,linestyle='-',linewidth=3,color='k')
        for tick in ax1.xaxis.get_major_ticks():
            tick.label1.set_fontsize(20)  # scale for publication needs
            tick.label1.set_fontname('Avenir')
        for tick in ax1.yaxis.get_major_ticks():
            tick.label1.set_fontsize(20)  # scale for publication needs
            tick.label1.set_fontname('Avenir')
        plt.ylabel(ylabel,size=22)
        plt.xlabel(xlabel,size=22)

        cycol = cycle(['#0F8985','#0E403E','#194E05','#57B036','#AA8B0B'])  # set of colors to iterate through.

        # set truncation point at high-q for the plot
        try:
            n=0
            for j in plotList:
                lowclip = np.where(j[0] >= truncation_q)[0][0]
                print('Data was truncated at q=%.2f inverse angstroms' % truncation_q)
                n+=1
                if n==1:
                    break
        except Exception:
            for j in plotList:
                lowclip = len(j[0]) - 1
                print(lowclip)
                print('Data did not go out to q=0.6 inverse angstroms \nso the entire data frame was plotted')
                n+=1
                if n==1:
                    break

        for i,j,k,z in zip(plotList,plotListLabels,rgList,I0List):
            plt.plot((i[0][:lowclip]*k),((i[0][:lowclip]*k)**2)*(i[1][:lowclip]/z),'-',label=j,
                     linewidth=linewidth,color=next(cycol),
                     linestyle='dashed')
        plt.legend(numpoints=1,fontsize=18,loc='best')
        fig.tight_layout()
        plt.savefig(savelabel + '.png',format='png',
                    bbox_inches='tight',dpi=500)
        plt.show()

                # save/etc



    def nPlot(self,pairList,labelList,savelabel,xlabel='No Label Provided',ylabel='No Label Provided',
              LogLin=True,LinLin=False,LogLog=False,linewidth=3,
              set_ylim=False,ylow=0.0001,yhigh=1,darkmode=False):
        '''
        :param pairList: list of lists (tuple), must be [[x1,y1],...[xn,yn]]
        :param labelList: list of length n, labeling the sets of tuples in pairList
        :param savelabel:
        :param xlabel:
        :param ylabel:
        :param linewidth:
        :return:
        '''

        if darkmode==True:
            plt.style.use('dark_background')
        else:
            mpl.rcParams.update(mpl.rcParamsDefault)

        fig=plt.figure(figsize=(10,8)) # set figure dimensions
        ax1=fig.add_subplot(1,1,1) # allows us to build more complex plots
        for tick in ax1.xaxis.get_major_ticks():
            tick.label1.set_fontsize(20) # scale for publication needs
            tick.label1.set_fontname('Avenir')
        for tick in ax1.yaxis.get_major_ticks():
            tick.label1.set_fontsize(20) # scale for publication needs
            tick.label1.set_fontname('Avenir')

        if LogLin == True and LinLin == True and LogLog == True: # kicks you out of the function if you set more then one mode to true
            print("Cannot set more than one mode equal to True")
            return
        elif LogLin == True and LinLin == True:
            print("Cannot set more than one mode equal to True")
            return
        elif LogLin == True and LogLog == True:
            print("Cannot set more than one mode equal to True")
            return
        elif LinLin == True and LogLog == True:
            print("Cannot set more than one mode equal to True")
            return

        n=0
        if LogLin==True:
            for i in pairList:
                plt.semilogy(i[0],i[1],
                            label=labelList[n],
                            linewidth=linewidth)
                n+=1
        elif LinLin==True:
            for i in pairList:
                plt.plot(i[0],i[1],
                            label=labelList[n],
                            linewidth=linewidth)
                n+=1
        elif LogLog==True:
            for i in pairList:
                plt.plot(i[0],i[1],
                            label=labelList[n],
                            linewidth=linewidth)
                n+=1
                ax1.set_yscale('log')
                ax1.set_xscale('log')


        plt.ylabel(ylabel,size=22)
        plt.xlabel(xlabel,size=22)
        plt.legend(numpoints=1,fontsize=18,loc='best')


        if set_ylim==True:
            ax1.set_ylim(ylow,yhigh)
        else:
            print('Using default y-limit range for the plot: %s'%savelabel)
        fig.tight_layout()

        plt.savefig(savelabel+'.png',format='png',bbox_inches='tight',dpi=300)
        plt.show()
        
    def nPlot_variColor(self,pairList,labelList,colorList,savelabel,xlabel='No Label Provided',ylabel='No Label Provided',
              LogLin=True,LinLin=False,LogLog=False,linewidth=3,
              set_ylim=False,ylow=0.0001,yhigh=1,darkmode=False):
        '''
        :param pairList: list of lists (tuple), must be [[x1,y1],...[xn,yn]]
        :param labelList: list of length n, labeling the sets of tuples in pairList
        :param savelabel:
        :param xlabel:
        :param ylabel:
        :param linewidth:
        :return:
        '''

        if darkmode==True:
            plt.style.use('dark_background')
        else:
            mpl.rcParams.update(mpl.rcParamsDefault)

        fig=plt.figure(figsize=(10,8)) # set figure dimensions
        ax1=fig.add_subplot(1,1,1) # allows us to build more complex plots
        for tick in ax1.xaxis.get_major_ticks():
            tick.label1.set_fontsize(20) # scale for publication needs
            tick.label1.set_fontname('Avenir')
        for tick in ax1.yaxis.get_major_ticks():
            tick.label1.set_fontsize(20) # scale for publication needs
            tick.label1.set_fontname('Avenir')

        if LogLin == True and LinLin == True and LogLog == True: # kicks you out of the function if you set more then one mode to true
            print("Cannot set more than one mode equal to True")
            return
        elif LogLin == True and LinLin == True:
            print("Cannot set more than one mode equal to True")
            return
        elif LogLin == True and LogLog == True:
            print("Cannot set more than one mode equal to True")
            return
        elif LinLin == True and LogLog == True:
            print("Cannot set more than one mode equal to True")
            return

        n=0
        if LogLin==True:
            for i,j in zip(pairList,colorList):
                plt.semilogy(i[0],i[1],
                            label=labelList[n],
                            color=j,
                            linewidth=linewidth,linestyle='dashed')
                n+=1
        elif LinLin==True:
            for i in pairList:
                plt.plot(i[0],i[1],
                            label=labelList[n],
                            linewidth=linewidth,linestyle='dashed')
                n+=1
        elif LogLog==True:
            for i in pairList:
                plt.plot(i[0],i[1],
                            label=labelList[n],
                            linewidth=linewidth,linestyle='dashed')
                n+=1
                ax1.set_yscale('log')
                ax1.set_xscale('log')


        plt.ylabel(ylabel,size=22)
        plt.xlabel(xlabel,size=22)
        plt.legend(numpoints=1,fontsize=18,loc='best')


        if set_ylim==True:
            ax1.set_ylim(ylow,yhigh)
        else:
            print('Using default y-limit range for the plot: %s'%savelabel)
        fig.tight_layout()

        plt.savefig(savelabel+'.png',format='png',bbox_inches='tight',dpi=300)
        plt.show()
        
    def nPlot_variColor_and_Range(self,pairList,labelList,colorList,savelabel,xlabel='No Label Provided',ylabel='No Label Provided',
              LogLin=True,LinLin=False,LogLog=False,linewidth=3,
              set_ylim=False,ylow=0.0001,yhigh=1,qmin=0,qmax=0,darkmode=False):
        '''
        :param pairList: list of lists (tuple), must be [[x1,y1],...[xn,yn]]
        :param labelList: list of length n, labeling the sets of tuples in pairList
        :param savelabel:
        :param xlabel:
        :param ylabel:
        :param linewidth:
        :return:
        '''

        ptypes = ('-','--')

        if darkmode==True:
            plt.style.use('dark_background')
        else:
            mpl.rcParams.update(mpl.rcParamsDefault)

        if qmax==0:
            n=0
            for i in pairList:
                qmax=len(i[0])
                n+=1
                if n>=1:
                    break

        fig=plt.figure(figsize=(10,8)) # set figure dimensions
        ax1=fig.add_subplot(1,1,1) # allows us to build more complex plots
        for tick in ax1.xaxis.get_major_ticks():
            tick.label1.set_fontsize(20) # scale for publication needs
            tick.label1.set_fontname('Avenir')
        for tick in ax1.yaxis.get_major_ticks():
            tick.label1.set_fontsize(20) # scale for publication needs
            tick.label1.set_fontname('Avenir')

        if LogLin == True and LinLin == True and LogLog == True: # kicks you out of the function if you set more then one mode to true
            print("Cannot set more than one mode equal to True")
            return
        elif LogLin == True and LinLin == True:
            print("Cannot set more than one mode equal to True")
            return
        elif LogLin == True and LogLog == True:
            print("Cannot set more than one mode equal to True")
            return
        elif LinLin == True and LogLog == True:
            print("Cannot set more than one mode equal to True")
            return

        n=0
        if LogLin==True:
            for i,j in zip(pairList,colorList):
                plt.semilogy(i[0][qmin:qmax],i[1][qmin:qmax],
                            label=labelList[n],
                            color=j,
                            linewidth=linewidth,linestyle='dashed')
                plt.minorticks_off()
                n+=1
        elif LinLin==True:
            for i,j,z in zip(pairList, colorList,itertools.cycle(ptypes)):
                plt.plot(i[0],i[1],
                            label=labelList[n],
                            linewidth=linewidth,linestyle=z,
                            color=j)
                n+=1
        elif LogLog==True:
            for i in pairList:
                plt.plot(i[0],i[1],
                            label=labelList[n],
                            linewidth=linewidth,linestyle='dashed')
                n+=1
                ax1.set_yscale('log')
                ax1.set_xscale('log')


        plt.ylabel(ylabel,size=22)
        plt.xlabel(xlabel,size=22)
        plt.legend(numpoints=1,fontsize=18,loc='best')


        if set_ylim==True:
            ax1.set_ylim(ylow,yhigh)
        else:
            print('Using default y-limit range for the plot: %s'%savelabel)
        fig.tight_layout()

        plt.savefig(savelabel+'.png',format='png',bbox_inches='tight',dpi=300)
        plt.show()

    def nPlot_variX(self,pairList,labelList,savelabel,xlabel='No Label Provided',ylabel='No Label Provided',
              LogLin=True,LinLin=False,LogLog=False,linewidth=3,
              set_ylim=False,ylow=0.0001,yhigh=1,darkmode=False):
        '''
        :param pairList: list of lists (tuple), must be [[x1,y1],...[xn,yn]]
        :param labelList: list of length n, labeling the sets of tuples in pairList
        :param savelabel:
        :param xlabel:
        :param ylabel:
        :param linewidth:
        :return:
        '''

        if darkmode==True:
            plt.style.use('dark_background')
        else:
            mpl.rcParams.update(mpl.rcParamsDefault)

        fig=plt.figure(figsize=(10,8)) # set figure dimensions
        ax1=fig.add_subplot(1,1,1) # allows us to build more complex plots
        for tick in ax1.xaxis.get_major_ticks():
            tick.label1.set_fontsize(20) # scale for publication needs
            tick.label1.set_fontname('Avenir')
        for tick in ax1.yaxis.get_major_ticks():
            tick.label1.set_fontsize(20) # scale for publication needs
            tick.label1.set_fontname('Avenir')

        if LogLin == True and LinLin == True and LogLog == True: # kicks you out of the function if you set more then one mode to true
            print("Cannot set more than one mode equal to True")
            return
        elif LogLin == True and LinLin == True:
            print("Cannot set more than one mode equal to True")
            return
        elif LogLin == True and LogLog == True:
            print("Cannot set more than one mode equal to True")
            return
        elif LinLin == True and LogLog == True:
            print("Cannot set more than one mode equal to True")
            return

        n=0
        if LogLin==True:
            for i in pairList:
                plt.semilogy(i[0],i[1],
                            label=labelList[n],
                            linewidth=linewidth,
                             linestyle='dotted')
                n+=1
                plt.semilogy(i[2],i[3],
                            label=labelList[n],
                            linewidth=linewidth)
                n+=1
        elif LinLin==True:
            for i in pairList:
                plt.plot(i[0],i[1],
                            label=labelList[n],
                            linewidth=linewidth,
                             linestyle='dotted')
                n+=1
                plt.plot(i[2],i[3],
                            label=labelList[n],
                            linewidth=linewidth)
                n+=1
        elif LogLog==True:
            for i in pairList:
                plt.plot(i[0],i[1],
                            label=labelList[n],
                            linewidth=linewidth,
                             linestyle='dotted')
                n+=1
                plt.plot(i[2],i[3],
                            label=labelList[n],
                            linewidth=linewidth)
                n+=1
                ax1.set_yscale('log')
                ax1.set_xscale('log')


        plt.ylabel(ylabel,size=22)
        plt.xlabel(xlabel,size=22)
        plt.legend(numpoints=1,fontsize=18,loc='best')


        if set_ylim==True:
            ax1.set_ylim(ylow,yhigh)
        else:
            print('Using default y-limit range for the plot: %s'%savelabel)
        fig.tight_layout()

        plt.savefig(savelabel+'.png',format='png',bbox_inches='tight',dpi=300)
        plt.show()

    def nPlot_variX_and_Color(self,
        pairList,labelList,colorList,savelabel,
        xlabel='No Label Provided',ylabel='No Label Provided',
        LogLin=True,LinLin=False,LogLog=False, kratky=False,
        linewidth=3,
        set_ylim=False,ylow=0.0001,yhigh=1,darkmode=False,
        lg=True,
        lg_size=14,
        lg_loc = 'best',
        linestyles = ['-','-','-']):
        '''
        :param pairList: list of lists (tuple), must be [[x1,y1],...[xn,yn]]
        :param labelList: list of length n, labeling the sets of tuples in pairList
        :param savelabel:
        :param xlabel:
        :param ylabel:
        :param linewidth:
        :return:
        '''
        
        plt.rc('font',**{'sans-serif': ['Avenir']})
        if darkmode==True:
            plt.style.use('dark_background')
            c1='#EFECE8'
        else:
            mpl.rcParams.update(mpl.rcParamsDefault)
            c1='k'

        fig=plt.figure(figsize=(10,8)) # set figure dimensions
        ax1=fig.add_subplot(1,1,1) # allows us to build more complex plots
        for tick in ax1.xaxis.get_major_ticks():
            tick.label1.set_fontsize(20) # scale for publication needs
            tick.label1.set_fontname('Avenir')
        for tick in ax1.yaxis.get_major_ticks():
            tick.label1.set_fontsize(20) # scale for publication needs
            tick.label1.set_fontname('Avenir')

        if LogLin == True and LinLin == True and LogLog == True: # kicks you out of the function if you set more then one mode to true
            print("Cannot set more than one mode equal to True")
            return
        elif LogLin == True and LinLin == True:
            print("Cannot set more than one mode equal to True")
            return
        elif LogLin == True and LogLog == True:
            print("Cannot set more than one mode equal to True")
            return
        elif LinLin == True and LogLog == True:
            print("Cannot set more than one mode equal to True")
            return

        cycol = cycle(linestyles)    

        n=0
        if LogLin==True:
            for i,j in zip(pairList,colorList):
                plt.semilogy(i[0],i[1],
                            label=labelList[n],
                            linewidth=linewidth,
                            color=j,
                            linestyle=next(cycol))
                n+=1
        elif LinLin==True:
            for i,j,z in zip(pairList,colorList,labelList):
                plt.plot(i[0],i[1],
                            label=z,
                            linewidth=linewidth,
                            color=j,
                            linestyle=next(cycol))
                n+=1
        elif LogLog==True:
            for i in pairList:
                plt.plot(i[0],i[1],
                            label=labelList[n],
                            linewidth=linewidth,
                             linestyle='dotted')
                n+=1
                plt.plot(i[2],i[3],
                            label=labelList[n],
                            linewidth=linewidth)
                n+=1
                ax1.set_yscale('log')
                ax1.set_xscale('log')
        elif kratky==True:
          for i in pairList:
              plt.plot(i[0],i[1]*(i[0]**2),
                          label=labelList[n],
                          linewidth=linewidth,
                          linestyle=next(cycol))
              n+=1

        plt.ylabel(ylabel,size=22)
        plt.xlabel(xlabel,size=22)
        if lg == True:
          plt.legend(numpoints=1,fontsize=lg_size,loc=lg_loc)


        if set_ylim==True:
            ax1.set_ylim(ylow,yhigh)
        # else:
            # print('Using default y-limit range for the plot: %s'%savelabel)
        fig.tight_layout()

        plt.savefig(savelabel+'.png',format='png',bbox_inches='tight',dpi=300)
        plt.show(fig);

    def nPlot_4Panel(self,
        pairList_1,labelList_1,colorList_1,
        pairList_2,labelList_2,colorList_2,
        pairList_3,labelList_3,colorList_3,
        pairList_4,labelList_4,colorList_4,
        savelabel,
        xlabel='No Label Provided',ylabel='No Label Provided',
        LogLin=True,LinLin=False,LogLog=False,linewidth=3,
        set_ylim=False,ylow=0.0001,yhigh=1,darkmode=False,
        lg_size=14,
        lg_loc = 'best'):
        '''
        :param pairList: list of lists (tuple), must be [[x1,y1],...[xn,yn]]
        :param labelList: list of length n, labeling the sets of tuples in pairList
        :param savelabel:
        :param xlabel:
        :param ylabel:
        :param linewidth:
        :return:
        '''

        if darkmode==True:
            plt.style.use('dark_background')
            c1='#EFECE8'
        else:
            mpl.rcParams.update(mpl.rcParamsDefault)
            c1='k'

        fig=plt.figure(figsize=(10,8)) # set figure dimensions
        axs=fig.add_subplot(2,2) # allows us to build more complex plots
        for tick in axs.xaxis.get_major_ticks():
            tick.label1.set_fontsize(20) # scale for publication needs
            tick.label1.set_fontname('Avenir')
        for tick in axs.yaxis.get_major_ticks():
            tick.label1.set_fontsize(20) # scale for publication needs
            tick.label1.set_fontname('Avenir')


        cycol = cycle(['-','-','-'])    

        for i,j,z,a,b,c,d,e,f,w,y,x in zip(pairList_1,colorList_1,labelList_1,
            pairList_2,colorList_2,labelList_2,
            pairList_3,colorList_3,labelList_3,
            pairList_4,colorList_4,labelList_4):
            axs[0, 0].plot(i[0],i[1],
                        label=z,
                        linewidth=linewidth,
                        color=j,
                        linestyle=next(cycol))
            axs[0, 1].plot(a[0],a[1],
                        label=c,
                        linewidth=linewidth,
                        color=b,
                        linestyle=next(cycol))
            axs[1, 0].plot(d[0],d[1],
                        label=f,
                        linewidth=linewidth,
                        color=e,
                        linestyle=next(cycol))
            axs[1, 1].plot(w[0],w[1],
                        label=x,
                        linewidth=linewidth,
                        color=y,
                        linestyle=next(cycol))


        plt.ylabel(ylabel,size=22)
        plt.xlabel(xlabel,size=22)
        plt.legend(numpoints=1,fontsize=lg_size,loc=lg_loc)


        if set_ylim==True:
            ax1.set_ylim(ylow,yhigh)

        fig.tight_layout()

        plt.savefig(savelabel+'.png',format='png',bbox_inches='tight',dpi=300)
        plt.show(fig);

    def nPointPlot_variX_and_Color(self,pairList,labelList,colorList,savelabel,xlabel='No Label Provided',ylabel='No Label Provided',
              LogLin=True,LinLin=False,LogLog=False,
              linewidth=3,
              marker_size=10,
              set_ylim=False,
              ylow=0.0001,yhigh=1,
              darkmode=False,
              lg_size=12,
              leg=True):
        '''
        :param pairList: list of lists (tuple), must be [[x1,y1],...[xn,yn]]
        :param labelList: list of length n, labeling the sets of tuples in pairList
        :param savelabel:
        :param xlabel:
        :param ylabel:
        :param linewidth:
        :return:
        '''

        if darkmode==True:
            plt.style.use('dark_background')
            c1='#EFECE8'
        else:
            mpl.rcParams.update(mpl.rcParamsDefault)
            c1='k'

        fig=plt.figure(figsize=(10,8)) # set figure dimensions
        ax1=fig.add_subplot(1,1,1) # allows us to build more complex plots
        for tick in ax1.xaxis.get_major_ticks():
            tick.label1.set_fontsize(20) # scale for publication needs
            tick.label1.set_fontname('Avenir')
        for tick in ax1.yaxis.get_major_ticks():
            tick.label1.set_fontsize(20) # scale for publication needs
            tick.label1.set_fontname('Avenir')

        if LogLin == True and LinLin == True and LogLog == True: # kicks you out of the function if you set more then one mode to true
            print("Cannot set more than one mode equal to True")
            return
        elif LogLin == True and LinLin == True:
            print("Cannot set more than one mode equal to True")
            return
        elif LogLin == True and LogLog == True:
            print("Cannot set more than one mode equal to True")
            return
        elif LinLin == True and LogLog == True:
            print("Cannot set more than one mode equal to True")
            return

        n=0
        if LogLin==True:
            for i,j in zip(pairList,colorList):
                plt.semilogy(i[0],i[1],
                            label=labelList[n],
                            linestyle=None,
                            linewidth=linewidth,
                            marker = 'o',
                            markersize = linewidth + 5,
                            color=j)
                n+=1
                # plt.semilogy(i[2],i[3],
                #             label=labelList[n],
                #             linewidth=linewidth,
                #             color=c1)
                # n+=1
        elif LinLin==True:
            for i,j in zip(pairList,colorList):
              plt.plot(i[0],i[1],
              label=labelList[n],
              linestyle=None,
              linewidth=linewidth,
              marker = 'o',
              markersize = marker_size,
              color=j)
        elif LogLog==True:
            for i in pairList:
                plt.plot(i[0],i[1],
                            label=labelList[n],
                            linewidth=linewidth,
                             linestyle='dotted')
                n+=1
                plt.plot(i[2],i[3],
                            label=labelList[n],
                            linewidth=linewidth)
                n+=1
                ax1.set_yscale('log')
                ax1.set_xscale('log')


        plt.ylabel(ylabel,size=22)
        plt.xlabel(xlabel,size=22)

        if leg == True:
          plt.legend(numpoints=1,fontsize=lg_size,loc='best')


        if set_ylim==True:
            ax1.set_ylim(ylow,yhigh)
        # else:
            # print('Using default y-limit range for the plot: %s'%savelabel)
        fig.tight_layout()

        plt.savefig(savelabel+'.png',format='png',bbox_inches='tight',dpi=300)
        plt.show()

    def nPointPlot_errorBars(self,pairList,labelList,colorList,savelabel,xlabel='No Label Provided',ylabel='No Label Provided',
              LogLin=True,LinLin=False,LogLog=False,
              linewidth=3,
              marker_size=10,
              set_ylim=False,
              ylow=0.0001,yhigh=1,
              darkmode=False,
              lg=True,
              lg_size=10,
              elinewidth=1,
              barsize=1,
              capsize=1):
        '''
        :param pairList: list of lists (tuple), must be [[x1,y1],...[xn,yn]]
        :param labelList: list of length n, labeling the sets of tuples in pairList
        :param savelabel:
        :param xlabel:
        :param ylabel:
        :param linewidth:
        :return:
        '''

        if darkmode==True:
            plt.style.use('dark_background')
            c1='#EFECE8'
        # else:
        #     mpl.rcParams.update(mpl.rcParamsDefault)
        #     c1='k'

        fig=plt.figure(figsize=(8,6)) # set figure dimensions
        ax1=fig.add_subplot(1,1,1) # allows us to build more complex plots
        ax1.tick_params(width=4,length=8)
        for tick in ax1.xaxis.get_major_ticks():
            tick.label1.set_fontsize(16) # scale for publication needs
            tick.label1.set_fontname('Avenir')
        for tick in ax1.yaxis.get_major_ticks():
            tick.label1.set_fontsize(16) # scale for publication needs
            tick.label1.set_fontname('Avenir')

        if LogLin == True and LinLin == True and LogLog == True: # kicks you out of the function if you set more then one mode to true
            print("Cannot set more than one mode equal to True")
            return
        elif LogLin == True and LinLin == True:
            print("Cannot set more than one mode equal to True")
            return
        elif LogLin == True and LogLog == True:
            print("Cannot set more than one mode equal to True")
            return
        elif LinLin == True and LogLog == True:
            print("Cannot set more than one mode equal to True")
            return

        n=0
        if LogLin==True:
            for i,j in zip(pairList,colorList):
                plt.semilogy(i[0],i[1],
                            label=labelList[n],
                            linestyle=None,
                            linewidth=linewidth,
                            marker = 'o',
                            markersize = linewidth + 5,
                            color=j)
                n+=1
        elif LinLin==True:
            n=0
            for i,j in zip(pairList,colorList):
              plt.plot(i[0],i[1],
              label=labelList[n],
              linestyle=None,
              linewidth=linewidth,
              marker = 'o',
              markersize = marker_size,
              color=j)
              plt.errorbar(x=i[0],
                y=i[1],
                yerr=i[2],
                xerr=None,
                label=None,
                elinewidth=elinewidth,
                capsize=capsize,
                capthick=barsize,
                color=j)
              n+=1
        elif LogLog==True:
            for i in pairList:
                plt.plot(i[0],i[1],
                            label=labelList[n],
                            linewidth=linewidth,
                             linestyle='dotted')
                n+=1
                plt.plot(i[2],i[3],
                            label=labelList[n],
                            linewidth=linewidth)
                n+=1
                ax1.set_yscale('log')
                ax1.set_xscale('log')


        plt.ylabel(ylabel,size=18)
        plt.xlabel(xlabel,size=18)
        if lg == True:
            plt.legend(numpoints=1,fontsize=lg_size,loc='best')


        if set_ylim==True:
            ax1.set_ylim(ylow,yhigh)
        # else:
            # print('Using default y-limit range for the plot: %s'%savelabel)
        fig.tight_layout()

        plt.savefig(savelabel+'.png',format='png',bbox_inches='tight',dpi=300)
        plt.show()



    def IFT_plot(self,IFT,savelabel1='tkRubisCO_0MPa_GNOM_PDDF',savelabel2='RegularizedFit_GNOM',
                 plotlabel1='Pair Distance Distribution',plotlabel2='Baseline',plotlabel3='Expt',plotlabel4='Regularized Fit',
                 darkmode=False):
        '''
        Plot the output of DATGNOM calculation
        '''

        Pr,R,Pr_err,Jexp,qshort,Jerr,Jreg,results,Ireg,qfull = IFT[0],IFT[1],IFT[2],IFT[3],IFT[4],IFT[5],IFT[6],IFT[7], \
                                                               IFT[8],IFT[9]
        self.twoPlot(X=R,Y1=Pr,Y2=[0] * len(Pr),savelabel=savelabel1,
                           plotlabel1=plotlabel1,plotlabel2=plotlabel2,
                           xlabel='r($\\AA$)',ylabel='P(r)',linewidth=4,darkmode=darkmode)
        self.twoPlot_variX(X1=qshort,Y1=Jexp,X2=qshort,Y2=Jreg,plotlabel1=plotlabel3,plotlabel2=plotlabel4,
                                 savelabel=savelabel2,xlabel='q $\\rm \\AA^{-1}$',ylabel='I(q)',LogLin=True,darkmode=darkmode)

    def vertical_stackPlot(self,X1=[],Y1=[],Y1err=[],X2=[],Y2=[],ylabel1='No label provided',ylabel2='No label provided',xlabel='No label provided',
                           Label1='',
                           Label2='',saveLabel='Vertical_Residuals',bottomPlot_yLabel='$\\rm ln(\\frac{I_{expt}(q)}{I_{model}(q)}) \cdot (\\frac{1}{\sigma_{expt}})$',
                           LinLin=True,linewidth=4,labelSize=20,darkmode=False,plot=True):
        '''
        X1:
        X2:
        Y1:
        Y1err:
        Y2:
        ylabel1:
        ylabel2:
        xlabel:
        Label1:
        Label2:
        saveLabel:
        bottomPlot_yLabel:


        Second plot is a residual plot, therefore X2/Y2 should be the model.. Also, X1 must = x2
        '''

        if darkmode==True:
            plt.style.use('dark_background')
            c1='#D2700F'
            c2='#22A2AA'
        else:
            mpl.rcParams.update(mpl.rcParamsDefault)
            c1='k'
            c2='#F17800'

        plt.rcParams['xtick.major.pad'] = 10
        plt.rcParams['ytick.major.pad'] = 10
        plt.rcParams['axes.linewidth'] = 2
        # plt.rcParams["figure.figsize"] = (20,10)
        fg = plt.figure(figsize=(10,7))
        ax = plt.subplot2grid((3,3),(0,0),rowspan=2,colspan=3)
        ax2 = plt.subplot2grid((3,3),(2,0),rowspan=1,colspan=3)
        plt.rc("axes",linewidth=linewidth)
        plt.rc('font',**{"sans-serif":["Avenir"]})
        if LinLin==True:
            ax.plot(X1,Y1,
                color=c1,
                marker='o',
                markersize=linewidth,
                linestyle='None',
                label=Label1)
            ax.plot(X2,Y2,
                color=c2,
                linestyle='-',
                linewidth=linewidth,
                label=Label2)
        else:
            ax.semilogy(X1,Y1,
                    color=c1,
                    linestyle='-',
                    linewidth=linewidth,
                    label=Label1)
            ax.semilogy(X2,Y2,
                    color=c2,
                    linestyle='-',
                    linewidth=linewidth,
                    label=Label2)
        ax.set_ylabel(ylabel1,size=labelSize)
        ax.legend(numpoints=1,fontsize=18,loc="best")
        ax.xaxis.set_tick_params(which='both',width=2)
        ax.set_xticklabels([])
        ax.yaxis.set_tick_params(which='both',width=2)
        for tick in ax.yaxis.get_major_ticks():
            tick.label1.set_fontsize(labelSize)
            tick.label1.set_fontname('Avenir')
        ax2.plot(X1,((Y1 - Y2)) / ((Y1err)),
                 color=c1,
                 linestyle='-',
                 linewidth=linewidth-1,
                 label=ylabel2)
        ax2.plot(X1,[0] * len(Y1),
                 color=c2,
                 linestyle='--',
                 linewidth=linewidth-1)
        ax2.set_xlabel(xlabel,size=labelSize) # 'q = $\\frac{4 \pi sin(\\theta)}{\\lambda}$ ($\\AA^{-1}$)'
        ax2.set_ylabel(bottomPlot_yLabel,size=labelSize)
        # $\\frac{\\frac{ln(I_{expt}(q))}{ln(I_{model}(q))}}{\sigma_{expt}}$
        # \\frac{ln(I_{expt}(q))}{ln(I_{model}(q))} \cdot \\frac{1}{\sigma_{expt}}
        ax2.legend(numpoints=1,fontsize=18,loc="best")
        # ax2.set_xlim(q[nmin],0.71)
        # ax.set_ylim(0.01,0.025)
        ax2.xaxis.set_tick_params(which='both',width=2)
        ax2.yaxis.set_tick_params(which='both',width=2)
        for tick in ax2.xaxis.get_major_ticks():
            tick.label1.set_fontsize(labelSize)
            tick.label1.set_fontname('Avenir')
        for tick in ax2.yaxis.get_major_ticks():
            tick.label1.set_fontsize(labelSize)
            tick.label1.set_fontname('Avenir')
        ax2.yaxis.set_major_locator(plt.MaxNLocator(5))
        fg.tight_layout()
        plt.savefig(saveLabel + 'VerticalStack_plot.png',
                    format='png',dpi=200,bbox_inches='tight')
        if plot==True:
            plt.show()


    def vertical_stackPlot_Guiner(self,X1=[],Y1=[],Y1err=[],X2=[],Y2=[],idxmin=0,idxmax=50,ylabel1='No label provided',ylabel2='No label provided',xlabel='No label provided',
                           Label1='',
                           Label2='',saveLabel='Vertical_Residuals',
                           bottomPlot_yLabel='$ln(\\frac{I_{expt}(q)}{I_{model}(q)}) \cdot (\\frac{1}{\sigma_{expt}})$',
                           LinLin=True,linewidth=4,
                           labelSize=16,
                           darkmode=False,
                           dpi=100,
                           plot=True):
        '''
        X1: q
        X2: q
        Y1: Raw Data
        Y1err:
        Y2: Model
        ylabel1:
        ylabel2:
        xlabel:
        Label1:
        Label2:
        saveLabel:
        bottomPlot_yLabel:


        Second plot is a residual plot, therefore X2/Y2 should be the model.. Also, X1 must = x2
        '''



        Y1_org=np.log(Y1[idxmin:idxmax])
        Y1err_org=np.log(Y1err[idxmin:idxmax])
        # X1=X1[idxmin:idxmax]**2
        X1_org=X2[idxmin:idxmax]**2
        X2=X2[idxmin:idxmax]**2
        # Y2=Y2[idxmin:idxmax]

        if idxmin >= 15 and idxmax <= len(Y1)-15:
            X1=X1[idxmin-15:idxmax+15]**2
            Y1=np.log(Y1[idxmin-15:idxmax+15])
            Y1err=np.log(Y1err[idxmin-15:idxmax+15])
        elif idxmin <= 15 and idxmax <= len(Y1)-15:
            X1=X1[0:idxmax+15]**2
            Y1=np.log(Y1[0:idxmax+15])
            Y1err=np.log(Y1err[0:idxmax+15])

        if darkmode==True:
            plt.style.use('dark_background')
            c1='#D2700F'
            c2='#22A2AA'
        else:
            mpl.rcParams.update(mpl.rcParamsDefault)
            c1='k'
            c2='#0E8987'

        plt.rcParams['xtick.major.pad'] = 10
        plt.rcParams['ytick.major.pad'] = 10
        plt.rcParams['axes.linewidth'] = 2
        # plt.rcParams["figure.figsize"] = (20,10)
        fg = plt.figure(figsize=(10,6))
        ax = plt.subplot2grid((3,3),(0,0),rowspan=2,colspan=3)
        ax2 = plt.subplot2grid((3,3),(2,0),rowspan=1,colspan=3)
        plt.rc("axes",linewidth=linewidth)
        plt.rc('font',**{"sans-serif":["Avenir"]})
        if LinLin==True:
            ax.plot(X1,Y1,
                color=c1,
                marker='o',
                markersize=linewidth,
                linestyle='None',
                label=Label1)
            ax.plot(X2,Y2,
                color=c2,
                linestyle='-',
                linewidth=linewidth,
                label=Label2)
        else:
            ax.semilogy(X1,Y1,
                    color=c1,
                    marker='o',
                    markersize=linewidth,
                    linestyle='None',
                    label=Label1)
            ax.semilogy(X2,Y2,
                    color=c2,
                    linestyle='-',
                    linewidth=linewidth,
                    label=Label2)

        ax.axvline(x=X2[0],
            linestyle='--',
            linewidth=4,
            color='#DA2700')
        ax.axvline(x=X2[-1],
            linestyle='--',
            linewidth=4,
            color='#DA2700')
        ax.set_ylabel(ylabel1,size=labelSize)
        ax.legend(numpoints=1,fontsize=14,loc="best")
        ax.xaxis.set_tick_params(which='both',width=2)
        ax.set_xticklabels([])
        ax.yaxis.set_tick_params(which='both',width=2)
        for tick in ax.yaxis.get_major_ticks():
            tick.label1.set_fontsize(labelSize)
            tick.label1.set_fontname('Avenir')
        ax2.plot(X1_org,((Y1_org - Y2)) / ((Y1err_org)),
                 color=c1,
                 linestyle='-',
                 linewidth=linewidth-1,
                 label=ylabel2)
        ax2.plot(X1_org,[0] * len(Y1_org),
                 color=c2,
                 linestyle='--',
                 linewidth=linewidth-1)
        ax2.set_xlabel(xlabel,size=labelSize) # 'q = $\\frac{4 \pi sin(\\theta)}{\\lambda}$ ($\\AA^{-1}$)'
        ax2.set_ylabel(bottomPlot_yLabel,size=labelSize)
        # $\\frac{\\frac{ln(I_{expt}(q))}{ln(I_{model}(q))}}{\sigma_{expt}}$
        # \\frac{ln(I_{expt}(q))}{ln(I_{model}(q))} \cdot \\frac{1}{\sigma_{expt}}
        ax2.legend(numpoints=1,fontsize=14,loc="best")
        # ax2.set_xlim(q[nmin],0.71)
        # ax.set_ylim(0.01,0.025)
        ax2.xaxis.set_tick_params(which='both',width=2)
        ax2.yaxis.set_tick_params(which='both',width=2)
        for tick in ax2.xaxis.get_major_ticks():
            tick.label1.set_fontsize(labelSize)
            tick.label1.set_fontname('Avenir')
        for tick in ax2.yaxis.get_major_ticks():
            tick.label1.set_fontsize(labelSize)
            tick.label1.set_fontname('Avenir')
        ax2.yaxis.set_major_locator(plt.MaxNLocator(5))
        fg.tight_layout()
        plt.savefig(saveLabel + 'VerticalStack_plot.png',
                    format='png',dpi=dpi,bbox_inches='tight')
        if plot==True:
            plt.show();
        else:
            plt.close(fg);

    def plotSAXS_Chromatogram(self, q, I, Frame,
                              saveLabel,
                              ylabel,
                              xlabel,
                              zlabel,
                              linewidth=3,
                              darkmode=False):
        '''
        Describe function here

        Getting I(q) for each frame might actually be tricky... 
        (1) Read config file
        (2) Apply to/read H5 files
        (3) Store them..
        (4) Plot them
        '''

        if darkmode==True:
            plt.style.use('dark_background')
            c1='#EFECE8'
        else:
            mpl.rcParams.update(mpl.rcParamsDefault)
            c1='k'

        fig=plt.figure(figsize=(10,8)) # set figure dimensions
        # ax1=fig.add_subplot(1,1,1) # allows us to build more complex plots
        ax1 = plt.axes(projection='3d')
        # for tick in ax1.xaxis.get_major_ticks():
        #     tick.label1.set_fontsize(20) # scale for publication needs
        #     tick.label1.set_fontname('Avenir')
        # for tick in ax1.yaxis.get_major_ticks():
        #     tick.label1.set_fontsize(20) # scale for publication needs
        #     tick.label1.set_fontname('Avenir')


        
        ax.plot_surface(q,I,Frame, rstride=1, cstride=1,
            cmap='viridis', edgecolor='none')

        ax.set_xlabel('x')
        ax.set_ylabel('y')
        ax.set_zlabel('z');

        # plt.ylabel(ylabel,size=22)
        # plt.xlabel(xlabel,size=22)
        # plt.legend(numpoints=1,fontsize=18,loc='best')

        fig.tight_layout()

        plt.savefig(savelabel+'.png',format='png',bbox_inches='tight',dpi=300)
        plt.show()

    def dual_yAxis(self,
                   pairList,
                   colorList,
                   labelList,
                   set_ylim1=False,
                   set_ylim2=False,
                   xlabel='X_Label',
                   ylabel1='Y_Label_#1',
                   ylabel2='Y_Label_#2',
                   linewidth=3,
                   savelabel='dual_yAxis',
                   markerLabel=None,
                   legend=None):
        '''
        doc string


        '''

        fig=plt.figure(figsize=(10,8)) # set figure dimensions
        ax1=fig.add_subplot(1,1,1) # allows us to build more complex plots
        for tick in ax1.xaxis.get_major_ticks():
            tick.label1.set_fontsize(20) # scale for publication needs
            tick.label1.set_fontname('Avenir')
        for tick in ax1.yaxis.get_major_ticks():
            tick.label1.set_fontsize(20) # scale for publication needs
            tick.label1.set_fontname('Avenir')

        cycol_color = cycle(['#0E403E','#194E05','#57B036','#AA8B0B','#0F8985'])
        cycol_linestyle = cycle(['solid','dotted','dashed','dashdot','dashdotdotted'])
        cycol_markerstyle = cycle(['o','P','X','^','s'])

        lns = {}
        c=0
        n=0
        for i,j,z in zip(pairList,colorList,labelList):
            lns['%s'%str(n)]=ax1.plot(i[0],i[1],
                        label=z,
                        linewidth=linewidth,
                        linestyle='-',
                        color=colorList[c])
            n+=1
            c+=1
        
        ax1.set_ylabel(ylabel1,size=22)
        ax1.set_xlabel(xlabel,size=22)

        ax2=ax1.twinx()

        for tick in ax2.xaxis.get_major_ticks():
            tick.label2.set_fontsize(20) # scale for publication needs
            tick.label2.set_fontname('Avenir')
        for tick in ax2.yaxis.get_major_ticks():
            tick.label2.set_fontsize(20) # scale for publication needs
            tick.label2.set_fontname('Avenir')

        if markerLabel == None:
          for i,j,z in zip(pairList,colorList,labelList):
            lns['%s'%str(n)]=ax2.plot(i[0],i[2],
                        label='',
                        linestyle='None',
                        markersize=linewidth,
                        marker=next(cycol_markerstyle),
                        color=colorList[c])
            n+=1
            c+=1
        else:
          for i,j,z in zip(pairList,colorList,labelList):
              lns['%s'%str(n)]=ax2.plot(i[0],i[2],
                          label=z,
                          linestyle='None',
                          markersize=linewidth,
                          marker=next(cycol_markerstyle),
                          color=colorList[c])
              n+=1
              c+=1

        ax2.set_ylabel(ylabel2,size=22)

        all_lns = []
        all_labs = []
        for key,value in lns.items():
          all_labs.append(value[0].get_label())
          all_lns.append(value[0])

        if legend == None:
          ax1.legend(all_lns, all_labs, loc=0, fontsize = 18)
        else:
          ax1.legend([(all_lns)], [all_labs], numpoints=2,
               handler_map={tuple: HandlerTuple(ndivide=None)})

        # if set_ylim==True:
        #     ax1.set_ylim(ylow,yhigh)
        # else:
            # print('Using default y-limit range for the plot: %s'%savelabel)
        fig.tight_layout()

        plt.savefig(savelabel+'.png',format='png',bbox_inches='tight',dpi=300)
        plt.show()














