#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jun  4 14:47:22 2020

@author: Rob
"""

from matplotlib import pyplot as plt

class PlotClass:
    
    def __init__(self,notify=True):
        '''
        Describe the class here.
        The initialization of the class sets a few global parameters
        (1) Axes width
        (2) Line width
        (3) Set font to Helvetica. Note, if when running locally you receive a traceback error
            along the lines of "falling back to Arial" you need to enable the base Helvetica font
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
        self.font=plt.rc('font',**{'sans-serif': ['Helvetica']})
        
        
    def basicPlot(self,X,Y,plotlabel='',savelabel='',xlabel='',ylabel='NOT PROVIDED'):
        '''
        As simple as it gets
        '''
        fig=plt.figure(figsize=(10,8)) # set figure dimensions
        ax1=fig.add_subplot(1,1,1) # allows us to build more complex plots
        for tick in ax1.xaxis.get_major_ticks():
            tick.label1.set_fontsize(20) # scale for publication needs
            tick.label1.set_fontname('Helvetica')
        for tick in ax1.yaxis.get_major_ticks():
            tick.label1.set_fontsize(20) # scale for publication needs
            tick.label1.set_fontname('Helvetica')
        plt.ylabel(ylabel,size=22)
        plt.xlabel(xlabel,size=22)
        plt.plot(X,Y,'-',label=plotlabel,
                 color='#E55334') # best to use HTML color codes: https://htmlcolorcodes.com
        plt.legend(numpoints=1,fontsize=18,loc='best')
        fig.tight_layout()
        plt.savefig(savelabel+'.png',format='png',
                    bbox_inches='tight',dpi=300)
        plt.show()
        
    def semilogyPlot(self,X,Y,plotlabel='',savelabel='',xlabel='',ylabel='',linewidth=4):
        '''
        As simple as it gets
        '''
        fig=plt.figure(figsize=(10,8)) # set figure dimensions
        ax1=fig.add_subplot(1,1,1) # allows us to build more complex plots
        for tick in ax1.xaxis.get_major_ticks():
            tick.label1.set_fontsize(20) # scale for publication needs
            tick.label1.set_fontname('Helvetica')
        for tick in ax1.yaxis.get_major_ticks():
            tick.label1.set_fontsize(20) # scale for publication needs
            tick.label1.set_fontname('Helvetica')
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
        
        
    def twoPlot(self,X,Y1,Y2,plotlabel1='',plotlabel2='',savelabel='',xlabel='',ylabel='',linewidth=2):
        '''
        As simple as it gets P2
        '''
        fig=plt.figure(figsize=(10,8)) # set figure dimensions
        ax1=fig.add_subplot(1,1,1) # allows us to build more complex plots
        for tick in ax1.xaxis.get_major_ticks():
            tick.label1.set_fontsize(20) # scale for publication needs
            tick.label1.set_fontname('Helvetica')
        for tick in ax1.yaxis.get_major_ticks():
            tick.label1.set_fontsize(20) # scale for publication needs
            tick.label1.set_fontname('Helvetica')
        plt.ylabel(ylabel,size=22)
        plt.xlabel(xlabel,size=22)
        plt.plot(X,Y1,'-',label=plotlabel1,
                 linewidth=linewidth+2,
                 color='#E55334')
        plt.plot(X,Y2,'o',label=plotlabel2,
                 color='#1283BC')
        plt.legend(numpoints=1,fontsize=18,loc='best')
        fig.tight_layout()
        plt.savefig(savelabel+'.png',format='png',
                    bbox_inches='tight',dpi=300)
        plt.show()

    def twoPlot_variX(self,X1,Y1,X2,Y2,plotlabel1='',plotlabel2='',savelabel='',xlabel='',ylabel='',linewidth=2,
                      LogLin=False,set_ylim=False,ylow=0.0001,yhigh=1):
        '''
        As simple as it gets P2
        '''
        if LogLin==False:
            fig=plt.figure(figsize=(10,8)) # set figure dimensions
            ax1=fig.add_subplot(1,1,1) # allows us to build more complex plots
            for tick in ax1.xaxis.get_major_ticks():
                tick.label1.set_fontsize(20) # scale for publication needs
                tick.label1.set_fontname('Helvetica')
            for tick in ax1.yaxis.get_major_ticks():
                tick.label1.set_fontsize(20) # scale for publication needs
                tick.label1.set_fontname('Helvetica')
            plt.ylabel(ylabel,size=22)
            plt.xlabel(xlabel,size=22)
            plt.plot(X1,Y1,'-',label=plotlabel1,
                    linewidth=linewidth,
                    color='#E55334')
            plt.plot(X2,Y2,'o',label=plotlabel2,
                    color='#1283BC')
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
            fig=plt.figure(figsize=(10,8)) # set figure dimensions
            ax1=fig.add_subplot(1,1,1) # allows us to build more complex plots
            for tick in ax1.xaxis.get_major_ticks():
                tick.label1.set_fontsize(20) # scale for publication needs
                tick.label1.set_fontname('Helvetica')
            for tick in ax1.yaxis.get_major_ticks():
                tick.label1.set_fontsize(20) # scale for publication needs
                tick.label1.set_fontname('Helvetica')
            plt.ylabel(ylabel,size=22)
            plt.xlabel(xlabel,size=22)
            plt.semilogy(X1,Y1,'-',label=plotlabel1,
                    linewidth=linewidth,
                    color='#E55334')
            plt.semilogy(X2,Y2,'o',label=plotlabel2,
                    color='#1283BC')
            plt.legend(numpoints=1,fontsize=18,loc='best')
            if set_ylim == True:
                ax1.set_ylim(ylow,yhigh)
            # else:
            #     print('Using default y-limit range')
            fig.tight_layout()
            plt.savefig(savelabel+'.png',format='png',
                    bbox_inches='tight',dpi=300)
            plt.show()


    def nPlot(self,pairList,labelList,savelabel,xlabel='No Label Provided',ylabel='No Label Provided',
              LogLin=True,LinLin=False,LogLog=False,linewidth=3,
              set_ylim=False,ylow=0.0001,yhigh=1):
        '''
        :param pairList: list of lists (tuple), must be [[x1,y1],...[xn,yn]]
        :param labelList: list of length n, labeling the sets of tuples in pairList
        :param savelabel:
        :param xlabel:
        :param ylabel:
        :param linewidth:
        :return:
        '''

        fig=plt.figure(figsize=(10,8)) # set figure dimensions
        ax1=fig.add_subplot(1,1,1) # allows us to build more complex plots
        for tick in ax1.xaxis.get_major_ticks():
            tick.label1.set_fontsize(20) # scale for publication needs
            tick.label1.set_fontname('Helvetica')
        for tick in ax1.yaxis.get_major_ticks():
            tick.label1.set_fontsize(20) # scale for publication needs
            tick.label1.set_fontname('Helvetica')

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

    def nPlot_variX(self,pairList,labelList,savelabel,xlabel='No Label Provided',ylabel='No Label Provided',
              LogLin=True,LinLin=False,LogLog=False,linewidth=3,
              set_ylim=False,ylow=0.0001,yhigh=1):
        '''
        :param pairList: list of lists (tuple), must be [[x1,y1],...[xn,yn]]
        :param labelList: list of length n, labeling the sets of tuples in pairList
        :param savelabel:
        :param xlabel:
        :param ylabel:
        :param linewidth:
        :return:
        '''

        fig=plt.figure(figsize=(10,8)) # set figure dimensions
        ax1=fig.add_subplot(1,1,1) # allows us to build more complex plots
        for tick in ax1.xaxis.get_major_ticks():
            tick.label1.set_fontsize(20) # scale for publication needs
            tick.label1.set_fontname('Helvetica')
        for tick in ax1.yaxis.get_major_ticks():
            tick.label1.set_fontsize(20) # scale for publication needs
            tick.label1.set_fontname('Helvetica')

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


    def IFT_plot(self,IFT,savelabel1='tkRubisCO_0MPa_GNOM_PDDF',savelabel2='RegularizedFit_GNOM',
                 plotlabel1='Pair Distance Distribution',plotlabel2='Baseline',plotlabel3='Expt',plotlabel4='Regularized Fit'):
        '''
        Plot the output of DATGNOM calculation
        '''
        Pr,R,Pr_err,Jexp,qshort,Jerr,Jreg,results,Ireg,qfull = IFT[0],IFT[1],IFT[2],IFT[3],IFT[4],IFT[5],IFT[6],IFT[7], \
                                                               IFT[8],IFT[9]
        self.twoPlot(X=R,Y1=Pr,Y2=[0] * len(Pr),savelabel=savelabel1,
                           plotlabel1=plotlabel1,plotlabel2=plotlabel2,
                           xlabel='r($\\AA$)',ylabel='P(r)',linewidth=4)
        self.twoPlot_variX(X1=qshort,Y1=Jexp,X2=qshort,Y2=Jreg,plotlabel1=plotlabel3,plotlabel2=plotlabel4,
                                 savelabel=savelabel2,xlabel='q $\\AA^{-1}$',ylabel='I(q)',LogLin=True)

    def vertical_stackPlot(self,X1=[],Y1=[],Y1err=[],X2=[],Y2=[],ylabel1='No label provided',ylabel2='No label provided',xlabel='No label provided',
                           Label1='',
                           Label2='',saveLabel='Vertical_Residuals',bottomPlot_yLabel='$ln(\\frac{I_{expt}(q)}{I_{model}(q)}) \cdot (\\frac{1}{\sigma_{expt}})$'):
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
        plt.rcParams['xtick.major.pad'] = 10
        plt.rcParams['ytick.major.pad'] = 10
        plt.rcParams['axes.linewidth'] = 2
        fg = plt.figure(figsize=(15,12))
        ax = plt.subplot2grid((3,3),(0,0),rowspan=2,colspan=3)
        ax2 = plt.subplot2grid((3,3),(2,0),rowspan=1,colspan=3)
        plt.rc("axes",linewidth=2)
        plt.rc('font',**{"sans-serif":["Helvetica"]})
        ax.plot(X1,Y1,
                color='k',
                marker='o',
                markersize=3,
                linestyle='None',
                label=Label1)
        ax.plot(X2,Y2,
                color='#7F817F',
                linestyle='-',
                linewidth=3,
                label=Label2)
        ax.set_ylabel(ylabel1,size=20)
        ax.legend(numpoints=1,fontsize=18,loc="best")
        ax.xaxis.set_tick_params(which='both',width=2)
        ax.set_xticklabels([])
        ax.yaxis.set_tick_params(which='both',width=2)
        for tick in ax.yaxis.get_major_ticks():
            tick.label1.set_fontsize(20)
            tick.label1.set_fontname('Helvetica')
        ax2.plot(X1,((Y1 - Y2)) / ((Y1err)),
                 color='k',
                 linestyle='-',
                 linewidth=2,
                 label=ylabel2)
        ax2.plot(X1,[0] * len(Y1),
                 color='k',
                 linestyle='--',
                 linewidth=2)
        ax2.set_xlabel(xlabel,size=20) # 'q = $\\frac{4 \pi sin(\\theta)}{\\lambda}$ ($\\AA^{-1}$)'
        ax2.set_ylabel(bottomPlot_yLabel,size=20)
        # $\\frac{\\frac{ln(I_{expt}(q))}{ln(I_{model}(q))}}{\sigma_{expt}}$
        # \\frac{ln(I_{expt}(q))}{ln(I_{model}(q))} \cdot \\frac{1}{\sigma_{expt}}
        ax2.legend(numpoints=1,fontsize=18,loc="best")
        # ax2.set_xlim(q[nmin],0.71)
        # ax.set_ylim(0.01,0.025)
        ax2.xaxis.set_tick_params(which='both',width=2)
        ax2.yaxis.set_tick_params(which='both',width=2)
        for tick in ax2.xaxis.get_major_ticks():
            tick.label1.set_fontsize(20)
            tick.label1.set_fontname('Helvetica')
        for tick in ax2.yaxis.get_major_ticks():
            tick.label1.set_fontsize(20)
            tick.label1.set_fontname('Helvetica')
        ax2.yaxis.set_major_locator(plt.MaxNLocator(5))
        fg.tight_layout()
        plt.savefig(saveLabel + 'VerticalStack_plot.png',
                    format='png',dpi=500,bbox_inches='tight')
        plt.show()




