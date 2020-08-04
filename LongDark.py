import numpy as np
import matplotlib.pyplot as plt
import matplotlib
matplotlib.use('Agg')
from scipy import integrate
import AnaUtils as au
import os, glob

wdir = '/media/easystore/KA0206_LongDR_2020Q2'
fdirs = glob.glob(wdir+'/*')
fdirs.sort()
fdirT = [os.stat(idir).st_mtime for idir in fdirs]
#fdirs,fdirT

def grabTimes(mdir,fdirT):
    runF = [os.path.splitext(iF)[0] for iF in glob.glob(mdir+'/*.bin')]
    runsort = [os.path.basename(iF).zfill(9) for iF in runF]
    runF = [x for _,x in sorted(zip(runsort,runF))]
    runT = np.array([os.stat(iFi+'.bin').st_mtime for iFi in runF])
    runT = runT - runT[-1] + fdirT
    print(mdir+' : '+str(len(runT)))
    return runF,runT

allRuns = [grabTimes(fdirs[iR],fdirT[iR]) for iR in range(len(fdirs))]
allF = np.array([iRun[0] for iRun in allRuns]).flatten()
allT = np.array([iRun[1] for iRun in allRuns]).flatten()
allT -= allT[0]

def DRun(fname,nbase=150,winS=5,winF=5,cut=4,pmt=0,qbins=1000,ret=False,plot=False):
    waves = au.ReadDDC10_BinWave(fname)
    waves[0],base = au.Subtract_Baseline(waves[0],nBase=nbase)
    #require baseline has no pulse. i.e. integral over baseline less than cut*rms
    bmask = np.absolute(integrate.simps(waves[0][:,pmt,:nbase]))<cut*integrate.simps(np.ones(nbase))*base[1][:,pmt]
    PromptPeak = np.argmax(waves[0][:,pmt,:].sum(axis=0))
    wStart = PromptPeak-winS
    wFin = PromptPeak+winF
    Qhist = au.winQHist(waves,ch=pmt,init=wStart,end=wFin,nBins=qbins,hrange=[-3,10],evMask=bmask[...,np.newaxis])
    if plot:
        #au.plotWaves(waves[0],pmt,500)
        plt.clf()
        plt.errorbar(Qhist['qHist'][1],Qhist['qHist'][0],marker='+',yerr=np.sqrt(Qhist['qHist'][2]))
        plt.xlabel('Q [pC]')
        plt.ylabel('Rate [Hz]')
        plt.savefig(fname+"_QHist.png")
        #plt.show()
    
    if ret:
        Qhist['waves'] = waves[0][:,pmt]
        Qhist['evMask'] = bmask
        Qhist['baserms'] = base[0][:,pmt]
        return Qhist,waves[1]
    else:
        return Qhist['qHist']

from scipy.optimize import curve_fit

def fitQ(Qhist,P,bounds=(-np.inf,np.inf),doErr=False):
    def gauss(x, x0, y0, sigma):
        p = [x0, y0, sigma]
        return p[1]* np.exp(-((x-p[0])/p[2])**2)
    def gauss2(x,x0,y0,s0,x1,y1,s1):
        p0 = gauss(x,x0,y0,s0)
        p1 = gauss(x,x1,y1,s1)
        return p0+p1
    def gauss3(x,x0,y0,s0,x1,y1,s1,x2,y2,s2):
        p0 = gauss(x,x0,y0,s0)
        p1 = gauss(x,x1,y1,s1)
        g2 = 2*x1 - x0 +x2
        p2 = gauss(x,g2,y2,s2)
        return p0+p1+p2
    ng = len(P)/3
    mx = Qhist[1]
    my = Qhist[0]
    merr = None
    abSig = None
    if doErr:
        args = Qhist[3]
        mx = mx[args]
        my = my[args]
        merr = np.sqrt(Qhist[2][args])
        abSig = True
    if ng==3:
        fit,tmp = curve_fit(gauss3,mx,my,p0=P,bounds=bounds,sigma=merr,absolute_sigma=abSig,maxfev=50000,ftol=1e-7,gtol=1e-7)
    elif ng==2:
        fit,tmp = curve_fit(gauss2,mx,my,p0=P,bounds=bounds,sigma=merr,absolute_sigma=abSig,maxfev=50000,ftol=1e-7,gtol=1e-7)
    elif ng==1:
        fit,tmp = curve_fit(gauss,mx,my,p0=P,bounds=bounds,sigma=merr,absolute_sigma=abSig,maxfev=50000,ftol=1e-7,gtol=1e-7)
    else:
        print('No valid fit function found')
        return None
    return fit,tmp

def fullAna(fname,nbase=150,winS=5,winF=5,cut=4,pmt=0,qbins=200,ret=False,plot=True):
    Qhist = DRun(fname,nbase,winS,winF,cut,pmt,qbins,ret,plot)
    """
    popt,pcov = fitQ(Qhist,[0,10,1,3,1,1])
    QNoise = min([popt[0],popt[3]])
    QSPE = max([popt[0],popt[3]])
    Qbin0 = np.argmax(Qhist[1]>QNoise)
    Qbin1 = np.argmax(Qhist[1]>QSPE)
    print(popt)
    Qvalley = np.argmax(Qhist[1]>(QSPE/5))
    """
    Qvalley = np.argmax(Qhist[1]>(2.2/5))
    DR = np.sum(Qhist[0][Qvalley:])
    DRErr = np.sqrt(np.sum(Qhist[2][Qvalley:]))
    print([fname,DR,DRErr])
    return DR,DRErr

allRuns = [fullAna(iF) for iF in allF]
allDR = np.array([iRun[0] for iRun in allRuns]).flatten()
allDRerr = np.array([iRun[1] for iRun in allRuns]).flatten()

plt.clf()
plt.errorbar(allT/3600.0,allDR,marker='+',yerr=allDRerr)
plt.xlabel('T [hrs]')
plt.ylabel('Rate [Hz]')
plt.savefig("DarkRate.png")
#plt.show()
import cPickle as pickle
with open(wdir+'/DR.npz','wb') as outfile:
    pickle.dump({'File':allF,'Time':allT,'DR':allDR,'DRerr':allDRerr}, outfile, pickle.HIGHEST_PROTOCOL)