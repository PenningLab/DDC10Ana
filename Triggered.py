import numpy as np
import matplotlib.pyplot as plt
from scipy import integrate
from scipy.optimize import curve_fit
import AnaUtils as au

def extRun(fname,nbase,winS,winF,cut=4,pmt=1,trigM=100,qbins=1000,qW = 0,ret=False,plot=False):
    waves = au.ReadDDC10_BinWave(fname)
    waves[0],base = au.Subtract_Baseline(waves[0],nBase=nbase)
    #require baseline has no pulse. i.e. integral over baseline less than cut*rms
    bmask = np.absolute(integrate.simps(waves[0][:,pmt,:nbase]))<cut*integrate.simps(np.ones(nbase))*base[1][:,pmt]
    if plot:
        au.plotWaves(waves[0],pmt,1000)
    TrigPeaks = au.peakHist(waves,chan=0,ret=True,doplot=plot)
    PromptPeak = au.peakHist(waves,chan=pmt,ret=True,doplot=plot)
    
    tmask = TrigPeaks[2]>trigM #trigger pulse must have amplitude > trigM
    
    Trigshift = np.average(PromptPeak[1]-TrigPeaks[1],weights=tmask*PromptPeak[2])
    wStart = TrigPeaks[1]+Trigshift-winS
    wFin = TrigPeaks[1]+Trigshift+winF
    evmask = bmask*tmask

    Qhist = au.winQHist(waves,ch=pmt,init=wStart.astype(int),end=wFin.astype(int),nBins=qbins,binW = qW,evMask=evmask[...,np.newaxis])
    if ret:
        Qhist['waves'] = waves[0][:,pmt]
        Qhist['trgT'] = TrigPeaks[1]
        Qhist['evMask'] = evmask
        Qhist['baserms'] = base[1][:,pmt]
        return Qhist,waves[1]
    else:
        return Qhist['qHist']

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
        fit,tmp = curve_fit(gauss3,mx,my,p0=P,bounds=bounds,sigma=merr,absolute_sigma=abSig,maxfev=5000,ftol=1e-7,gtol=1e-7)
    if ng==2:
        fit,tmp = curve_fit(gauss2,mx,my,p0=P,bounds=bounds,sigma=merr,absolute_sigma=abSig,maxfev=5000,ftol=1e-7,gtol=1e-7)
    else:
        fit,tmp = curve_fit(gauss,mx,my,p0=P,bounds=bounds,sigma=merr,absolute_sigma=abSig,maxfev=5000,ftol=1e-7,gtol=1e-7)
    return fit,tmp
