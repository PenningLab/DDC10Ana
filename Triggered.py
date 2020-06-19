import numpy as np
import matplotlib.pyplot as plt
from scipy import integrate
from scipy.optimize import curve_fit
import AnaUtils as au

def extRun(fname,nbase,winS,winF,cut=4,pmt=1,trigM=100,qbins=1000):
    waves = au.ReadDDC10_BinWave(fname)
    waves[0],base = au.Subtract_Baseline(waves[0],nBase=nbase)
    #require baseline has no pulse. i.e. integral over baseline less than cut*rms
    bmask = np.absolute(integrate.simps(waves[0][:,pmt,:nbase]))<cut*integrate.simps(np.ones(nbase))*base[1][:,pmt]

    TrigPeaks = au.peakHist(waves,chan=0,ret=True)[1:]
    wStart = TrigPeaks[0]+winS
    wFin = TrigPeaks[0]+winF
    tmask = TrigPeaks[1]>trigM #trigger pulse must have amplitude > trigM
    evmask = bmask*tmask

    Qhist = au.winQHist(waves,ch=1,init=wStart.astype(int),end=wFin.astype(int),nBins=qbins,evMask=evmask[...,np.newaxis])['qHist']

    return Qhist

def fitQ(Qhist,P,bounds=(-np.inf,np.inf)):
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
        p2 = gauss(x,x2,y2,s2)
        return p0+p1+p2
    ng = len(P)/3
    if ng==3:
        fit,tmp = curve_fit(gauss3,Qhist[1],Qhist[0],p0=P,bounds=bounds)
    if ng==2:
        fit,tmp = curve_fit(gauss2,Qhist[1],Qhist[0],p0=P,bounds=bounds)
    else:
        fit,tmp = curve_fit(gauss,Qhist[1],Qhist[0],p0=P,bounds=bounds)
    return fit,tmp
