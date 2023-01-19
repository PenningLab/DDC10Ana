import numpy as np
import matplotlib.pyplot as plt
from scipy import integrate
from scipy.optimize import curve_fit
from scipy.stats import norm,poisson,expon,chisquare
import AnaUtils as au
import collections

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

def extRunCh(fname,nbase,winS,winF,cut=4,pmt=1,trigM=100,qbins=1000,qW = 0,ret=False,plot=False):
    waves = au.ReadDDC10_BinWave(fname)
    waves[0],base = au.Subtract_Baseline(waves[0],nBase=nbase)
    #require baseline has no pulse. i.e. integral over baseline less than cut*rms
    bmask = np.absolute(integrate.simps(waves[0][:,pmt,:nbase]))<cut*integrate.simps(np.ones(nbase))*base[1][:,pmt]
    if plot:
        au.plotWaves(waves[0],pmt,1000)
    #TrigPeaks = au.peakHist(waves,chan=0,ret=True,doplot=plot)
    #PromptPeak = au.peakHist(waves,chan=pmt,ret=True,doplot=plot)
    
    #tmask = TrigPeaks[2]>trigM #trigger pulse must have amplitude > trigM
    sumWaves = waves[0][:,pmt,:].sum(0)/float(waves[1]['numEvents'])
    if plot:
        plt.plot(sumWaves)
        plt.xlabel('Sample (10ns)')
        plt.ylabel('V')
        plt.show()
    
    Trigshift = np.argmax(sumWaves)
    wStart = Trigshift-winS
    wFin = Trigshift+winF
    evmask = bmask
    print(Trigshift)

    Qhist = au.winQHist(waves,ch=pmt,init=wStart.astype(int),end=wFin.astype(int),nBins=qbins,binW = qW,evMask=evmask[...,np.newaxis])
    if ret:
        Qinfo = waves[1]
        #Qhist['waves'] = waves[0][:,pmt]
        #Qinfo['trgT'] = TrigPeaks[1]
        #Qhist['evMask'] = evmask
        #Qinfo['baserms'] = base[1][:,pmt]
        Qinfo['dof'] = Qhist['dof']
        return Qhist['qHist'],Qinfo
    else:
        return Qhist['qHist']



def gpn(q,n,q0,q1,s0,s1,u,yq=0,ys=1,Ny=0):
    if n==0:
        return norm.pdf(q,q0,np.abs(s0))*float(poisson.pmf(0,u)) + Ny*expon.pdf(q,yq,ys)
    else:
        sn = s0*s0 + (n*s1*s1)
        gan = norm.pdf(q,q0+n*q1,np.sqrt(sn))*float(poisson.pmf(n,u))
        return gan+gpn(q,n-1,q0,q1,s0,s1,u,yq,ys)

def fitQP(Qhist,P,N=50,doErr=False,dof=0):
    P = collections.OrderedDict(P)
    
    ng = len(P)
    mx = Qhist[1]
    mN = Qhist[0].sum()*(mx[1]-mx[0])
    my = Qhist[0]/mN
    merr = None
    abSig = None
    if doErr:
        args = Qhist[3]
        mx = mx[args]
        my = my[args]
        merr = np.sqrt(Qhist[2][args]/(mN*mN))
        abSig = True

    lambdgpn = lambda q,q0,q1,s0,s1,u,yq,ys: gpn(q,N,q0,q1,s0,s1,u,yq,ys,1)
    if 'yq' not in P and 'ys' not in P:
        lambdgpn = lambda q,q0,q1,s0,s1,u: gpn(q,N,q0,q1,s0,s1,u)
    elif 'yq' not in P:
        P['yq'] = 0
        lambdgpn = lambda q,q0,q1,s0,s1,u,yq,ys: gpn(q,N,q0,q1,s0,s1,u,yq,ys,1)
    elif 'ys' not in P:
        P['ys'] = 1
        lambdgpn = lambda q,q0,q1,s0,s1,u,yq,ys: gpn(q,N,q0,q1,s0,s1,u,yq,ys,1)
    
    fit,tmp = curve_fit(lambdgpn,mx,my,p0=list(P.values()),sigma=merr,absolute_sigma=abSig,maxfev=10000,ftol=1e-8,gtol=1e-8)
    chi2r = np.nonzero(mx>0)[0]
    mchi2 = chisquare(my[chi2r],gpn(mx[chi2r],N,*fit),ddof=dof)
    #print(fit)
    params = P.copy()
    params.update(zip(params,fit))
    paramerr = params.copy()
    paramerr.update(zip(paramerr,np.diag(tmp)))
    params['chi2'] = mchi2
    params['norm'] = mN
    return params,paramerr
