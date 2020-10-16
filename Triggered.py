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
    sumWaves = waves[0][:,pmt,:].sum(0)/float(waves[1]['numEvents'])
    
    plt.plot(sumWaves)
    plt.xlabel('Sample (10ns)')
    plt.ylabel('V')
    plt.show()
    
    Trigshift = np.average(PromptPeak[1]-TrigPeaks[1],weights=tmask*PromptPeak[2])
    wStart = TrigPeaks[1]+Trigshift-winS
    wFin = TrigPeaks[1]+Trigshift+winF
    evmask = bmask*tmask

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

def gauss(x, x0, y0, sigma):
    p = [x0, y0, sigma]
    return p[1]* np.exp(-0.5*((x-p[0])/p[2])**2)
def gauss2(x,x0,y0,s0,x1,y1,s1):
    p0 = gauss(x,x0,y0,s0)
    p1 = gauss(x,x1,y1,s1)
    return p0+p1
def gauss3(x,x0,y0,s0,x1,y1,s1,y2,s2):
    p01 = gauss2(x,x0,y0,s0,x1,y1,s1)
    g2 = 2*(x1 - x0)
    p2 = gauss(x,g2,y2,s2)
    return p01+p2
def gauss4(x,x0,y0,s0,x1,y1,s1,y2,s2,y3,s3):
    p012 = gauss3(x,x0,y0,s0,x1,y1,s1,y2,s2)
    g3 = 3*(x1 -x0)
    p3 = gauss(x,g3,y3,s3)
    return p012+p3
    

def fitQ(Qhist,P,bounds=(-np.inf,np.inf),doErr=False,dof=0):
    
    ng = len(P)
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
    if ng==10:
        fit,tmp = curve_fit(gauss4,mx,my,p0=P,bounds=bounds,sigma=merr,absolute_sigma=abSig,maxfev=5000,ftol=1e-7,gtol=1e-7)
        mchi2 = chisquare(my,gauss4(mx,fit[0],fit[1],fit[2],fit[3],fit[4],fit[5],fit[6],fit[7],fit[8],fit[9]),ddof=dof)
    elif ng==8:
        fit,tmp = curve_fit(gauss3,mx,my,p0=P,bounds=bounds,sigma=merr,absolute_sigma=abSig,maxfev=5000,ftol=1e-7,gtol=1e-7)
        mchi2 = chisquare(my,gauss3(mx,fit[0],fit[1],fit[2],fit[3],fit[4],fit[5],fit[6],fit[7]),ddof=dof)
    elif ng==6:
        fit,tmp = curve_fit(gauss2,mx,my,p0=P,bounds=bounds,sigma=merr,absolute_sigma=abSig,maxfev=5000,ftol=1e-7,gtol=1e-7)
        mchi2 = chisquare(my,gauss2(mx,fit[0],fit[1],fit[2],fit[3],fit[4],fit[5]),ddof=dof)
    else:
        fit,tmp = curve_fit(gauss,mx,my,p0=P,bounds=bounds,sigma=merr,absolute_sigma=abSig,maxfev=5000,ftol=1e-7,gtol=1e-7)
        mchi2 = chisquare(my,gauss(mx,fit[0],fit[1],fit[2]),ddof=dof)
    return fit,tmp,mchi2

def gpn(q,n,q0,q1,s0,s1,u):
    if n==0:
        return norm.pdf(q,q0,np.abs(s0))*float(poisson.pmf(0,u))# + expon.pdf(q,yq,ys)
    else:
        sn = s0*s0 + (n*s1*s1)
        gan = norm.pdf(q,q0+n*q1,np.sqrt(sn))*float(poisson.pmf(n,u))
        return gan+gpn(q,n-1,q0,q1,s0,s1,u)

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
        
    fit,tmp = curve_fit(lambda q,q0,q1,s0,s1,u: gpn(q,N,q0,q1,s0,s1,u),mx,my,p0=list(P.values()),sigma=merr,absolute_sigma=abSig,maxfev=5000,ftol=1e-7,gtol=1e-7)
    mchi2 = chisquare(my,gpn(mx,N,*fit),ddof=dof)
    #print(fit)
    params = P.copy()
    params.update(zip(params,fit))
    paramerr = params.copy()
    paramerr.update(zip(paramerr,np.diag(tmp)))
    params['chi2'] = mchi2
    params['norm'] = mN
    return params,paramerr
