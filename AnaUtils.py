import numpy as np
import uproot
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from scipy.stats import norm,poisson,expon,chisquare
from scipy import integrate
import awkward as awk
from pylab import rcParams
rcParams['figure.figsize'] = 15, 11
import collections

ClockDDC10 = 6e8
adccperVolt = 8192
resistance_ohm = 50
sampleWidth_ns = 10

def ReadDDC10_BinWave(fName, doTime=True):
    waveInfo = {}
    with open(fName+'.bin','rb') as fp:
        try:
            Header = np.fromfile(fp,dtype=np.uint32,count=4)

            waveInfo['numEvents'] = int(Header[0])
            waveInfo['numSamples'] = int(Header[1])
            waveInfo['chMap'] = np.array([1 if digit=='1' else 0 for digit in bin(Header[2])[2:]])
            waveInfo['numChan'] = np.sum(waveInfo['chMap'])
            waveInfo['file'] = fName
            byteOrderPattern = hex(int(np.fromfile(fp,dtype=np.uint32,count=1)))
            print(waveInfo)
        except ValueError as e:
            print(e)
            return None
        
        try:
            waveArr = np.empty((waveInfo['numEvents']*waveInfo['numChan']*(waveInfo['numSamples']+6)),dtype=np.int16)
            waveArr[:-2] = np.fromfile(fp,dtype=np.int16)
        except ValueError as e:
            print(e)
            return None
    
    waveArr = np.reshape(waveArr.astype(dtype=np.float64)/adccperVolt,(waveInfo['numEvents'],waveInfo['numChan'],(waveInfo['numSamples']+6)))[...,2:-4]

    if doTime:
        with open(fName+'.log','r') as fl:
            waveInfo['liveTimes_s'] = np.loadtxt(fl,delimiter=',',skiprows=5,max_rows=waveInfo['numEvents'],usecols=(2),dtype=np.float64)/(ClockDDC10)
            waveInfo['totliveTime_s'] = np.sum(waveInfo['liveTimes_s'])
            
    return [waveArr,waveInfo]


def Subtract_Baseline(waveArr,nBase=50):
    baseWave = waveArr[...,:nBase]
    sumax = len(waveArr.shape)-1
    waveBaseline = np.sum(baseWave,axis = sumax)/nBase
    waveBaserms = np.sqrt(np.sum(baseWave*baseWave,axis = sumax)/nBase - waveBaseline*waveBaseline)
    subtwaveArr = waveArr - waveBaseline[...,np.newaxis]
    
    return subtwaveArr,(waveBaseline,waveBaserms)

from collections.abc import Iterable
def winQHist(wave,ch,init=175,end=250,nBins=10000,hrange=None,sub=False,evMask=True,nBase=50,doLive=True,binW=0):
    if sub:
        wave[0],baseD = Subtract_Baseline(wave[0],nBase)
    sumax = len(wave[0][:,ch,:].shape)-1
    wmask=1
    if isinstance(init,Iterable):
        wmask1 = np.indices(wave[0][:,ch].shape)[1]>init[...,np.newaxis]
        wmask *= wmask1
    else:
        wmask1 = np.indices(wave[0][:,ch].shape)[1]>init
        wmask *= wmask1
    if isinstance(end,Iterable): 
        wmask1 = np.indices(wave[0][:,ch].shape)[1]<end[...,np.newaxis]
        wmask *= wmask1
    else:
        wmask1 = np.indices(wave[0][:,ch].shape)[1]<end
        wmask *= wmask1
    qArr = 1e3*integrate.simps(evMask*wmask*wave[0][:,ch])*sampleWidth_ns/resistance_ohm
    ret = {'qData':qArr}
    ret['dof'] = len(ret['qData'])-1
    if isinstance(hrange,Iterable):
        bRange = hrange[1]-hrange[0]
    else:
        bRange = np.amax(qArr)-np.amin(qArr)
    if binW>0:
        nBins = int(bRange/binW)
    tmpQ = list(np.histogram(qArr,bins=nBins,range=hrange))
    tmpQ[0] = tmpQ[0].astype(float)
    bWidth = (tmpQ[1][-1] - tmpQ[1][0])/float(nBins)
    bTot = tmpQ[0].sum()
    bNorm = bTot*bWidth
    tmpQ[1] = (tmpQ[1][1:]+tmpQ[1][:-1])/2.0
    tmpQ.append(tmpQ[0]*np.square(1.0/bNorm))
    tmpQ[0] *= 1.0/bNorm
    if doLive:
        tmpQ[0] *= bTot/float(wave[1]['totliveTime_s'])
        tmpQ[2] *= np.square(bTot/float(wave[1]['totliveTime_s']))
    tmpQ.append(np.nonzero(tmpQ[2])[0])
    
    ret['qHist'] = list(tmpQ)
    return ret

import matplotlib as mpl
from pylab import rcParams
rcParams['figure.figsize'] = 15, 11
mpl.rc('axes.formatter', useoffset=False)
def peakHist(waveArr,chan=0,yrange=None,yscale=1,ret=False,doplot=True):
    peakT = np.argmax(waveArr[0][:,chan,:],axis=1)
    peakV = waveArr[0][np.arange(0,waveArr[1]['numEvents']),chan,peakT]*1e3
    pHist = np.histogram2d(peakT,peakV,bins=[waveArr[1]['numSamples'],int(adccperVolt/yscale)])
    if doplot:
        plt.pcolormesh(pHist[1][:-1],pHist[2][:-1],np.transpose(pHist[0])/waveArr[1]['totliveTime_s'],norm=mpl.colors.LogNorm())
        cbar = plt.colorbar()
        plt.xlabel("peak Time (samples)")
        plt.ylabel("peak Amplitude (mV)")
        if isinstance(yrange,(tuple,list)):
            plt.ylim(yrange)
        maxT = np.argmax(np.sum(pHist[0],axis=1))
        plt.xlim(pHist[1][maxT]-100,pHist[1][maxT]+100)
        plt.show()
        plt.plot(pHist[1][:-1],np.sum(pHist[0],axis=1))
        plt.yscale('log')
        plt.xlabel("peak Time (samples)")
        plt.show()
        plt.plot(pHist[2][:-1],np.sum(pHist[0],axis=0))
        plt.yscale('log')
        plt.xlabel("peak Amplitude (mV)")
        plt.show()
    if ret:
        return pHist,peakT,peakV
    else:
        return pHist

def plotWaves(waveArr,chan=0,nWaves=100):
    plt.figure()
    for i in range(min(nWaves,len(waveArr))):
        plt.plot(waveArr[i,chan,:],marker='+')
    plt.xlabel('samples (10ns)')
    plt.ylabel('V')
    plt.show()
    return plt.gcf()

def gpn(q,n,q0,q1,s0,s1,u):
    if n==0:
        return norm.pdf(q,q0,np.abs(s0))*float(poisson.pmf(0,u))
    else:
        sn = s0*s0 + (n*s1*s1)
        gan = norm.pdf(q,q0+n*q1,np.sqrt(sn))*float(poisson.pmf(n,u))
        return gan+gpn(q,n-1,q0,q1,s0,s1,u)
def gpn2(q,n,q0,q1,s0,s1,u,Na=1):
    return Na*gpn(q,n,q0,q1,s0,s1,u)

def g2(q,q0,q1,s0,s1):
		return norm.pdf(q,q0,np.abs(s0)) + norm.pdf(q,q1,np.abs(s1))

def fitQ(Qhist,P,doErr=False,dof=0):
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
    lambdg = lambda q,q0,q1,s0,s1: g2(q,q0,q1,s0,s1)
    
    fit,tmp = curve_fit(lambdg,mx,my,p0=list(P.values()),bounds=([-1,0,0,0],np.inf),sigma=merr,absolute_sigma=abSig,maxfev=10000,ftol=1e-8,gtol=1e-8)
    mchi2 = chisquare(my,g2(mx,*fit),ddof=dof)
    #print(fit)
    params = P.copy()
    params.update(zip(params,fit))
    paramerr = params.copy()
    paramerr.update(zip(paramerr,np.diag(tmp)))
    params['chi2'] = mchi2
    params['norm'] = mN
    return params,paramerr

def fitQP(Qhist,P,N=50,doErr=False,dof=0):
    P['Na'] = 1
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
    lambdgpn = lambda q,q0,q1,s0,s1,u,Na: gpn2(q,N,q0,q1,s0,s1,u,Na)
    
    fit,tmp = curve_fit(lambdgpn,mx,my,p0=list(P.values()),bounds=([-np.inf,-np.inf,0,0,0,-np.inf],np.inf),sigma=merr,absolute_sigma=abSig,maxfev=10000,ftol=1e-8,gtol=1e-8)
    mchi2 = chisquare(my,gpn2(mx,N,*fit),ddof=dof)
    #print(fit)
    params = P.copy()
    params.update(zip(params,fit))
    paramerr = params.copy()
    paramerr.update(zip(paramerr,np.diag(tmp)))
    params['chi2'] = mchi2
    params['norm'] = mN
    return params,paramerr
    
#Pulse Finding 
#create kernel for Laplacian of Gaussian edge finding filter
def LoGkernel(sigma=1,scale=5,norm=False):
    sigma2 = sigma*sigma
    size = int(scale*sigma)
    x = (np.arange(size) - (size-1)/2.0)
    kernel = (x*x/sigma2 - 1)/sigma2
    N = 1/(sigma*np.sqrt(2*np.pi)) if norm else 1
    x2 = N*np.exp(-x*x/(2*sigma2))
    print(x)
    print(kernel)
    LoG = kernel*x2
    return LoG

#filter performed using signal.fftconvolve() which allows convolution along axis, i.e. all filtered wave forms generated in one line utiizing scipy optimizations

#calculate minimum maximum and gradient of sliding window of Filtered Waveform
def mmg_rolling(a, window):
    axis =-1
    shape = a.shape[:axis] + (a.shape[axis] - window + 1, window)
    strides = a.strides + (a.strides[-1],)
    rolling = np.lib.stride_tricks.as_strided(a, shape=shape, strides=strides)
    grad = (rolling[...,-1]-rolling[...,0])/float(window-1.0)
    return np.max(rolling,axis=axis),np.min(rolling,axis=axis),grad

#zero crossings of filtered waveform are candidate edges, gradient discriminates rising edge v falling edge candidates
#-1 elements are rising edge candidates and +1 elements are falling edge candidates
def zero_crossing(mLoG,mWave,thresh,window=3):
    maxL,minL,gradL = mmg_rolling(mLoG,window)
    #print(grad2L,gradL)
    pzero = (mLoG[...,1:-1]>0)
    zeroCross = np.zeros(shape=mLoG.shape).astype(np.int)
    zeroCross[...,1:-1] = pzero*(minL<0) + (1-pzero)*(maxL>0)
    diffL = maxL-minL
    zeroCross[...,1:-1] = zeroCross[...,1:-1]*(diffL > thresh)
    zeroCross[...,1:-1] = ((gradL>thresh/window).astype(np.int)-(gradL<-thresh/window).astype(np.int))*zeroCross[...,1:-1]
    
    #zeroCross[zeroCross==0] = np.nan
    lEd = awk.fromiter([np.nonzero(zi)[0] for zi in zeroCross<0])
    rEd = awk.fromiter([np.nonzero(zi)[0] for zi in zeroCross>0])
    return (lEd,rEd),np.pad(gradL,[(0,)]*(gradL.ndim-1)+[(1,)],'constant',constant_values=(0))


#Now sort through candidate edges
#Maybe integrate from every left edge to every right edge (one sided search)
#Find minima between left right edge pairs
#Set a min threshold for the integral
#All edge combinations with integral above threshold kept
#Sort kept ranges by size
#Starting from smallest range integrate beyond boundary (both sides) until fractional change in integral is <1% (set as new left and right edge)
#Now look for overlapping pulses and merge