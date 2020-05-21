import numpy as np
import uproot
import matplotlib.pyplot as plt
from scipy import integrate
from pylab import rcParams
rcParams['figure.figsize'] = 15, 11

ClockDDC10 = 6e8
adccperVolt = 8192
resistance_ohm = 50
sampleWidth_ns = 10

def ReadDDC10_BinWave(fName, doTime=True):
    waveInfo = {}
    fp = open(fName+'.bin','rb')
    Header = np.fromfile(fp,dtype=np.uint32,count=4)
    
    waveInfo['numEvents'] = int(Header[0])
    waveInfo['numSamples'] = int(Header[1])
    waveInfo['chMap'] = np.array([1 if digit=='1' else 0 for digit in bin(Header[2])[2:]])
    waveInfo['numChan'] = np.sum(waveInfo['chMap'])
    byteOrderPattern = hex(int(np.fromfile(fp,dtype=np.uint32,count=1)))
    
    tmpArr = np.concatenate((np.fromfile(fp,dtype=np.int16),np.zeros(2)))
    fp.close()
    
    waveArr = np.reshape(tmpArr.astype(dtype=np.float64)/adccperVolt,(waveInfo['numEvents'],waveInfo['numChan'],(waveInfo['numSamples']+6)))[...,2:-4]

    if doTime:
        with open(fName+'.log','r') as fl:
            waveInfo['liveTimes_s'] = np.loadtxt(fl,delimiter=',',skiprows=5,max_rows=waveInfo['numEvents'],usecols=(2),dtype=np.float64)/(ClockDDC10)
            waveInfo['totliveTime_s'] = np.sum(waveInfo['liveTimes_s'])
            
    return [waveArr,waveInfo]


def Subtract_Baseline(waveArr,nBase=150):
    baseWave = waveArr[...,:nBase]
    sumax = len(waveArr.shape)-1
    waveBaseline = np.sum(baseWave,axis = sumax)/nBase
    waveBaserms = np.sqrt(np.sum(baseWave*baseWave,axis = sumax)/nBase - waveBaseline*waveBaseline)
    subtwaveArr = waveArr - waveBaseline[...,np.newaxis]
    
    return subtwaveArr,(waveBaseline,waveBaserms)

def BaseCut(waveArr,baseD,cut=6):
    mask = waveArr>6*baseD[1]

def winQHist(wave,ch=0,init=175,end=250,nBins=10000,hrange=[-10,250],sub=True):
    if sub:
        wave[0],baseD = Subtract_Baseline(wave[0])
    sumax = len(wave[0][:,ch,:].shape)-1
    qArr = 1e3*integrate.simps(wave[0][:,ch,init:end],axis=sumax)*sampleWidth_ns/resistance_ohm
    ret = {'qData':qArr}
    ret['qHist'] = np.histogram(qArr,bins=nBins,range=hrange,weights=np.ones(shape=qArr.shape)/wave[1]['totliveTime_s'])
    return ret

import matplotlib as mpl
from pylab import rcParams
rcParams['figure.figsize'] = 15, 11
def peakHist(waveArr,chan=0,yrange=(-1,1),yscale=1):
    peakT = np.argmax(np.absolute(waveArr[0][:,chan,:]),axis=1)
    peakV = np.array([waveArr[0][ev,chan,peakT[ev]] for ev in range(waveArr[1]['numEvents'])])*1e3
    pHist = np.histogram2d(peakT,peakV,bins=[waveArr[1]['numSamples'],int(adccperVolt/yscale)],range=[[0,waveArr[1]['numSamples']],[-1e3,1e3]])
    plt.pcolormesh(pHist[1][:-1],pHist[2][:-1],np.transpose(pHist[0])/waveArr[1]['totliveTime_s'],norm=mpl.colors.LogNorm())
    cbar = plt.colorbar()
    plt.xlabel("peak Time (samples)")
    plt.ylabel("peak Amplitude (mV)")
    plt.ylim(yrange)
    plt.show()
    plt.plot(pHist[1][:-1],np.sum(pHist[0],axis=1))
    plt.xlabel("peak Time (samples)")
    plt.show()
    plt.plot(pHist[2][:-1],np.sum(pHist[0],axis=0))
    plt.xlabel("peak Amplitude (mV)")
    plt.show()
    return pHist

def plotWaves(waveArr,chan=0,nWaves=100):
    plt.figure()
    for i in range(min(nWaves,len(waveArr))):
        plt.plot(waveArr[i,chan,:],marker='+')
    plt.xlabel('samples (10ns)')
    plt.ylabel('V')
    plt.show()