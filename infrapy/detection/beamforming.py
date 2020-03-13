from obspy.core import UTCDateTime
from obspy.core import trace
from obspy.core import Stream
from obspy.core.util import AttribDict
from obspy.signal.array_analysis import *
import obspy.signal.invsim as inv
#import obspy.signal.util.utlGeoKm

import matplotlib.pyplot as pl
import matplotlib.mlab as mlab

from obspy.signal.util import util_geo_km, next_pow_2
import numpy as np
from ..utils.cart2pol import cart2pol
from ..utils.obspy_conversion import stream2array

import cmath
import math
import itertools

from scipy import interpolate
from scipy import signal
from numpy import linalg as LA
import scipy.linalg as scli
import sys

from scipy.ndimage.filters import maximum_filter
from scipy.ndimage.morphology import generate_binary_structure, binary_erosion
from sklearn import mixture

import math

import warnings
warnings.filterwarnings("ignore")
warnings.filterwarnings('default', category=PendingDeprecationWarning)

#Modified   AIC   and   MDL   Model Selection  Criteria  for  short  Data  Records
def AIC(eigenV,n,m):
    valT=np.zeros(m)
    for k in range(m):
        ll=eigenV[k+1:]
        a1=np.prod(ll)
        a2=np.sum(ll)/(m-k)
        Ln=np.power(a1,1.0/(m-k))/a2
        valT[k]=-2*np.log(Ln)+2*k*(2*m-k)
        print(valT[k])
    return valT


def MDL(eigenV,n,m):
    valT=np.zeros(m)
    for k in range(m):
        ll=eigenV[k+1:]
        a1=np.prod(ll)
        a2=np.sum(ll)/(m-k)
        Ln=np.power(a1,1.0/(m-k))/a2
        valT[k]=-(m-k)*np.log(Ln)+k*(2*m-k)/2
        print(valT[k])
    return valT



def detect_peaks(image,size=3,num_peaks=None):
    """
    Takes an image and detect the peaks using the local maximum filter.
    Returns a boolean mask of the peaks (i.e. 1 when
    the pixel's value is the neighborhood maximum, 0 otherwise)
    """
    if num_peaks==None:
        num_peaks=1

    local_max = maximum_filter(image,size=size)
    #local_max is a mask that contains the peaks we are
    #looking for, but also the background.
    #In order to isolate the peaks we must remove the background from the mask.

    #we create the mask of the background
    local_maxima = (image == local_max) * ( image!= np.zeros(image.shape))
    POSM=np.argwhere(local_maxima==True)


    [nn,mm]=POSM.shape
    val=[]
    pos=[]
    for pp in range(nn):
        val.append(image[POSM[pp,0],POSM[pp,1]])
        pos.append([POSM[pp,0],POSM[pp,1]])

    val=np.asarray(val)
    pos=np.asarray(pos)
    valS=val[val.argsort()[::-1]]
    posS=pos[val.argsort()[::-1]]
    valF=np.zeros((num_peaks))*np.nan
    posF=np.zeros((num_peaks,2))*np.nan
    try:
        if nn >= num_peaks:
            valF=valS[0:num_peaks]
            posF=posS[0:num_peaks]
        else:
            if nn>0:
                valF[0:nn]=valS[0:nn]
                posF[0:nn]=posS[0:nn]
    except Exception as ex1:
        print('problem with the num of peaks')
        embed()
        exit()
    return posF,valF

def procPEAKS(peaks,slow):
    try:
        if (peaks[1]>len(slow) or (peaks[1]<0) or (peaks[0]>len(slow)) or (peaks[0]<0)) or (np.isnan(peaks[0]))or (np.isnan(peaks[1])):
            return np.nan,np.nan,np.nan,np.nan
        bz, slofk = cart2pol(slow[int(peaks[1])],slow[int(peaks[0])])
        bz = bz*180.0/math.pi
        if (bz < 0.0):
            bz = bz + 360
        slowx=slow[int(peaks[0])]
        slowy=slow[int(peaks[1])]
        return bz,slofk,slowx,slowy
    except Exception as ex1:
        print('problem inside procPEAKS',ex1)
        embed()
        exit()


def _fkBARTLET(stream,sps,slow,mult_vectors,fN,x,y,timeSTAMP,func,freqN):
    ss_ch,ss_len=stream.shape
    win_lf=float(ss_len/sps)
    win_l=int(win_lf)
    nchan=ss_ch
    #fN_S,freqN_S,mult_vectors=_mult_vectors(win_lf,sps,slow,x,y,nchan,fmin=freqN[0],fmax=freqN[-1])
    try:
        streamF=np.fft.rfft(stream,axis=1)
        streamF=np.conj(streamF)
        streamF = streamF[:,fN.astype('int')]
        tot_mul=np.multiply(mult_vectors,streamF)
        inte1=np.sum(tot_mul,axis=2)
        inte2=np.absolute(inte1)
        inte3=np.power(inte2,2)
        FK=np.sum(inte3,axis=2)
    except Exception as ex1:
        print(ex1)
        embed()
        exit()
        FK=[]
    return FK

def _fkCAPON(stream,sps,slow,mult_vectors,fN,x,y,timeSTAMP,func,freqN):
    # remember it is 'lonlat'
    streamF=np.fft.rfft(stream,axis=1)
    streamF = streamF[:,fN.astype('int')]

    nn,mm=streamF.shape
    xx,yy,s1,s2=mult_vectors.shape
    C=np.zeros((xx,yy,mm),dtype=np.complex)
    for f_i in range(mm):
        SS=np.matrix(streamF[:,f_i])
        covM=np.dot(np.conj(SS.T),SS)
        DAMP=0.0001*np.matrix(np.identity(nn))
        icovM=np.linalg.pinv(covM+DAMP)
        mvec=mult_vectors[:,:,:,f_i]
        for x1 in range(xx):
            for y1 in range(yy):
                sv=np.matrix(mvec[x1,y1,:])
                CC=np.dot(np.dot(np.conj(sv),icovM),sv.T)
                C[x1,y1,f_i]=np.ravel(CC)[0]
    VAL=1/np.sum(np.abs(C),axis=2)
    # my mult_vectors are the complemetnary direction of z(theta), we need to inver the direction of the FK plane
    C_R=VAL[::-1,::-1]
    return C_R



def _mult_vectors(wind_len,sps,slow,x,y,nchan,fmin,fmax):
    '''
    these are the steering vectors

    to-do: we can improve this multiplication
    output
    fN are the index of the frquncy we are interested
    freqsN are the actual frequencies we ar einterested

    '''
    band = [0, int(sps/2.0)]
    dfreq = 1./float(wind_len)
    lfreq = float(math.ceil(min(band)/dfreq) + 1)
    hfreq = float(math.floor(max(band)/dfreq) + 1)
    f = list(range(int(lfreq),int(hfreq)+1))
    f=np.asarray(f,dtype='float')
    f_1 = [fi-1 for fi in f]
    freqs=np.linspace(0, sps/2.0,len(f))
    fN=f[(freqs>=fmin)&(freqs<=fmax)]
    freqsN=freqs[(freqs>=fmin)&(freqs<=fmax)]
    #pdx = np.dot(np.matrix(x).transpose(),np.matrix(-slow))
    #pdy = np.dot(np.matrix(y).transpose(),np.matrix(-slow))
    #wl=np.matrix(fN)*dfreq
    pdx = np.zeros(shape=(len(x),len(slow)))
    for aa in range(len(x)):
        for bb in range(len(slow)):
            pdx[aa][bb]=x[aa]*-slow[bb]
    pdy = np.zeros(shape=(len(y),len(slow)))
    for aa in range(len(y)):
        for bb in range(len(slow)):
            pdy[aa][bb]=y[aa]*-slow[bb]
    wl=fN*dfreq
    epdxf = np.exp(complex(0,-1) * 2 * math.pi *  np.outer(pdx.transpose().ravel(),wl.transpose()))
    epdyf = np.exp(complex(0,-1) * 2 * math.pi *  np.outer(pdy.transpose().ravel(),wl.transpose()))
    numslow = float(len(slow))
    [MS,NS]=np.shape(epdxf)
    MS=float(MS)
    NS=float(NS)
    nslow=len(slow)
    chans_listINT=np.arange(0,int(nchan),dtype='int')
    mult_vectors=np.zeros((int(nslow), int(nslow),int(nchan),int(NS)),dtype=complex)
    for sx in range(0,nslow):
        sxx = int(sx)*int(nchan) + chans_listINT
        for sy in range(0,int(nslow)):
            syy = int(sy)*int(nchan) + chans_listINT
            mult_vectors[int(sx),int(sy),:,:] = np.multiply(epdxf[sxx,:],epdyf[syy,:])

    return fN,freqsN,mult_vectors

def _array_response_short(sx,sy,mult_vectors,wind_len,sps,slow,x,y,nchan,fmin,fmax):
    band = [0, int(sps/2.0)]
    dfreq = 1./float(wind_len)
    lfreq = float(math.ceil(min(band)/dfreq) + 1)
    hfreq = float(math.floor(max(band)/dfreq) + 1)
    f = list(range(int(lfreq),int(hfreq)+1))
    f=np.asarray(f,dtype='float')
    f_1 = [fi-1 for fi in f]
    freqs=np.linspace(0, sps/2.0,len(f))
    fN=f[(freqs>=fmin)&(freqs<=fmax)]
    freqsN=freqs[(freqs>=fmin)&(freqs<=fmax)]
    wl=np.matrix(fN)*dfreq
    pdx_sx = np.dot(np.matrix(x).transpose(),np.matrix(-sx))
    pdy_sy = np.dot(np.matrix(y).transpose(),np.matrix(-sy))
    epdxf_sx = np.exp(complex(0,-1) * 2 * math.pi *  np.outer(pdx_sx.transpose().ravel(),wl.transpose()))
    epdyf_sy = np.exp(complex(0,-1) * 2 * math.pi *  np.outer(pdy_sy.transpose().ravel(),wl.transpose()))
    sx_sy=np.multiply(epdxf_sx,epdyf_sy)
    respUAX=np.multiply(mult_vectors,sx_sy)
    inte1=np.sum(respUAX,axis=2)
    inte2=np.absolute(inte1)
    inte3=np.power(inte2,2)
    FK=np.sum(inte3,axis=2)
    return FK

def _array_response(sx,sy,wind_len,sps,slow,x,y,nchan,fmin,fmax):
    band = [0, int(sps/2.0)]
    dfreq = 1./float(wind_len)
    lfreq = float(math.ceil(min(band)/dfreq) + 1)
    hfreq = float(math.floor(max(band)/dfreq) + 1)
    f = list(range(int(lfreq),int(hfreq)+1))
    f=np.asarray(f,dtype='float')
    f_1 = [fi-1 for fi in f]
    freqs=np.linspace(0, sps/2.0,len(f))
    fN=f[(freqs>=fmin)&(freqs<=fmax)]
    freqsN=freqs[(freqs>=fmin)&(freqs<=fmax)]
    pdx = np.zeros(shape=(len(x),len(slow)))
    for aa in range(len(x)):
        for bb in range(len(slow)):
            pdx[aa][bb]=x[aa]*-slow[bb]
    pdy = np.zeros(shape=(len(y),len(slow)))
    for aa in range(len(y)):
        for bb in range(len(slow)):
            pdy[aa][bb]=y[aa]*-slow[bb]
    #pdx = np.dot(np.matrix(x).transpose(),np.matrix(-slow))
    #pdy = np.dot(np.matrix(y).transpose(),np.matrix(-slow))
    #wl=np.matrix(fN)*dfreq
    wl=fN*dfreq
    epdxf = np.exp(complex(0,-1) * 2 * math.pi *  np.outer(pdx.transpose().ravel(),wl.transpose()))
    epdyf = np.exp(complex(0,-1) * 2 * math.pi *  np.outer(pdy.transpose().ravel(),wl.transpose()))

    pdx_sx = np.dot(np.matrix(x).transpose(),np.matrix(-sx))
    pdy_sy = np.dot(np.matrix(y).transpose(),np.matrix(-sy))
    epdxf_sx = np.exp(complex(0,-1) * 2 * math.pi *  np.outer(pdx_sx.transpose().ravel(),wl.transpose()))
    epdyf_sy = np.exp(complex(0,-1) * 2 * math.pi *  np.outer(pdy_sy.transpose().ravel(),wl.transpose()))
    sx_sy=np.multiply(epdxf_sx,epdxf_sy)
    numslow = float(len(slow))
    [MS,NS]=np.shape(epdxf)
    MS=float(MS)
    NS=float(NS)
    nslow=len(slow)
    chans_listINT=np.arange(0,int(nchan),dtype='int')

    mult_vectors=np.zeros((int(nslow), int(nslow),int(nchan),int(NS)),dtype=complex)
    for sx in range(0,nslow):
        sxx = int(sx)*int(nchan) + chans_listINT
        for sy in range(0,int(nslow)):
            syy = int(sy)*int(nchan) + chans_listINT
            mult_vectors[int(sx),int(sy),:,:] = np.multiply(epdxf[sxx,:],epdyf[syy,:])
    respUAX=np.multiply(mult_vectors,sx_sy)
    inte1=np.sum(respUAX,axis=2)
    inte2=np.absolute(inte1)
    inte3=np.power(inte2,2)
    FK=np.sum(inte3,axis=2)
    return FK

def _fkMUSIC(streamOR,sps,slow,mult_vectors,fN,x,y,timeSTAMP,func,freqN,number_sources,number_div=None):
    # remember it is 'lonlat'
    # streamOR array data
    # sps sampling rate (Hz)
    # mult_vectors this help improving the multiplicaiton especailly when using multiprocessing
    # fN
    stream=streamOR.copy()
    ss_ch,ss_len=stream.shape
    num_win=float(number_div)
    win_lf=float((ss_len/sps)/num_win)
    win_l=int(win_lf)
    nchan=ss_ch
    fN_S,freqN_S,mult_vectors_S=_mult_vectors(win_lf,sps,slow,x,y,nchan,fmin=freqN[0],fmax=freqN[-1])
    #covT=[]
    xx,yy,s1,mm=mult_vectors_S.shape

    covT=[]
    for w_i in range(int(num_win)):
        streamS=stream[:,int(sps)*win_l*w_i:int(sps*win_lf*(w_i+1))]
        nn_s,mm_s=streamS.shape
        val_TA=signal.tukey(mm_s, alpha=0.25, sym=True)
        for ii in range(nn_s):
            streamS[ii,:]=streamS[ii,:]*val_TA
        streamF=np.fft.rfft(streamS,axis=1)
        streamF = streamF[:,fN_S.astype('int')]
        nn,mm=streamF.shape
        cov_W=[]

        for f_i in range(mm):
            SS=np.matrix(streamF[:,f_i]).T
            covM=np.dot(SS,np.conj(SS).T)
            #covM=(1-damp)*covM-damp*np.eye(nn)
            cov_W.append(covM)
        covT.append(np.asarray(cov_W))
    covT=np.asarray(covT)
    covT_AV=np.mean(covT,axis=0)
    C=np.zeros((xx,yy,mm))
    C2=np.zeros((xx,yy,mm))
    for f_i in range(mm):
        dat=covT_AV[f_i]
        temp = dat + (1.0e-3 )*np.trace(dat)/nchan
        w,v=scli.eig(temp)
        argwS=np.argsort(w,kind='mergesort')
        argwS=argwS[::-1]
        wS=w[argwS]
        vS=v[:,argwS]
        vSN=vS[:,number_sources:]
        wSN=wS[number_sources:]
        mvec=mult_vectors_S[:,:,:,f_i]
        ss=np.conj(mvec)
        ss2=np.dot(np.asarray(ss),np.asarray(vSN))
        ss3=np.power(np.abs(ss2),2)
        ss4=np.sum(ss3,axis=2)
        C[:,:,f_i]=1.0/ss4

    try:
        VAL=np.mean(C,axis=2)
    except Exception as ex1:
        print('problem with sum')
        VAL=0*np.sum(np.abs(C),axis=2)
    FK=VAL
    return FK


def xcorr_beam(dat_st):
    sig1=dat_st[0]-np.mean(dat_st[0])
    sig1N=sig1/np.std(sig1)
    ss=[]
    coXC=[]
    for ii in range(dat_st.shape[0]):
        sig2=dat_st[ii]-np.mean(dat_st[ii])
        sig2N=sig2/np.std(sig2)
        tt=np.correlate(sig1N,sig2N,'full')/float(len(sig1))
        sigN=np.roll(dat_st[ii],-(len(dat_st[ii])-np.argmax(tt)))
        ss.append(sigN)
        print(np.max(tt))
        coXC.append(np.max(tt))
    coXC=np.asarray(coXC)
    ind=np.argsort(coXC)
    ind=ind[::-1]
    ss=np.asarray(ss)
    ssSORT=ss[ind]
    return ss,ssSORT,coXC

def _fkMUSIC_dynamic(stream,sps,slow,mult_vectors,fN,x,y,timeSTAMP,func,freqN,number_sources,number_div=None):
    # remember it is 'lonlat'
    overlap=0.5
    ss_ch,ss_len=stream.shape
    # one is added to accomodate for the overlapping
    num_win=float(number_div)
    step_l=(ss_len/sps)/(num_win+1)
    win_lf=float(2*(ss_len/sps)/(num_win+1))
    win_l=int(win_lf)
    nchan=ss_ch
    fN_S,freqN_S,mult_vectors_S=_mult_vectors(win_lf,sps,slow,x,y,nchan,fmin=freqN[0],fmax=freqN[-1])
    xx,yy,s1,mm=mult_vectors_S.shape
    covT=[]


    pl.figure()
    ax1=pl.subplot(1,2,1)
    ax2=pl.subplot(1,2,2)

    for w_i in range(int(num_win)):
        st_p=int(step_l*w_i*sps)
        en_p=int(st_p+win_lf*sps)
        streamSA=stream[:,st_p:en_p].copy()
        nn_s,mm_s=streamSA.shape
        val_TA=signal.hann(mm_s,sym=True)
        for ii in range(nn_s):
            streamSA[ii,:]=streamSA[ii,:]-np.mean(streamSA[ii,:])
            streamSA[ii,:]=streamSA[ii,:]*val_TA
        streamFO=np.fft.rfft(streamSA,axis=1)
        streamF = streamFO[:,fN_S.astype('int')]
        nn,mm=streamF.shape
        cov_W=[]
        for f_i in range(mm):
            SS=np.matrix(streamF[:,f_i])
            covM=np.dot(SS.T,np.conj(SS))
            covM=covM+np.eye(nn_s)*np.trace(covM)/(1E3*nn_s)
            cov_W.append(covM)
        covT.append(np.asarray(cov_W))

    covT=np.asarray(covT)
    covT_AV=np.mean(covT,axis=0)


    C=np.zeros((xx,yy,mm),dtype=np.complex)
    covT_UNI=np.mean(covT_AV,axis=0)
    wO,v=scli.eigh(covT_UNI)
    w=np.sort(wO)[::-1]

    # this section gets the number of sources
    sumV=np.sum(np.abs(w))
    cusum_PER=100*np.cumsum(np.abs(w))/sumV
    n_sources=len(cusum_PER[cusum_PER<85])
    print('at the num of source',n_sources)

    number_sources=n_sources
    #print(np.abs(w[0])/np.abs(w[-1]))
    #AIC(w,win_lf*sps,10)

    if(w[0]/w[-1]>1E2):
        fw=[]
        for f_i in range(mm):
            w,v=scli.eigh(covT_AV[f_i])
            argwS=np.argsort(w,kind='mergesort')
            argwS=argwS[::-1]
            wS=w[argwS]
            vS=v[:,argwS]

            sumV=np.sum(np.abs(wS))
            cusum_PER=100*np.cumsum(np.abs(wS))/sumV
            number_sources=len(cusum_PER[cusum_PER<85])
            print('at the num of source',n_sources,'int:',round(freqN_S[f_i],4),'_',wS[0]/wS[-1])
            vSN=vS[:,number_sources:]
            wSN=wS[number_sources:]
            mvec=mult_vectors_S[:,:,:,f_i]
            ss=np.conj(mvec)
            ss2=np.dot(np.asarray(ss),np.asarray(vSN))
            ss3=np.power(np.abs(ss2),2)
            ss4=np.sum(ss3,axis=2)
            C[:,:,f_i]=1.0/ss4
            posM,valM=detect_peaks(1.0/ss4,size=5,num_peaks=number_sources)

            for ii in range(number_sources):
                RES=procPEAKS(posM[ii],slow)
                print(round(RES[0],1))
                fw.append(round(RES[0],1))
            continue



            vSN=vS[:,number_sources:]
            wSN=wS[number_sources:]
            mvec=mult_vectors_S[:,:,:,f_i]
            ss=np.conj(mvec)
            ss2=np.dot(np.asarray(ss),np.asarray(vSN))
            ss3=np.power(np.abs(ss2),2)
            ss4=np.sum(ss3,axis=2)
            C[:,:,f_i]=1.0/ss4
        try:
            VAL=np.mean(np.abs(C),axis=2)
        except Exception as ex1:
            print('problem with sum')
            VAL=0*np.sum(np.abs(C),axis=2)

        ax2.pcolor(np.log10(VAL))
        vgrad = np.gradient(VAL)
        fulgrad = np.sqrt(vgrad[0]**2 + vgrad[1]**2)
        ax1.pcolor(fulgrad, vmin = np.amin(fulgrad),vmax = np.amax(fulgrad))
        pl.title(n_sources)
        pl.savefig(str(timeSTAMP)+'.png')
        pl.close()
            #pdb.set_trace()
    else:
        print('noise')
        VAL=0*np.sum(np.abs(C),axis=2)
    FK=VAL


    return FK,number_sources







def svdAV_wv(stream,sps,slow,mult_vectors,fN,x,y,timeSTAMP,func,freqN,number_div=None):
    # remember it is 'lonlat'
    ss_ch,ss_len=stream.shape
    num_win=float(number_div)
    win_lf=float((ss_len/sps)/num_win)
    win_l=int(win_lf)
    nchan=ss_ch
    covT=[]
    fN_S,freqN_S,mult_vectors_S=_mult_vectors(win_lf,sps,slow,x,y,nchan,fmin=freqN[0],fmax=freqN[-1])
    for w_i in range(int(num_win)):
        streamS=stream[:,int(sps)*win_l*w_i:int(sps*win_lf*(w_i+1))]
        nn_s,mm_s=streamS.shape
        val_TA=signal.tukey(mm_s, alpha=0.1, sym=True)
        for ii in range(nn_s):
            streamS[ii,:]=streamS[ii,:]*val_TA
        streamF=np.fft.rfft(streamS,axis=1)
        streamF = streamF[:,fN_S.astype('int')]
        nn,mm=streamF.shape
        cov_W=[]
        for f_i in range(mm):
            SS=np.matrix(streamF[:,f_i])
            covM=np.dot(np.conj(SS.T),SS)
            cov_W.append(covM)
        covT.append(np.asarray(cov_W))
    covT=np.asarray(covT)
    covT_AV=np.mean(covT,axis=0)
    aver_CO=np.mean(covT_AV,axis=0)
    w,v=scli.eigh(aver_CO)
    argwS=np.argsort(w,kind='mergesort')
    argwS=argwS[::-1]
    wS_T=w[argwS]
    vS_T=v[:,argwS]
    return wS_T, vS_T

def _fkCAPON_AV(stream,sps,slow,mult_vectors,fN,x,y,timeSTAMP,func,freqN,number_div=None):
    # remember it is 'lonlat'

    ss_ch,ss_len=stream.shape
    num_win=float(number_div)
    win_lf=float((ss_len/sps)/num_win)
    win_l=int(win_lf)
    nchan=ss_ch
    fN_S,freqN_S,mult_vectors_S=_mult_vectors(win_lf,sps,slow,x,y,nchan,fmin=freqN[0],fmax=freqN[-1])
    #covT=[]
    xx,yy,s1,mm=mult_vectors_S.shape

    covT=[]
    for w_i in range(int(num_win)):
        streamS=stream[:,int(sps)*win_l*w_i:int(sps*win_lf*(w_i+1))]
        nn_s,mm_s=streamS.shape
        val_TA=signal.tukey(mm_s, alpha=0.1, sym=True)
        for ii in range(nn_s):
            streamS[ii,:]=streamS[ii,:]*val_TA
        streamF=np.fft.rfft(streamS,axis=1)
        streamF = streamF[:,fN_S.astype('int')]
        nn,mm=streamF.shape
        cov_W=[]
        for f_i in range(mm):
            SS=np.matrix(streamF[:,f_i])
            covM=np.dot(np.conj(SS.T),SS)
            cov_W.append(covM)
        covT.append(np.asarray(cov_W))
    covT=np.asarray(covT)
    covT_AV=np.mean(covT,axis=0)
    C=np.zeros((xx,yy,mm),dtype=np.complex)
    for f_i in range(mm):
        w,v=scli.eigh(covT_AV[f_i])
        argwS=np.argsort(w,kind='mergesort')
        argwS=argwS[::-1]
        wS=w[argwS]
        vS=v[:,argwS]
        vSN=vS[:,:]
        wSN=wS[:]
        mvec=mult_vectors_S[:,:,:,f_i]
        ss=np.conj(mvec)
        ss2=np.dot(np.asarray(ss),np.asarray(vSN))
        ss3=np.power(np.abs(ss2),2)/wSN
        #ss3=np.power(np.abs(ss2),2)
        ss4=np.sum(ss3,axis=2)
        C[:,:,f_i]=1.0/ss4
    try:
        VAL=np.sum(np.abs(C),axis=2)
    except Exception as ex1:
        print('problem with sum')
        VAL=0*np.sum(np.abs(C),axis=2)
        #pdb.set_trace()
    FK=VAL[::-1,::-1]
    return FK

def fkfromOStream(St,wlen,overlap,freqmin,freqmax,slow=None):
    if (slow==None) or (len(slow)==0):
        slow=np.arange(-3.6,3.6,0.05)
    X,Y=getXY_array(St)
    lat2MeanKM=np.asarray(Y,dtype='float')
    lon2MeanKM=np.asarray(X,dtype='float')
    num_p=St[0].stats.npts
    sps=St[0].stats.sampling_rate
    win_l=wlen
    x=lon2MeanKM
    y=lat2MeanKM
    nchan = float(len(lat2MeanKM))
    fN,freqsN,mult_vectors=_mult_vectors(win_l,sps,slow,x,y,nchan,fmin=freqmin,fmax=freqmax)
    curr_pos=0
    data=stream2array(St)
    [m,n]=data.shape
    timeSTAMP=St[0].stats.starttime
    resT=[]
    tt=0
    while (curr_pos+int(win_l*sps))<=num_p:
        dataAUX=data[:,curr_pos:curr_pos+int(win_l*sps)]
        res=fkBARLET(dataAUX,sps,slow,mult_vectors,fN,x,y,timeSTAMP,None,freqsN,None)
        resT.append([res,tt])
        curr_pos=curr_pos+int((win_l-overlap)*sps)
        tt=tt+(win_l-overlap)
    return resT

def tapering(st_aux):
    [ch,samp]=st_aux.shape
    try:
        cosVal=inv.cosTaper(int(samp),p=0.1)
    except:
        cosVal=inv.cosine_taper(int(samp),p=0.1)
    for ii in range(ch):
        st_aux[ii,:]=st_aux[ii,:]*cosVal
    return st_aux


def fkPROC(method,stream,sps,slow,mult_vectors,fN,x,y,timeSTAMP,func,freqN,aux,num_sources=None):
    import warnings
    warnings.filterwarnings("ignore")
    if num_sources==None:
        num_sources=1
    try:
        STO=stream.copy()

        if method=='bartlett':
            ST=tapering(STO.copy())
            FK=_fkBARTLET(ST,sps,slow,mult_vectors,fN,x,y,timeSTAMP,func,freqN)
        elif method=='music':
            ST=STO
            FK=_fkMUSIC(ST,sps,slow,mult_vectors,fN,x,y,timeSTAMP,func,freqN,num_sources,number_div=16)
        elif method=='music16':
            ST=STO
            #print 'inside fkPROC, music16'
            FK=_fkMUSIC(ST,sps,slow,mult_vectors,fN,x,y,timeSTAMP,func,freqN,num_sources,number_div=16)
        elif method=='music32':
            ST=STO
            #print 'inside fkPROC, music32'
            FK=_fkMUSIC(ST,sps,slow,mult_vectors,fN,x,y,timeSTAMP,func,freqN,num_sources,number_div=32)
        elif method=='music32D':
            ST=STO
            #print 'inside fkPROC, music32 dynamic'
            #embed()
            FK,num_sources=_fkMUSIC_dynamic(ST,sps,slow,mult_vectors,fN,x,y,timeSTAMP,func,freqN,num_sources,number_div=32)
        elif method=='music16D':
            ST=STO
            #print 'inside fkPROC, music32 dynamic'
            #embed()
            FK,num_sources=_fkMUSIC_dynamic(ST,sps,slow,mult_vectors,fN,x,y,timeSTAMP,func,freqN,num_sources,number_div=16)
        elif method=='capon16':
            ST=STO
            #print 'inside fkPROC, music32'
            FK=_fkCAPON_AV(ST,sps,slow,mult_vectors,fN,x,y,timeSTAMP,func,freqN,num_sources,number_div=16)
        elif method=='capon32':
            ST=STO
            #print 'inside fkPROC, music32'
            FK=_fkCAPON_AV(ST,sps,slow,mult_vectors,fN,x,y,timeSTAMP,func,freqN,num_sources,number_div=32)
        else:
            print('problem with method selection')
            exit()
        #bz,slofk,slowx,slowy=procFK(FK,slow)
    except Exception as ex1:
        print('Problem at selection method',ex1)
        embed()
        exit()

    try:
        posM,valM=detect_peaks(FK,size=5,num_peaks=num_sources)
        streamF_T=np.fft.rfft(stream,axis=1)
        streamF=streamF_T*0
        streamF[:,fN.astype('int')] = streamF_T[:,fN.astype('int')]
        stream_filter=np.fft.irfft(streamF)
        resN=[]
        ss_ch,ss_len=stream.shape
        win_lf=float((ss_len/sps))
        for ns_i in range(num_sources):
            if np.isnan(posM[ns_i][0]) or np.isnan(posM[ns_i][1]):
                #print 'no min found'
                resN.append([0,0,-1,0,0,0,0,0,0])
            else:
                # work on errorbars
                #FK_resp=_array_response_short(posM[ns_i][0],posM[ns_i][1],mult_vectors,win_lf,sps,slow,x,y,ss_ch,fmin=freqN[0],fmax=freqN[-1])
                #contours = measure.find_contours(FK_resp/FK_resp.max(), 0.5)
                RES=procPEAKS(posM[ns_i],slow)
                beamB=tdelay(stream_filter.copy(),sps, RES[0], RES[1],x,y)
                x_beam = np.sum(beamB, axis=0)/len(beamB)
                x_beam_rms = np.sqrt(np.sum(x_beam**2)/len(x_beam))
                fstat = bfstat2(beamB)
                corr= corrp(beamB)
                sx=RES[2]
                sy=RES[3]
                resN.append([RES[0],RES[1],fstat,corr[0],corr[1],corr[2],x_beam_rms,sx,sy])
    except Exception as ex1:
        print("there is a problem with the processing at estimation level (inside fkPROC):",ex1)
        embed()
        exit()
    #pdb.set_trace()
    try:
        if func is not None:
            sys.path.append('./')
            from func_fk import func_fk
            func_fk(stream,sps,slow,mult_vectors,fN,x,y,timeSTAMP,func,freqN,aux=aux)
        return resN
    except Exception as ex1:
        print("there is a problem with the processing at outside function (inside fkPROC):",ex1)
        embed()
        exit()


def getXY_array(stream):
    '''
    Returns the site coordinates for a specific array in the format required
    by fk
    '''
    lo=[]
    la=[]
    for ii in range(len(stream)):
        try:
            stream[ii].stats.coordinates
            coord=0
        except Exception as ex1:
            #print('No coordinates')
            coord=1
        if coord==0:
            lo.append(stream[ii].stats.coordinates.longitude)
            la.append(stream[ii].stats.coordinates.latitude)
        else:
            lo.append(stream[ii].stats.sac['stlo'])
            la.append(stream[ii].stats.sac['stla'])
    lo=np.asarray(lo)
    la=np.asarray(la)
    #embed()
    X = []; Y = []
    loc = np.array([la,lo])
    loc = np.mean(loc,axis=1)
    #print loc
    for i in range(0,len(la)):
        try:
            tempDX,tempDY=utlGeoKm(np.mean(lo), np.mean(la), lo[i], la[i])
        except :
            from obspy.signal.util import util_geo_km
            tempDX,tempDY=util_geo_km(np.mean(lo), np.mean(la), lo[i], la[i])
        X.append(tempDX)
        Y.append(tempDY)
    return np.asarray(X), np.asarray(Y)

def dist_az(X1,X2):
    '''
    Returns the angular distance and direction between two locations
    '''

    lat1 = float(X1[0,]) * math.pi/180.
    lat2 = float(X2[0,]) * math.pi/180.
    lon_dist = (float(X1[1,]) - float(X2[1,])) * math.pi/180.
    lat1 = math.atan2 ( math.cos(lat1), 0.99330552180201 * math.sin(lat1) )
    lat2 = math.atan2 ( math.cos(lat2), 0.99330552180201 * math.sin(lat2) )

    s1 = math.sin(lat1)
    s2 = math.sin(lat2)
    c1 = math.cos(lat1)
    c2 = math.cos(lat2)
    sd = math.sin(lon_dist)
    cd = math.cos(lon_dist)
    D = math.acos(c1*c2 + s1*s2*cd) * 180./math.pi
    A = math.atan2(-s2*sd, s1*c2 - c1*s2*cd) * 180./math.pi
    A = A + 360.*(A<0)
    return D, A

def tdelay(data,sps, bazimuth, slowness,x,y):
    '''
    Returns the time-delay-and-sum beam of data for a particular window (FK sliding window)

    Inputs:
    - data is the non-time-aligned waveform in a particular window (FK sliding window), this may be already in real values
    - samprate is the sampling rate
    - x is the array of x coordinates
    - y is the array of y coordinates
    - azimuth is the azimuth for calculating the beam
    - slowness is the slowness for calculating the beam

    Outputs:
    - beam is the time-aligned waveform in a particular window (FK sliding window)
    '''
    '''
    X,Y=getXY_array(data)
    y=np.asarray(Y)
    x=np.asarray(X)
    x = np.array(x); y = np.array(y)
    '''
    [M,N]=np.shape(data)
    dist = np.sqrt(x**2 + y**2)
    bz = 90. - (180./np.pi)*np.arctan2(y,x)
    delay = np.round((dist*slowness) * np.cos((bz-bazimuth)*np.pi/180.)*sps)

    [M,N]=np.shape(data)
    beam=[]

    for m in range(0,M):
        beam.append(np.roll(data[m],int(delay[m])))
    beam=np.asarray(beam)
    return beam

def tdelayT(data,sps, bazimuth, slowness,x,y,slowx,slowy):
    '''
    Returns the time-delay-and-sum beam of data for a particular window (FK sliding window)

    Inputs:
    - data is the non-time-aligned waveform in a particular window (FK sliding window), this may be already in real values
    - samprate is the sampling rate
    - x is the array of x coordinates
    - y is the array of y coordinates
    - azimuth is the azimuth for calculating the beam
    - slowness is the slowness for calculating the beam

    Outputs:
    - beam is the time-aligned waveform in a particular window (FK sliding window)
    '''

    dist = np.sqrt(x**2 + y**2)
    bz = 90. - (180./np.pi)*np.arctan2(y,x)
    delay = np.round((dist*slowness) * np.cos((bz-bazimuth)*np.pi/180.)*sps)

    [M,N]=np.shape(data)
    beam = np.zeros((M,N))
    beam=[]
    for m in range(0,M):
        beam.append(np.roll(data[m],delay[m]))
    beam=np.asarray(beam)
    embed()
    return beam

def bfstat2(beam):
    '''
    Returns the F-statistic using the formalism of Blandford (1974) for a particular window
    (FK sliding window)

    Inputs:
    - beam is the time-aligned waveform in a particular window (FK sliding window)

    Outputs:
    - F is the Blandford F-statistic
    '''
    # normalizing the beam is important for the fstats.

    nchan = float(beam.shape[0])
    inte1=np.sum(beam,axis=0)
    inte2=np.power(inte1,2)
    Num= np.sum(inte2)
    term1=inte1/nchan
    term1=term1
    term1_0 = term1

    for i in range(1,int(nchan)):
        term1 = np.vstack((term1,term1_0))
    try:
        Den = np.sum(np.sum((beam-term1)**2,axis=0))
    except Exception as ex1:
        print('587')
        embed()
        exit()

    F = (nchan-1)*Num/(nchan*Den)
    #embed()
    return F

def bfstatT(beam):
    '''
    Returns the F-statistic using the formalism of Blandford (1974) for a particular window
    (FK sliding window)

    Inputs:
    - beam is the time-aligned waveform in a particular window (FK sliding window)

    Outputs:
    - F is the Blandford F-statistic
    '''
    # normalizing the beam is important for the fstats.

    nchan = float(beam.shape[0])
    for i in range(nchan):
         beam[i,:]=beam[i,:]/np.nanmax(beam[i,:])
    inte1=np.sum(beam,axis=0)
    inte2=np.power(inte1,2)
    Num= np.sum(inte2)
    term1=inte1/nchan
    term1=term1/np.max(term1)
    term1_0 = term1

    for i in range(1,nchan):
        term1 = np.vstack((term1,term1_0))
    try:
        Den = np.sum(np.sum((beam-term1)**2,axis=0))
    except Exception as ex1:
        print('587')
        embed()
        exit()

    F = (nchan-1)*Num/(nchan*Den)
    embed()
    return F

def bfstat(beam):
    '''
    Returns the F-statistic using the formalism of Blandford (1974) for a particular window
    (FK sliding window)

    Inputs:
    - beam is the time-aligned waveform in a particular window (FK sliding window)

    Outputs:
    - F is the Blandford F-statistic
    '''
    # normalizing the beam is important for the fstats.

    nchan = float(beam.shape[0])
    for i in range(nchan):
         beam[i,:]=beam[i,:]/np.nanmax(beam[i,:])
    inte1=np.sum(beam,axis=0)
    inte2=np.power(inte1,2)
    Num= np.sum(inte2)
    term1=inte1/nchan
    term1_0 = term1

    for i in range(1,nchan):
        term1 = np.vstack((term1,term1_0))
    try:
        Den = np.sum(np.sum((beam-term1)**2,axis=0))
    except Exception as ex1:
        print('587')
        embed()
        exit()

    F = (nchan-1)*Num/(nchan*Den)
    #print F,Num,Den
    #embed()
    return F

def corrp(beam):
    '''
    Computes the maximum average cross correlation from all triplets of elements in an
    array for a particular window (FK sliding window)

    Inputs:
    - beam is the time-aligned waveform in a particular window (FK sliding-window)

    Outputs:
    - C is the minimum, maximum, mean average cross-correlation
    '''
    num_sensors = beam.shape[0]
    sensors = list(range(0,num_sensors))
    combs = list(itertools.combinations(sensors,3))
    num_combs = len(combs)

    C_t = np.zeros(num_combs)
    for i in range(0,num_combs):
        C12 = np.corrcoef(beam[combs[i][0],:],beam[combs[i][1],:])
        C13 = np.corrcoef(beam[combs[i][0],:],beam[combs[i][2],:])
        C23 = np.corrcoef(beam[combs[i][1],:],beam[combs[i][2],:])
        C_t[i] = (C12[0,1] + C13[0,1] + C23[0,1])/3.0

    if len(C_t[C_t>0])==0:
        #pdb.set_trace()
        return [0,0,0]
    else:
        return [np.min(C_t[C_t>0]),np.max(C_t[C_t>0]),np.mean(C_t[C_t>0])]

def corrpT(beam):
    '''
    Computes the maximum average cross correlation from all triplets of elements in an
    array for a particular window (FK sliding window)

    Inputs:
    - beam is the time-aligned waveform in a particular window (FK sliding-window)

    Outputs:
    - C is the minimum, maximum, mean average cross-correlation
    '''
    num_sensors = beam.shape[0]
    sensors = list(range(0,num_sensors))
    combs = list(itertools.combinations(sensors,3))
    num_combs = len(combs)

    C_t = np.zeros(num_combs)
    for i in range(0,num_combs):
        C12 = np.corrcoef(beam[combs[i][0],:],beam[combs[i][1],:])
        C13 = np.corrcoef(beam[combs[i][0],:],beam[combs[i][2],:])
        C23 = np.corrcoef(beam[combs[i][1],:],beam[combs[i][2],:])
        C_t[i] = (C12[0,1] + C13[0,1] + C23[0,1])/3.0
        if C_t[i]<0.1:
            embed()
    embed()
    if len(C_t[C_t>0])==0:
        return [0,0]
    else:
        return [np.min(C_t[C_t>0]),np.max(C_t[C_t>0]),np.mean(C_t[C_t>0])]
