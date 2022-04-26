# -*- coding: utf-8 -*-
"""
Created on Mon Apr  4 09:13:47 2022

@author: oaklin keefe
"""
### Function for Calculating Dissipation: Epsilon
from scipy import signal

nfft = 512
pow1 = 5/3
pow2 = 2/3
alphaU = 0.52
freq_id_1 = max(1,min(1.5*Ustream/zu),4)
freq_id_2 = min(freq_id_1+2,5)
fs = fr1 #this comes from each sonic so we need to figure out how to adjust this
x= np.array(Ustream)

f, Pxx = signal.welch(
                        x, 
                        fs, 
                        window = 'hamming', 
                        nperseg=1024,
                        noverlap = None,
                        nfft= nfft/4, 
                        detrend = 'linear',
                        return_onesided=True,
                        scaling='density',
                        axis=-1,
                        average='mean'
                      ) 
#f = nd array; array of sample frequencies
#Pxx = nd array; Power spectral density (aka, power spectrum of input, x)

imag_component = complex(Pxx) #jim calls this 'i'
if len(imag_component)== 0:
    # clear imag_component
    #below, Jim uses 'i' instead of imag_component
    imag_component = np.argwhere((fr>freq_id_1)&(fr<freq_id_2))
    fSf = Pxx(imag_component)*(fr(imag_component)**pow1)*((2*math.pi/Ustream)**pow2)
    epsilon = (np.nanmean(fSf/alphaU))**1.5
    #clear fSf Pxx fr imag_component
else:
    epsilon = NaN

Cnfft = len(Tvrot)
f, Pxy = signal.csd(
                        x, 
                        y, 
                        fs, 
                        window = 'hamming',
                        nperseg=1024,
                        noverlap=None,
                        nfft = Cnfft,
                        detrend = 'constant',
                        return_onesided=True,
                        scaling='density',
                        axis=-1,
                        average='mean'
                        
                    )
#f = nd array; array of sample frequencies
#Pxx = nd array; Cross spectral density (aka, cross power spectrum of inputs x,y)

Cnfft=length(Tvrot);
[Cspc,fr] = csd2(urotr(:,3),urotr(:,1),Cnfft,fs,hamming(Cnfft),0,'none');
f0=fr(2);
frhigh=max(min(wspd/zu,4.5),1);  %Use fz/U<=1
Cuw = real(Cspc);Cuw(1)=0;clear Cspc
[Cspc,fr] = csd2(urotr(:,3),urotr(:,2),Cnfft,fs,hamming(Cnfft),0,'none');
Cvw = real(Cspc);Cvw(1)=0;clear Cspc
[Cspc,fr] = csd2(urotr(:,3),Tvrot,Cnfft,fs,hamming(Cnfft),0,'none');
CwT = real(Cspc);CwT(1)=0;clear Cspc

j=find(fr>=f0 & fr<=frhigh);       %Total
uwspec = sum(Cuw(j))*f0;
vwspec = sum(Cvw(j))*f0;
wTspec = sum(CwT(j))*f0; 
clear j fr Tvrot urot urotr Cspc Cuw Cvw CwT

