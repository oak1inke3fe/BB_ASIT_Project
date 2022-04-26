# -*- coding: utf-8 -*-
"""
Created on Fri Mar 25 09:09:20 2022

@author: oaklin keefe
"""

# Computing u'w' flux from sonics
## inputs
### Tower Sonic #N (N = number of sonic in project)
# frN       1x1           double
# diagN     1x1728201     double
# diagidN   1x1728201     double
# TvN       1x1728201     double
# UVWN      3x1728201     double
# YdayN     1x1728201     double 

### Mast Sonic #N (N = number of sonic in project)
# frN       1x1           double
# diagN     1x1728201     double 
# diagidN   1x1728201     double
# TvN       1x1728201     double
# UVWN      3x1728201     double
# YdayN     1x1728201     double

## outputs
# uw_NMAX 
# (N= the number of sonics you have; this will be the tallest sonic)
# ...
# uw_N 
# (N = number of sonic above the lowest sonic)
# ...
# uw_1 
# (lowest sonic)

import numpy as np
import scipy.signal as signal
import math
P=6
BPS, UVW = despikeSimplest(UVW,P)
BPS, Tv = despikeSimplest(Tv,P)
U, alpha, beta = alignwind(UVW)
U_h = math.sqrt(U[:,1]*U[:,1]+U[:,2]*U[:,2])
Ustreamwise=math.sqrt(U[:,1]*U[:,1]+U[:,2]*U[:,2]+U[:,3]*U[:,3]);
Ustream_avg=np.nanmean(Ustreamwise)
wspd=np.nanmean(U[1,:])
Urot=signal.detrend(U)
w = Urot[:,3]
u = Urot[:,1]
u_sqr_avg = np.nanmean(Urot[:,1]**2)
v_sqr_avg = np.nanmean(Urot[:,2]**2)
w_sqr_avg = np.nanmean(Urot[:,3]**2)
tke=(Urot[:,1]**2)+(Urot[:,2]**2)+(Urot[:,3]**2)
e=tke/2
uw = np.nanmean(u*w)
we = np.nanmean(w*e)
wTv = np.nanmean(w*Tv)
u_avg = np.nanmean(u)

ux=np.nanmean(UVW[1,:])
uy=np.nanmean(UVW[2,:])
uz=np.nanmean(UVW[3,:])
wind_dir=math.atan3(uy,ux)
wind_dir_avg = np.nanmean(wind_dir)

#epsilon
nfft=512
pow1=5/3
pow2=2/3
alphau=0.52
frid1=max(1,min(1.5*Ustream_avg/zu,4))
frid2=min(frid1+2,5)

f, Pxx_den = signal.welch(x, fs, nperseg=1024)
# [Suu,fr] = psd2(Ustreamwise,nfft,fs,hamming(nfft),nfft/4,'linear')
i=find(imag(Suu))
if isempty(i)
    clear i
    i=find(fr>frid1 & fr<frid2)
    fSf=Suu(i).*(fr(i).^pow1).*((2*pi./Ustream).^pow2)
    epsilon=(nanmean(fSf./alphau)).^1.5
    clear fSf Suu fr i
else
    epsilon=NaN;