# -*- coding: utf-8 -*-
"""
Created on Tue Mar 29 15:39:09 2022

@author: Oaklin Keefe
"""
# imports
import os
import numpy as np
import pandas as pd
import math
import natsort
from scipy import interpolate
from scipy import signal
import csv
import time

# FUNCTIONS

def despikeSimplest(df,P):
    Y = df.to_numpy()
    N = len(Y) #N number of rows in matrix Y
    col = len(Y[0]) #M number of columns in matrix Y
    t = range(0,N-1)
    BPS = np.zeros((N,1)) # column vector of N rows x 1 column
    for i in range(0, col):
        for iter in range(0,col): #iterate through 3 times
            M = np.nanmedian(Y[i,:])
            S = np.nanstd(Y[i,:])
            j = np.argwhere(Y[i,:]<M+P*S & Y[i,:]>M-P*S & ~np.isnan(Y[i,:]))
            if iter ==1:
                bps = N-len(j) #bad points
            else:
                continue
            if len(j) != 0:
                Y[i,:] = np.interp(t[j],Y[i,j],t,'nearest', 'extrap')
            else:
                continue
        BPS[i]=bps
    return BPS, Y

# Function for aligning the U,V,W coordinaes to the mean wind direction
### function start
#######################################################################################
def alignwind(wind_df):
    try: 
        wind_df = wind_df.replace('NAN', np.nan)
        wind_df['u'] = wind_df['u'].astype(float)
        wind_df['v'] = wind_df['v'].astype(float)
        wind_df['w'] = wind_df['w'].astype(float)
        wind_df = wind_df[(wind_df['u'] > -10)&(wind_df['u'] < 10)]
        Ub = wind_df['u'].mean()
        Vb = wind_df['v'].mean()
        Wb = wind_df['w'].mean()
        Sb = math.sqrt((Ub**2)+(Vb**2))
        beta = math.atan2(Wb,Sb)
        alpha = math.atan2(Vb,Ub)
        x1 = wind_df.index
        x = np.array(x1)
        Ur = wind_df['u']*math.cos(alpha)*math.cos(beta)+wind_df['v']*math.sin(alpha)*math.cos(beta)+wind_df['w']*math.sin(beta)
        Ur_arr = np.array(Ur)
        Vr = wind_df['u']*(-1)*math.sin(alpha)+wind_df['v']*math.cos(alpha)
        Vr_arr = np.array(Vr)
        Wr = wind_df['u']*(-1)*math.cos(alpha)*math.sin(beta)+wind_df['v']*(-1)*math.sin(alpha)*math.sin(beta)+wind_df['w']*math.cos(beta)     
        Wr_arr = np.array(Wr)
        T_arr = np.array(wind_df['T'])
        u_arr = np.array(wind_df['u'])
        v_arr = np.array(wind_df['v'])
        w_arr = np.array(wind_df['w'])

        fu_p = interpolate.interp1d(x, Ur_arr,fill_value="extrapolate") #interpolate u_p
        fv_p = interpolate.interp1d(x, Vr_arr,fill_value="extrapolate") #interpolate v_p
        fw_p = interpolate.interp1d(x, Wr_arr,fill_value="extrapolate") #interpolate w_p
        ft_p = interpolate.interp1d(x, T_arr,fill_value="extrapolate") #interpolate T (aka T_p)
        fu = interpolate.interp1d(x, u_arr,fill_value="extrapolate") #interpolate u
        fv = interpolate.interp1d(x, v_arr,fill_value="extrapolate") #interpolate v
        fw = interpolate.interp1d(x, w_arr,fill_value="extrapolate") #interpolate w

        sonic_xnew = np.arange(0, 24000)
        sonic_ynew_up = fu_p(sonic_xnew)   # use interpolation function returned by `interp1d`
        sonic_ynew_vp = fv_p(sonic_xnew)
        sonic_ynew_wp = fw_p(sonic_xnew)
        sonic_ynew_tp = ft_p(sonic_xnew)
        sonic_ynew_u = fu(sonic_xnew)
        sonic_ynew_v = fv(sonic_xnew)
        sonic_ynew_w = fw(sonic_xnew)

        S = Sb
        b = beta*180/math.pi
        a = alpha*180/math.pi
        df_aligned = pd.DataFrame({'base_index':sonic_xnew,'u':sonic_ynew_u, 'v':sonic_ynew_v, 'w':sonic_ynew_u,'u_p':sonic_ynew_up, 'v_p':sonic_ynew_vp, 'w_p':sonic_ynew_wp, 'T_p':sonic_ynew_tp})

    except:
        pass
    return S,b,a,df_aligned
#######################################################################################
### function end


N = 6        # N = number of sonics
M = 10       # M = number of 20-minute files
Oo = 1000    # O = number of lines per 20-min file

df_UW = pd.DataFrame()
df_WT = pd.DataFrame()
df_WE = pd.DataFrame()
df_U = pd.DataFrame()
df_wDir = pd.DataFrame()
df_U_p_sqr = pd.DataFrame()
df_V_p_sqr = pd.DataFrame()
df_W_p_sqr = pd.DataFrame()
iterator = 1
for n in range(1,N):
    U = []
    Ustream = []
    ALPHA = []
    uw = []
    UW =[]
    wT = []
    WT = []
    tke = []
    we = []
    WE = []
    u_p_sqr = []
    U_p_sqr = []
    v_p_sqr = []
    V_p_sqr = []
    w_p_sqr = []
    W_p_sqr = []
    for m in range(1,M):
        P=6
        BPS, UVW = despikeSimplest(UVW,P)
        BPS, Tv = despikeSimplest(Tv,P)
        df_UVWT = pd.Dataframe({'u':UVW[:,0], 'v':UVW[:,1], 'w':UVW[:,2], 'T':Tv[:,0]})
        S, b, a, df_aligned = alignwind(df_UVWT)
        wind_sp = np.nanmean(df_aligned['u_p'])
        urot = signal.detrend(df_aligned)
        
        u_prime = df_aligned['u_p']
        v_prime = df_aligned['v_p']
        w_prime = df_aligned['w_p']
        T_prime = df_aligned['T_p']
        u = np.array(df_aligned['u'])
        v = np.array(df_aligned['v'])
        w = np.array(df_aligned['w'])
        for i in range(1,Oo):
            uw_i = u_prime[i]*w_prime[i]
            uw.append(uw_i)
            wT_i = w_prime[i]*T_prime[i]
            wT.append(wT_i)
            # tke_i = 0.5*((u_prime[i]**2)+(v_prime[i]**2)+(w_prime[i]**2))
            tke_i = 0.5*(((signal.detrend(df_aligned['u_p']))**2)+((signal.detrend(df_aligned['v_p']))**2)+((signal.detrend(df_aligned['w_p']))**2))
            tke.append(tke_i)
            we_i = tke_i*w_prime[i]
            we.append(we_i)
            u_p_sqr_i = u_prime[i]**2
            u_p_sqr.append(u_p_sqr_i)
            v_p_sqr_i = v_prime[i]**2
            v_p_sqr.append(v_p_sqr_i)
            w_p_sqr_i = w_prime[i]**2
            w_p_sqr.append(w_p_sqr_i)
            streamwise_i = math.sqrt((u**2)+(v**2)+(w**2))
        ALPHA_i = np.array(a)
        ALPHA.append(ALPHA_i)
        U_i = u.mean()
        U.append(U_i)
        Ustream_i = streamwise_i.mean()
        Ustream.append(Ustream_i)
        uw_arr = np.array(uw)
        UW_i = uw_arr.mean()
        UW.append(UW_i)
        uw.clear()
        wT_arr = np.array(wT)
        WT_i = wT_arr.mean()
        WT.append(WT_i)
        wT.clear()
        we_arr = np.array(we)
        WE_i = we_arr.mean()
        WE.append(WE_i)
        we.clear()
        u_p_sqr_arr = np.array(u_p_sqr)
        U_p_sqr_i = u_p_sqr_arr.mean()
        U_p_sqr.append(U_p_sqr_i)
        u_p_sqr.clear()
        v_p_sqr_arr = np.array(v_p_sqr)
        V_p_sqr_i = v_p_sqr_arr.mean()
        V_p_sqr.append(V_p_sqr_i)
        v_p_sqr.clear()
        w_p_sqr_arr = np.array(w_p_sqr)
        W_p_sqr_i = w_p_sqr_arr.mean()
        W_p_sqr.append(W_p_sqr_i)
        w_p_sqr.clear()
        

    UW_col_name = 'UW_sonic' + str(iterator)
    WT_col_name = 'WT_sonic' + str(iterator)
    WE_col_name = 'WE_sonic' + str(iterator)
    U_col_name = 'U_sonic' + str(iterator)
    wDir_col_name = 'wDir_sonic' + str(iterator)
    U_p_sqr_col_name = 'UpSqr_sonic' + str(iterator)
    V_p_sqr_col_name = 'VpSqr_sonic' + str(iterator)
    W_p_sqr_col_name = 'WpSqr_sonic' + str(iterator)
    
    df_UW.loc[:,UW_col_name] = UW
    df_WT.loc[:,WT_col_name] = WT
    df_WE.loc[:,WE_col_name] = WE
    df_U.loc[:,U_col_name] = U
    df_wDir.loc[:,wDir_col_name] = ALPHA
    df_U_p_sqr.loc[:,U_p_sqr_col_name] = U_p_sqr
    df_V_p_sqr.loc[:,V_p_sqr_col_name] = V_p_sqr
    df_W_p_sqr.loc[:,W_p_sqr_col_name] = W_p_sqr  
    iterator += 1 
