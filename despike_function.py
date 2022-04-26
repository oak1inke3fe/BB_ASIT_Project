# -*- coding: utf-8 -*-
"""
Created on Wed Mar 23 12:52:14 2022

@author: oaklin keefe
"""

# purpose: to remove outlier points (aka: bad points)

# inputs
# Y     NxM matrix (N = # of rows, M = number of columns)
# P     float; arbitrary value of how many points makes a bad point?

# vars created inside function
# N     float, number of rows in matrix Y
# col   float, number of columns in matrix Y
# t     row vector, 1xN-1
# bps   float, number of bad points in select column of input matrix Y
# j     column vector, each row is number of good points per column of input matrix Y

# outputs
# BPS   column vector; each row is number of bad points per column of input matrix Y

import numpy as np


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

#detrending
# from scipy import signal
# x_detrended = signal.detrend(x)
# urot = signal.detrend(u)