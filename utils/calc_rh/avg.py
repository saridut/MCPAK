#!/usr/bin/env python

'''
Block averages time-series data.

'''

import sys
import os
import math
import numpy as np
import pandas as pd

fn_ts = 'tmp.txt'
tmin = 0
tmax = None
nblks = 12 #Number of blocks

df = pd.read_csv(fn_ts, skipinitialspace=True)
df_ss = df[df.index > tmin]
if tmax is not None:
    df_ss = df_ss[df_ss.index < tmax]
mean = df_ss.mean(axis=0) 
#std = df_ss.std(axis=0)

bs = df_ss.shape[0]//nblks #Block size
ncols = df_ss.shape[1] #Number of columns (except the index col 'time')
blkdata = np.zeros((nblks,ncols))

#Block data
n = df_ss.shape[0]
for i in range(nblks):
    jbeg = i*bs; jend = jbeg + bs
    blkdata[i,:] = df_ss.iloc[jbeg:jend,:].mean(axis=0)

std = np.std(blkdata, axis=0)

for col in range(ncols):
    print('{0:8s}  {1:14.6g},  {2:14.6g}'.format(df.columns[col], mean[col], std[col]))
    #print('{0:.6g},{1:.6g},'.format(mean[col], std[col]), end='')
print()
