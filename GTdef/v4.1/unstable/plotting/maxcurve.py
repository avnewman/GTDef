#!/usr/bin/env python3
# v1.5.1 redone to normalize by own parameter, not oposite.
# AVN 22 Jan 2021

import pandas as pd
from numpy import max,min,abs,gradient
import numpy as np
import sys

fin=sys.argv[1]

#invout=pd.read_table('100kmchkbd_40kmstat_inv_inv.out',sep=r"\s+", header=1, names=['beta', 'kappa', 'data_num', 'slip_num', 'ndf', 'rss', 'rms', 'wrrs', 'wrms', 'chi2', 'rchi2', 'r1d', 'r2d', 'strain'])
invout=pd.read_table(fin,sep=r"\s+", header=1, comment='#', names=['beta', 'kappa', 'data_num', 'slip_num', 'ndf', 'rss', 'rms', 'wrrs', 'wrms', 'chi2', 'rchi2', 'r1d', 'r2d', 'strain'])
#invout=invout[ ~ invout['beta'].str.contains("#")] # remove lines that begin with "#" after reading header in read_table
#get exponent of first rss value 
rss0=int(abs(np.floor(np.log10(invout.loc[0,'rss']))))
dp=rss0+2 # add to more decimal places to get final
#print(dp)
invout['rss']=invout[['rss']].round(decimals=dp) # round rss to desired dp
invout=invout.drop_duplicates(subset=['rss'],keep='last') # removes duplicates in rss
#print(invout['rss'])

x=invout.r2d.values
y=invout.rms.values
#old
#ynorm=y*(max(x)-min(x))
#xnorm=x*(max(y)-min(y))
#new
ynorm=y/(max(y)-min(y))
xnorm=x/(max(x)-min(x))
k=invout.kappa.values


yp=gradient(ynorm,xnorm) #y-prime
ypp=gradient(yp,xnorm) #y-dble prime
curve=np.true_divide(abs(ypp),(1+yp**2)**1.5)  # unsigned max curvature
curve[ ~ np.isfinite(curve)]=0  # set singularities to zero.
invout["curvature"]=curve

imax=curve.argmax()

invout_sm= invout[['kappa','r2d','rms','curvature']]
print("# ",imax,k[imax],x[imax],y[imax], max(curve))
print("# Curvature from rms/range(roughness) and roughness/range(rms) ")
print("# ",invout_sm.columns.values)
print(invout_sm.to_string(index=False,header=False))

## interpolated solution
#from scipy.interpolate import interp1d
#from numpy import linspace
# nx=100
# f=interp1d(x,y,kind='cubic')
# xi=linspace(min(x),max(x),num=nx, endpoint=True)
# dx=(max(x)-min(x)/nx)
# yi=f(xi)
# yip=gradient(yi,dx)
# yipp=gradient(yip,dx)

