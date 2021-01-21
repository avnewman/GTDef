#!/usr/bin/env python3
# v1.5 redone to normalize by own parameter, not oposite.
# AVN 21 Jan 2021
import pandas as pd
from numpy import max,min,abs,gradient
import sys

fin=sys.argv[1]

#invout=pd.read_table('100kmchkbd_40kmstat_inv_inv.out',sep=r"\s+", header=1, names=['beta', 'kappa', 'data_num', 'slip_num', 'ndf', 'rss', 'rms', 'wrrs', 'wrms', 'chi2', 'rchi2', 'r1d', 'r2d', 'strain'])
invout=pd.read_table(fin,sep=r"\s+", header=1, names=['beta', 'kappa', 'data_num', 'slip_num', 'ndf', 'rss', 'rms', 'wrrs', 'wrms', 'chi2', 'rchi2', 'r1d', 'r2d', 'strain'])

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
curve=abs(ypp)/(1+yp**2)**1.5  # unsigned max curvature
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
