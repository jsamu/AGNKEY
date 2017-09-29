#!/dark/usr/anaconda/bin/python

from astropy.io import fits as pyfits
import numpy as np
import pylab as pl
import glob
import os

pl.ion()
os.chdir('./proc/')
color = 'cmkgbyr'*20
imglist = glob.glob('*1438*.fits')
imglist.sort()
for jj,img in enumerate(imglist):
#    hdulist = pyfits.open(img)
#    data = hdulist[0].data
#    header = hdulist[0].header
#    xx = np.arange(len(data))

    ar, hd = pyfits.getdata(img, header=True)
    xx = np.arange(len(ar[0][0]))
    ll = hd['CRVAL1']+ hd['CD1_1'] * xx
    fl1 = ar[0][0]
    pl.plot(ll,ar[0][0]+jj*2E-16,'-',color=color[jj])
    print img
    raw_input('go on')

