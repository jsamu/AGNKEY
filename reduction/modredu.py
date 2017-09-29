#!/dark/usr/anaconda/bin/python
import numpy as np
import glob
import shutil
import pylab as pl
from pyraf import iraf
from astropy.io import fits as pyfits
import os
from iraf import noao
from iraf import imred
from iraf import ccdred
from iraf import specred


#Function to overscan, trim, bias, and normalize STD, SCI, and LAMP

def zcomb(imglist,name,_gain,_ron):
    iraf.zerocombine(imglist,output=name,combine='average',reject='minmax',ccdtype='',process='no',gain=_gain,rdnoise=_ron,Stdout=1)

def fcomb(imglist,name,_gain,_ron):
    iraf.flatcombine(imglist,output=name,combine='average',reject='avsigclip',ccdtype='',process='no',subsets='yes',delete='no',clobber='no',gain=_gain,rdnoise=_ron,Stdout=1)

def proc(fits,name,ovrsc,trim,zcor,fcor,_zero,_flat):
    iraf.ccdproc(fits,output=name,ccdtype='',noproc='no',fixpix='no',overscan='yes',trim='yes',zerocor=zcor,darkcor='no',flatcor=fcor,illumcor='no',fringecor='no',readcor='no',scancor='no',biassec=ovrsc,trimsec=trim,zero=_zero,flat=_flat,Stdout=1)

def imrep(img,_area,_value):
    iraf.imreplace(images=img+_area,value=_value)

def response(cal,resp,_order):
    iraf.response(calibration=cal,normalization=cal,response=resp,interactive='no',function='legendre',order=_order,graphics='stdgraph')

def exsci(fits,name,_ron,_gain):
    iraf.specred.apall(fits,output=name,interactive='no',find='yes',recenter='yes',resize='yes',edit='no',trace='yes',fittrace='no',extract='yes',extras='no',review='no',line=600,nsum=20,lower=-5.,upper=5.,b_function='chebyshev',b_order=1,b_sample='-20:-10,10:20',b_naverage=3.,b_niterate=0,b_low_reject=3.,b_high_reject=3.,b_grow=0,width=10.,radius=10.,nfind=1,background='fit',readnoise=_ron,gain=_gain)

def extract(fits,name,ref,dsum,_ron,_gain):
    iraf.apsum(fits,output=name,references=ref,nsum=dsum,interactive='no',find='no',edit='no',trace='no',fittrace='no',extract='yes',review='no',readnoise=_ron,gain=_gain)

def identify(fits,linelist,width,radius):
    iraf.identify(fits,coordlist=linelist,fwidth=width,cradius=radius)

def refspec(speclist,ref):
    iraf.refspectra(speclist,references=ref,sort='',group='')

def dispcor(spec,name):
    iraf.dispcor(spec,output=name)

def standard(speclist,_std,extinc,_caldir,observ,starname,_airmass):
    iraf.standard(speclist,_std,extinction=extinc,caldir=_caldir,observatory=observ,star_name=starname,airmass=_airmass)

def sens(extinct,observ):
    iraf.sens(standards='std',sensitivity='sens',extinction=extinct,observatory=observ)

def calib(spec,redspec,extinct,observ,_airmass):
    iraf.calib(spec,redspec,extinction=extinct,observatory=observ,airmass=_airmass)


"""

#Sort images by chip number and by image type

if not os.path.isdir('CHIP2'):
    os.makedirs('CHIP2')
else:
    print 'CHIP2 directory already present'

if not os.path.isdir('CHIP1'):
    os.makedirs('CHIP1')
else:
    print 'CHIP1 directory already present'


imglist = glob.glob('FORS2*.fits')
for img in imglist:
    hdu = pyfits.getheader(img)
    if 'EXTNAME' in hdu:
        if hdu['EXTNAME'] == 'CHIP2':
            os.system('mv '+img+' CHIP2/')
        elif hdu['EXTNAME'] == 'CHIP1':
            os.system('mv '+img+' CHIP1/')
        else:
            print 'EXTNAME not available'
"""

os.chdir('CHIP1')



imglist = glob.glob('FORS2*.fits')
dictionary = {}


for img in imglist:
    hdu = pyfits.getheader(img)
    if 'OBJECT' in hdu:
        obj = hdu['OBJECT']
        if obj not in dictionary:
            dictionary[obj] = [img]
        else:
            dictionary[obj].append(img)



#Combine bias images into master bias image

if os.path.isfile('masterbias.fits'):
    os.remove('masterbias.fits')

f = open('biaslist','w')
for img in dictionary['BIAS']:
    f.write(img+'\n')
f.close()

imglist = '@biaslist'
name = 'masterbias.fits'
hdb = pyfits.getheader(dictionary['BIAS'][0])
_gain = hdb['HIERARCH ESO DET OUT1 GAIN']
_ron = hdb['HIERARCH ESO DET OUT1 RON']

zcomb(imglist,name,_gain,_ron)



#Combine flat images into master flat image

if os.path.isfile('masterflat.fits'):
    os.remove('masterflat.fits')

f = open('flatlist','w')
for img in dictionary['FLAT,LAMP']:
    f.write(img+'\n')
f.close()

imglist = '@flatlist'
name = 'masterflat.fits'
hdf = pyfits.getheader(dictionary['FLAT,LAMP'][0])
_gain = hdf['HIERARCH ESO DET OUT1 GAIN']
_ron = hdf['HIERARCH ESO DET OUT1 RON']

fcomb(imglist,name,_gain,_ron)



#Selection of overscan and trim regions

ovrsc = '[1:4,*]'
trim = '[200:2046,50:200]'

if os.path.isfile('tmasterbias.fits'):
    os.remove('tmasterbias.fits')

fits = 'masterbias.fits'
name = 'tmasterbias.fits'
zcor = 'no'
fcor = 'no'
_zero = ''
_flat = ''
proc(fits,name,ovrsc,trim,zcor,fcor,_zero,_flat)




if os.path.isfile('tmasterflat.fits'):
    os.remove('tmasterflat.fits')

fits = 'masterflat.fits'
name = 'tmasterflat.fits'
zcor = 'yes'
fcor = 'no'
_zero = 'tmasterbias.fits'
_flat = ''
proc(fits,name,ovrsc,trim,zcor,fcor,_zero,_flat)




#Normalization of flat field master image

if os.path.isfile('ntmasterflat.fits'):
    os.remove('ntmasterflat.fits')

cal = 'tmasterflat.fits'
resp = 'ntmasterflat.fits'
order = 100
response(cal,resp,order)



#Replacement of pixels in ntmasterflat

img = 'ntmasterflat'
_area = '[1:100,*]'
_value = 1.0
imrep(img,_area,_value)


#Applying overscan, trim, bias, and normalization to science, standard, and lamp frames

alist = glob.glob('AGN*.fits')
for f in alist:
    os.remove(f)
slist = glob.glob('STD*.fits')
for f in slist:
    os.remove(f)
llist = glob.glob('LAMP*.fits')
for f in llist:
    os.remove(f)


zcor = 'yes'
fcor = 'yes'
_zero = 'tmasterbias.fits'
_flat = 'ntmasterflat.fits'

for key in dictionary.keys():
    if 'AGN' in key:
        for img in dictionary[key]:
            hds = pyfits.getheader(img)
            if 'EXPTIME' in hds:
                if hds['EXPTIME'] >= 50:
                    num = img.rsplit('.',2)[1]
                    name = key+'_'+num
                    proc(img,name,ovrsc,trim,zcor,fcor,_zero,_flat)
            else:
                print 'No exposure time info'

    elif 'STD' in key:
        for img in dictionary[key]:
            num = img.rsplit('.',2)[1]
            name = key+'_'+num
            proc(img,name,ovrsc,trim,zcor,fcor,_zero,_flat)

    elif 'WAVE' in key:
        for img in dictionary[key]:
            num = img.rsplit('.',2)[1]
            name = 'LAMP_'+num
            proc(img,name,ovrsc,trim,zcor,fcor,_zero,_flat)



#Extract spectrum of science, standard, and lamp frames

if os.path.isdir('database'):
    shutil.rmtree('database')

scilist = glob.glob('AGN*.fits')
for img in scilist:
    hds = pyfits.getheader(img)
    _gain = hds['HIERARCH ESO DET OUT1 GAIN']
    _ron = hds['HIERARCH ESO DET OUT1 RON']
#    dline = 600
#    dsum = 20
    name = img.rsplit('.fits',1)[0]+'_EX'
    exsci(img,name,_ron,_gain)

"""
stdlist = glob.glob('STD*.fits')
lamplist = glob.glob('LAMP*.fits')
sllist = stdlist + lamplist
sci = scilist[0]
ref = sci.rsplit('.fits',1)[0]
for img in sllist:
    hdl = pyfits.getheader(img)
    _gain = hdl['HIERARCH ESO DET OUT1 GAIN']
    _ron = hdl['HIERARCH ESO DET OUT1 RON']
    dsum = 10
    name = img.rsplit('.fits',1)[0]+'_EX'
    extract(img,name,ref,dsum,_ron,_gain)



#Identify the spectral lines in LAMP_EX

lexlist = glob.glob('LAMP*EX.fits')
for img in lexlist:
    linelist = 'Lines_HgCdHeNeAr600.dat'
    width = 10.0
    radius = 15.0
    identify(img,linelist,width,radius)



#Correlate spectral lines with dispersion axes and apply sensitivity function to science frame

scispec = glob.glob('AGN*EX.fits')
stdspec = glob.glob('STD*EX.fits')
speclist = scispec + stdspec
lexlist = glob.glob('LAMP*EX.fits')
ref = lexlist[0]
ref = ref.rsplit('.fits',1)[0]

print ref

for spec in speclist:
    spec = spec.rsplit('.fits',1)[0]
    refspec(spec,ref)
    name = spec.rsplit('_EX',1)[0]+'_EL'
    dispcor(spec,name)


if os.path.isfile('std'):
    os.remove('std')
senslist = glob.glob('sens*.fits')
for f in senslist:
    os.remove(f)


extinct = 'extinction_lasilla.dat'
_caldir = 'onedstds$ctionewcal/'
observ = 'paranal'


stdlist = glob.glob('STD*EL.fits')

for spec in stdlist:
    _std = 'std'
    hds = pyfits.getheader(spec)
    starname = 'eg274'#hds['HIERARCH ESO OBS NAME']
    _airmass = np.average([hds['HIERARCH ESO TEL AIRM START'],hds['HIERARCH ESO TEL AIRM END']])
    standard(spec,_std,extinct,_caldir,observ,starname,_airmass)
    sens(extinct,observ)

scired = glob.glob('AGN*EL.fits')

for spec in scired:
    hdo = pyfits.getheader(spec)
    _airmass = np.average([hdo['HIERARCH ESO TEL AIRM START'],hdo['HIERARCH ESO TEL AIRM END']])
    name = spec.rsplit('_EL',1)[0]+'_LF'
    calib(spec,name,extinct,observ,_airmass)
"""


#    print img
#    iraf.display(img,frame=1,fill='yes')
#    raw_input(' go on ')

