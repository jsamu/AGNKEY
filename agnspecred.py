#!/dark/usr/anaconda/bin/python
import numpy as np
import glob
import shutil
import datetime
import os
import sys
#import getopt
#from optparse import OptionParser
import pylab as pl
from pyraf import iraf
from astropy.io import fits as pyfits
from iraf import noao
from iraf import imred
from iraf import ccdred
from iraf import specred
#import cosmics


def imgsort():
    if not os.path.isdir('CHIP1'):
        os.makedirs('CHIP1')
    else:
        pass
        
    if not os.path.isdir('CHIP2'):
        os.makedirs('CHIP2')
    else:
        pass

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





def process():
    ovrsc = '[1:4,*]'
    trim = '[200:2046,50:200]'
    imglist = glob.glob('FORS2*.fits')
#    objdict = {}
#    objdict.clear()
    for img in imglist:
        hdu = pyfits.getheader(img)
        if 'OBJECT' in hdu:
            obj = hdu['OBJECT']
            if obj not in objdict:
                objdict[obj] = [img]
            else:
                objdict[obj].append(img)

    if 'STD' not in objdict:
        os.system('cp /dark/jsamuel/agn/standards/FORS2.2016-05-12T09:15:14.678.fits ./')
        objdict['STD'] = ['FORS2.2016-05-12T09:15:14.678.fits']
        os.system('cp /dark/jsamuel/agn/standards/FORS2.2016-05-12T04:24:38.547.fits ./')
        objdict['STD'].append('FORS2.2016-05-12T04:24:38.547.fits')

    stars = {'specphot-LTT7379':'l7379','specphot-LDS749B':'lds749b','specphot-EG274':'eg274','specphot-G138-31':'g13831','specphot-C-32d9927':'cd32','specphot-LTT9491':'l9491','specphot-LTT7987':'l7987'}

    i = 0
    for key in objdict.keys():
        if 'STD' in key:
            for img in objdict[key]:
                i = i + 1
                numstars = len(objdict['STD'])
                hds = pyfits.getheader(img)
                _starname = stars[hds['HIERARCH ESO OBS NAME']]
                if _starname in ['lds749b','g13831']:
                    print 'Bad standard, copying from 2016-05-12'
                    if not os.path.isdir('badstd'):
                        os.mkdir('badstd')
                    os.system('mv '+img+' ./badstd')
                    if i >= numstars:
                        objdict.pop('STD')
                        os.system('cp /dark/jsamuel/agn/standards/FORS2.2016-05-12T09:15:14.678.fits ./')
                        objdict['STD'] = ['FORS2.2016-05-12T09:15:14.678.fits']
                        os.system('cp /dark/jsamuel/agn/standards/FORS2.2016-05-12T04:24:38.547.fits ./')
                        objdict['STD'].append('FORS2.2016-05-12T04:24:38.547.fits')



    

    if os.path.isfile('biaslist'):
        os.remove('biaslist')

    if os.path.isfile('masterbias.fits'):
        os.remove('masterbias.fits')

    f = open('biaslist','w')
    for img in objdict['BIAS']:
        f.write(img+'\n')
    f.close()

    imglist = '@biaslist'
    name = 'masterbias.fits'
    hdb = pyfits.getheader(objdict['BIAS'][0])
    _gain = hdb['HIERARCH ESO DET OUT1 GAIN']
    _ron = hdb['HIERARCH ESO DET OUT1 RON']
    iraf.zerocombine(imglist,output=name,combine='average',reject='minmax',ccdtype='none',process='no',gain=_gain,rdnoise=_ron,Stdout=1)



    if os.path.isfile('flatlist'):
        os.remove('flatlist')

    if os.path.isfile('sciflatlist'):
        os.remove('sciflatlist')

    if os.path.isfile('stdflatlist'):
        os.remove('stdflatlist')

    if os.path.isfile('masterflat.fits'):
        os.remove('masterflat.fits')

    if os.path.isfile('scimasterflat.fits'):
        os.remove('scimasterflat.fits')

    if os.path.isfile('stdmasterflat.fits'):
        os.remove('stdmasterflat.fits')



        
    f = open('sciflatlist','w')
    for img in objdict['FLAT,LAMP']:
        hdu = pyfits.getheader(img)
        if hdu['HIERARCH ESO DPR TECH'] == 'SPECTRUM':
            f.write(img+'\n')
    f.close()

    j = 0
    f = open('stdflatlist','w')
    for img in objdict['FLAT,LAMP']:
        hdu = pyfits.getheader(img)
        if hdu['HIERARCH ESO DPR TECH'] == 'MOS':
            f.write(img+'\n')
            j = j + 1
    f.close()





    imglist = '@sciflatlist'
    name = 'scimasterflat.fits'
    hdf = pyfits.getheader(objdict['FLAT,LAMP'][0])
    _gain = hdf['HIERARCH ESO DET OUT1 GAIN']
    _ron = hdf['HIERARCH ESO DET OUT1 RON']
    iraf.flatcombine(imglist,output=name,combine='average',reject='avsigclip',ccdtype='none',process='no',subsets='yes',delete='no',clobber='no',gain=_gain,rdnoise=_ron,Stdout=1)


    if j == 0:
        imglist = '@sciflatlist'
    elif j >= 1:
        imglist = '@stdflatlist'
    name = 'stdmasterflat.fits'
    hdf = pyfits.getheader(objdict['FLAT,LAMP'][0])
    _gain = hdf['HIERARCH ESO DET OUT1 GAIN']
    _ron = hdf['HIERARCH ESO DET OUT1 RON']
    iraf.flatcombine(imglist,output=name,combine='average',reject='avsigclip',ccdtype='none',process='no',subsets='yes',delete='no',clobber='no',gain=_gain,rdnoise=_ron,Stdout=1)





    if os.path.isfile('tmasterbias.fits'):
        os.remove('tmasterbias.fits')

    fits = 'masterbias.fits'
    name = 'tmasterbias.fits'
    zcor = 'no'
    fcor = 'no'
    _zero = ''
    _flat = ''
    iraf.ccdproc(fits,output=name,ccdtype='',noproc='no',fixpix='no',overscan='yes',trim='yes',zerocor=zcor,darkcor='no',flatcor=fcor,illumcor='no',fringecor='no',readcor='no',scancor='no',biassec=ovrsc,trimsec=trim,zero=_zero,flat=_flat,Stdout=1)


    if os.path.isfile('tmasterflat.fits'):
        os.remove('tmasterflat.fits')

    if os.path.isfile('tscimasterflat.fits'):
        os.remove('tscimasterflat.fits')

    if os.path.isfile('tstdmasterflat.fits'):
        os.remove('tstdmasterflat.fits')



    fits = 'scimasterflat.fits'
    name = 'tscimasterflat.fits'
    zcor = 'yes'
    fcor = 'no'
    _zero = 'tmasterbias.fits'
    _flat = ''
    iraf.ccdproc(fits,output=name,ccdtype='',noproc='no',fixpix='no',overscan='yes',trim='yes',zerocor=zcor,darkcor='no',flatcor=fcor,illumcor='no',fringecor='no',readcor='no',scancor='no',biassec=ovrsc,trimsec=trim,zero=_zero,flat=_flat,Stdout=1)

    fits = 'stdmasterflat.fits'
    name = 'tstdmasterflat.fits'
    zcor = 'yes'
    fcor = 'no'
    _zero = 'tmasterbias.fits'
    _flat = ''
    iraf.ccdproc(fits,output=name,ccdtype='',noproc='no',fixpix='no',overscan='yes',trim='yes',zerocor=zcor,darkcor='no',flatcor=fcor,illumcor='no',fringecor='no',readcor='no',scancor='no',biassec=ovrsc,trimsec=trim,zero=_zero,flat=_flat,Stdout=1)
                




    if os.path.isfile('ntmasterflat.fits'):
        os.remove('ntmasterflat.fits')

    if os.path.isfile('ntscimasterflat.fits'):
        os.remove('ntscimasterflat.fits')

    if os.path.isfile('ntstdmasterflat.fits'):
        os.remove('ntstdmasterflat.fits')

    cal = 'tscimasterflat.fits'
    resp = 'ntscimasterflat.fits'
    _order = 100
    iraf.response(calibration=cal,normalization=cal,response=resp,interactive='no',function='legendre',order=_order,graphics='stdgraph')

    cal = 'tstdmasterflat.fits'
    resp = 'ntstdmasterflat.fits'
    _order = 100
    iraf.response(calibration=cal,normalization=cal,response=resp,interactive='no',function='legendre',order=_order,graphics='stdgraph')




##########

    ovrsc = '[1:4,*]'
    trim = '[200:2046,50:200]'
    
    fimg = 'ntscimasterflat'
    _area = '[1:150,*]'
    _value = 1.0
    iraf.imreplace(images=fimg+_area,value=_value)

    fimg = 'ntstdmasterflat'
    _area = '[1:150,*]'
    _value = 1.0
    iraf.imreplace(images=fimg+_area,value=_value)
    

    alist = glob.glob('*AGN*.fits')
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
    sciflat = 'ntscimasterflat.fits'
    stdflat = 'ntstdmasterflat.fits'

    for key in objdict.keys():
        if 'AGN' in key:
            for img in objdict[key]:
                hds = pyfits.getheader(img)
                if 'EXPTIME' in hds:
                    if hds['EXPTIME'] >= 50.0:
                        num = img.rsplit('.',2)[1]
                        name = key+'_'+num
                        '''
                        iraf.ccdproc(img,output='trim_'+name,ccdtype='',noproc='no',fixpix='no',overscan='yes',trim='yes',zerocor='no',darkcor='no',flatcor='no',illumcor='no',fringecor='no',readcor='no',scancor='no',biassec=ovrsc,trimsec=trim,zero='',flat='',Stdout=1)

                        _gain = hds['HIERARCH ESO DET OUT1 GAIN']
                        _ron = hds['HIERARCH ESO DET OUT1 RON']
                        cosmics.lacos('trim_'+name+'.fits', output='c_'+name+'.fits', gain=_gain, readn=_ron, xorder=9, yorder=9, sigclip=4.5, sigfrac=0.5, objlim=1, verbose=True, interactive=False)

                        iraf.ccdproc('c_'+name,output=name,ccdtype='',noproc='no',fixpix='no',overscan='no',trim='no',zerocor=zcor,darkcor='no',flatcor=fcor,illumcor='no',fringecor='no',readcor='no',scancor='no',biassec='',trimsec='',zero=_zero,flat=sciflat,Stdout=1)
                        '''

                        iraf.ccdproc(img,output=name,ccdtype='',noproc='no',fixpix='no',overscan='yes',trim='yes',zerocor=zcor,darkcor='no',flatcor=fcor,illumcor='no',fringecor='no',readcor='no',scancor='no',biassec=ovrsc,trimsec=trim,zero=_zero,flat=sciflat,Stdout=1)

                    else:
                        pass


        elif 'STD' in key:
            for img in objdict[key]:
                num = img.rsplit('.',2)[1]
                name = key+'_'+num
                iraf.ccdproc(img,output=name,ccdtype='',noproc='no',fixpix='no',overscan='yes',trim='yes',zerocor=zcor,darkcor='no',flatcor=fcor,illumcor='no',fringecor='no',readcor='no',scancor='no',biassec=ovrsc,trimsec=trim,zero=_zero,flat=stdflat,Stdout=1)

        elif 'WAVE' in key:
            for img in objdict[key]:
                num = img.rsplit('.',2)[1]
                name = 'LAMP_'+num
                iraf.ccdproc(img,output=name,ccdtype='',noproc='no',fixpix='no',overscan='yes',trim='yes',zerocor=zcor,darkcor='no',flatcor=fcor,illumcor='no',fringecor='no',readcor='no',scancor='no',biassec=ovrsc,trimsec=trim,zero=_zero,flat=stdflat,Stdout=1)







def extract():
    if os.path.isdir('database'):
        shutil.rmtree('database')

    scilist = glob.glob('AGN*.fits')
    for img in scilist:
        hds = pyfits.getheader(img)
        _gain = hds['HIERARCH ESO DET OUT1 GAIN']
        _ron = hds['HIERARCH ESO DET OUT1 RON']
        name = img.rsplit('.fits',1)[0]+'_EX'
        iraf.specred.apall(img,output=name,interactive='no',find='no',recenter='no',resize='no',edit='no',trace='yes',fittrace='no',extract='yes',extras='yes',review='no',line=600,nsum=20,lower=-5.,upper=5.,b_function='chebyshev',b_order=1,b_sample='-20:-10,10:20',b_naverage=3.,b_niterate=0,b_low_reject=3.,b_high_reject=3.,b_grow=0,width=10.,radius=10.,nfind=1,background='fit',clean='yes',readnoise=_ron,gain=_gain)


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
        iraf.apsum(img,output=name,references=ref,nsum=dsum,interactive='no',find='no',edit='no',trace='no',fittrace='no',extract='yes',extras='yes',review='no',clean='yes',readnoise=_ron,gain=_gain)







def ident():
    lexlist = glob.glob('LAMP*EX.fits')
    linelist = 'Lines_HgCdHeNeAr600.dat'
    width = 10.0
    radius = 15.0
    idref = lexlist[0].rsplit('.fits',1)[0]
    iraf.identify(idref,coordlist=linelist,fwidth=width,cradius=radius)
    lexlist = lexlist[1:]
    for img in lexlist:
        iraf.reidentify(idref,img,cradius=radius)
    return idref


def reid(refpath,idref):
    if not os.path.isdir('database'):
        os.mkdir('database')
    os.system('cp '+refpath+' ./database')
    lexlist = glob.glob('LAMP*EX.fits')
    linelist = 'Lines_HgCdHeNeAr600.dat'
    width = 10.0
    radius = 15.0
    for img in lexlist:
        iraf.reidentify(idref,img,cradius=radius)


def correlate():
    scilist = glob.glob('AGN*EX.fits')
    stdlist = glob.glob('STD*EX.fits')
    lamps = glob.glob('LAMP*EX.fits')

    for img in lamps:
        hdu = pyfits.getheader(img)
        if hdu['HIERARCH ESO DPR TECH'] == 'SPECTRUM':
            specref = img.rsplit('.fits',1)[0]
        elif hdu['HIERARCH ESO DPR TECH'] == 'MOS':
            mosref = img.rsplit('.fits',1)[0]



    for spec in scilist:
        spec = spec.rsplit('.fits',1)[0]
        iraf.refspectra(spec,references=specref,sort='',group='',confirm='no',assign='yes')
        name = spec.rsplit('_EX',1)[0]+'_EL'
        iraf.dispcor(spec,output=name)

    for spec in stdlist:
        spec = spec.rsplit('.fits',1)[0]
        iraf.refspectra(spec,references=mosref,sort='',group='',confirm='no',assign='yes')
        name = spec.rsplit('_EX',1)[0]+'_EL'
        iraf.dispcor(spec,output=name)




        
def standard():
    extinct = 'extinction_lasilla.dat'
    observ = 'paranal'

    if os.path.isfile('std'):
        os.remove('std')
    senslist = glob.glob('sens*.fits')
    for f in senslist:
        os.remove(f)

#    os.system('cp /dark/jsamuel/VLTdata/sens.0001.fits ./')



    stars = {'specphot-LTT7379':'l7379','specphot-LDS749B':'lds749b','specphot-EG274':'eg274','specphot-G138-31':'g13831','specphot-C-32d9927':'cd32','specphot-LTT9491':'l9491','specphot-LTT7987':'l7987'}

    stdlist = glob.glob('STD*EL.fits')

    for spec in stdlist:
        _std = 'std'
        hds = pyfits.getheader(spec)
        _starname = stars[hds['HIERARCH ESO OBS NAME']]
        if _starname in ['l7379','eg274','cd32','l7987','l9491']:
            _caldir = '/dark/jsamuel/agn/standards/stds/'
        else:
            print '%s not found in caldir list' % (_starname)
        _airmass = np.average([hds['HIERARCH ESO TEL AIRM START'],hds['HIERARCH ESO TEL AIRM END']])
        iraf.standard(input=spec,output=_std,extinction=extinct,caldir=_caldir,observatory=observ,star_name=_starname,airmass=_airmass,interact='no')


    iraf.sensfunc(standards='std',sensitivity='sens',extinction=extinct,observatory=observ,interactive='no')





def calibrate():
    scilist = glob.glob('AGN*LF.fits')
    for f in scilist:
        os.remove(f)
    scired = glob.glob('AGN*EL.fits')
    extinct = 'extinction_lasilla.dat'
    observ = 'paranal'
    for spec in scired:
        hdo = pyfits.getheader(spec)
        _airmass = np.average([hdo['HIERARCH ESO TEL AIRM START'],hdo['HIERARCH ESO TEL AIRM END']])
        _exptime = hdo['EXPTIME']
        name = spec.rsplit('_EL',1)[0]+'_LF'
        iraf.calibrate(spec,name,extinction=extinct,observatory=observ,airmass=_airmass,exptime=_exptime)







###MAIN PROGRAM

imgdir = '/dark/jsamuel/scripts/VLT2/'
os.chdir(imgdir)

#iraf.setinstrument(instrument='rca4m',site='kpno',directory='ccddb$',review='no')
'''
imglist = glob.glob('FORS2*.fits')
datelist = []
for img in imglist:
    hdi = pyfits.getheader(img)
    _datetime = hdi['DATE']
    _obsdate = _datetime.rsplit('T',1)[0]
    if _obsdate not in datelist:
        datelist.append(_obsdate)
    if os.path.isdir(_obsdate):
        os.system('mv '+img+' '+_obsdate)
    else:
        os.mkdir(_obsdate)
        os.system('mv '+img+' '+_obsdate)
'''

i = 0
refpath = ''
refpath2 = ''
datelist = os.listdir(imgdir)
datelist.sort()

for _dir in datelist:
    toforget = ['flatcombine','zerocombine','ccdproc','specred.apall','identify','reidentify','specred.standard','dispcor','refspectra','response','apsum','sensfunc','calibrate']
    for t in toforget:
        iraf.unlearn(t)
#    iraf.set(use_new_imt='no')
    iraf.flpr()
    i = i + 1
    print 'Working on data from %s' % (_dir)
    os.chdir(imgdir+_dir)
    imgsort()
    os.chdir('CHIP1')
    os.system('cp /dark/jsamuel/agn/extinction_lasilla.dat ./')
    os.system('cp /dark/jsamuel/agn/Lines_HgCdHeNeAr600.dat ./')
    objdict = {}
    process()
    extract()
    if i <= 1:
        idref = ident()
#        refpath = imgdir+_dir+'/CHIP1/'+idref+'.fits'
        refpath = imgdir+_dir+'/CHIP1/database/id'+idref
        print idref
        print refpath
        print refpath2
    elif i > 1:
        reid(refpath,idref)
    correlate()
    standard()
    calibrate()

