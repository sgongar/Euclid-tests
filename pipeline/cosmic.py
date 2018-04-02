# !/usr/bin/python
# -*- coding: utf-8 -*-

"""Python script for cosmic rays removal process

Versions:
- 0.1

Todo:
    * Improve log messages

"""

from astropy.io import fits
from astroscrappy import detect_cosmics

dir_ = '/home/sgongora/Documents/CarpetaCompartida/luca_data/VIS_SC3/CCDs'
# file_ = 'EUC_VIS_SWL-DET-001-000000-0000000__20170630T011437.3Z_00.00'
# file_ = 'EUC_VIS_SWL-DET-002-000000-0000000__20170630T011642.0Z_00.00'
# file_ = 'EUC_VIS_SWL-DET-003-000000-0000000__20170630T011848.6Z_00.00'
file_ = 'EUC_VIS_SWL-DET-004-000000-0000000__20170630T012050.1Z_00.00'

indat, header = fits.getdata('{}/{}.fits'.format(dir_, file_), header=True)

crmask_, cleanarr = detect_cosmics(indat, sigclip=4.5, sigfrac=0.3, objlim=5.0,
                                   gain=3.1, readnoise=4, satlevel=64535.0,
                                   pssl=0.0, niter=4, sepmed=True,
                                   cleantype='meanmask', fsmode='median',
                                   psfmodel='gauss', psffwhm=0.18,
                                   psfsize=7, psfk=None, psfbeta=4.765,
                                   verbose=True)

fits.writeto('{}/{}_copy.fits'.format(dir_, file_), cleanarr, header)
