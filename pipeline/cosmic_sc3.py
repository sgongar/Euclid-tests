# !/usr/bin/python
# -*- coding: utf-8 -*-

"""Python script for cosmic rays removal process

Versions:
- 0.1

Todo:
    * Improve log messages

"""

from os import listdir
from astropy.io import fits
from astroscrappy import detect_cosmics
from multiprocessing import Process


def cosmic():
    fits_dir = '/pcdisk/holly/sgongora/Dev/CCDs'
    output_dir = '/pcdisk/holly/sgongora/Dev/CCDs_filt'
    cores_number = 8

    active_cosmic = []
    fits_files = listdir(fits_dir)

    for cosmic_idx in range(0, len(fits_files), cores_number):
        try:
            cosmic_j = []
            for proc in range(0, cores_number, 1):
                idx = cosmic_idx + proc  # index
                print(idx, len(fits_files))
                fits_file = fits_files[idx]
                cosmic_p = Process(target=cosmic_thread,
                                   args=(fits_dir, output_dir, fits_file,))
                cosmic_j.append(cosmic_p)
                cosmic_p.start()

                active_cosmic = list([job.is_alive() for job in cosmic_j])
            while True in active_cosmic:
                active_cosmic = list([job.is_alive() for job in cosmic_j])
                pass
        except IndexError:
            print('Extraction finished')

    print('Extraction process of fits images finished')

    return True


def cosmic_thread(fits_dir, output_dir, fits_file):
    """
    """
    data_, header = fits.getdata('{}/{}'.format(fits_dir, fits_file),
                                 header=True)

    (cr_mask,
     cleaned_array) = detect_cosmics(data_, sigclip=4.5, sigfrac=0.3,
                                     objlim=5.0, gain=3.1, readnoise=4,
                                     satlevel=64535.0, pssl=0.0, niter=4,
                                     sepmed=True, cleantype='meanmask',
                                     fsmode='median', psfmodel='gauss',
                                     psffwhm=0.18, psfsize=7, psfk=None,
                                     psfbeta=4.765, verbose=True)

    fits.writeto('{}/{}'.format(output_dir, fits_file), cleaned_array, header)
