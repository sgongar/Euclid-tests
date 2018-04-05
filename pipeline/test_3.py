from astropy.io import fits
from astropy.time import Time




time_1 = '2021-06-26T09:00:00.00000'
fits_file = 'EUC_VIS_SWL-DET-001-000000-0000000__20170630T011437.3Z_00.00.fits'
for idx in range(0, 36, 1):
    print(idx)
    """
    fits_name = 'EUC_VIS_SWL-DET-001-000000-0000000__20170630T011437.3Z_00.00_{}.fits'.format(idx)
    data, header = fits.getdata(fits_name, header=True)
    header['DATE-OBS'] = time_1
    # Change format
    t = Time(time_1)
    t.format = 'mjd'
    header['MJD-OBS'] = float(str(t))

    fits.writeto('{}_copy'.format(fits_name), data, header)

    fits_name = 'EUC_VIS_SWL-DET-001-000000-0000000__20170630T011437.3Z_00.00_f{}.fits'.format(idx)
    data, header = fits.getdata(fits_name, header=True)
    header['DATE-OBS'] = time_1
    # Change format
    t = Time(time_1)
    t.format = 'mjd'
    header['MJD-OBS'] = float(str(t))

    fits.writeto('{}_copy'.format(fits_name), data, header)
    """

    fits_name = 'EUC_VIS_SWL-DET-001-000000-0000000__20170630T011437.3Z_00.00_w{}.fits'.format(idx)
    data, header = fits.getdata(fits_name, header=True)
    header['DATE-OBS'] = time_1
    # Change format
    t = Time(time_1)
    t.format = 'mjd'
    header['MJD-OBS'] = float(str(t))

    fits.writeto('{}_copy'.format(fits_name), data, header)