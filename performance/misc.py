#!/usr/bin/python
# -*- coding: utf-8 -*-

"""

Content:
    * get_fits_limits
    * get_fits_d
    * get_os
    * conf_map
    * extract_settings_luca
    * extract_settings_sc3

Todo:
    * Improve log messages
    * Improve usability
"""
from multiprocessing import cpu_count
from ConfigParser import ConfigParser
from os import listdir
from platform import platform
from logging import getLogger, config

from astropy.io import fits
from astropy.coordinates import SkyCoord
from astropy.units import degree
from astropy.wcs import WCS

from errors import BadSettings

__author__ = "Samuel Góngora García"
__copyright__ = "Copyright 2018"
__credits__ = ["Samuel Góngora García"]
__version__ = "0.1"
__maintainer__ = "Samuel Góngora García"
__email__ = "sgongora@cab.inta-csic.es"
__status__ = "Development"


def create_sextractor_dict(conf_num, cat_conf):
    """

    @param conf_num:
    @param cat_conf:

    @return analysis_d:
    """

    if cat_conf:
        configurations = [2, 0.1, 5, 4, 'models/gauss_2.0_5x5.conv']
        len_conf = 1
    else:
        mode = {'type': 'sextractor'}  # harcoded
        configurations, len_conf = create_configurations(mode)

    analysis_l = []

    if type(configurations[0]) is list:
        for configuration in range(len(configurations)):
            temp_list = [configurations[configuration][0],
                         configurations[configuration][1],
                         configurations[configuration][2],
                         configurations[configuration][2],
                         configurations[configuration][3],
                         configurations[configuration][4]]
            analysis_l.append(temp_list)
        analysis_d = {'deblend_nthresh': analysis_l[conf_num][0],
                      'deblend_mincount': analysis_l[conf_num][1],
                      'detect_thresh': analysis_l[conf_num][2],
                      'analysis_thresh': analysis_l[conf_num][3],
                      'detect_minarea': analysis_l[conf_num][4],
                      'filter': analysis_l[conf_num][5]}
    else:
        analysis_d = {'deblend_nthresh': configurations[0],
                      'deblend_mincount': configurations[1],
                      'detect_thresh': configurations[2],
                      'analysis_thresh': configurations[2],
                      'detect_minarea': configurations[3],
                      'filter': configurations[4]}

    return analysis_d, len_conf


def create_configurations(mode):
    """ creates a list of configuration lists merging
        different input parameters
    :param mode: can be 'sextractor' or 'scamp'
    :return:
    """
    if mode['type'] == 'sextractor':
        l_deblending = [30]
        l_mincount = [0.01]
        l_threshold = [1.5]

        l_area = [4]
        l_filter_name = ['models/gauss_2.0_5x5.conv']

        configurations = []
        for deblending in l_deblending:
            for mincount in l_mincount:
                for threshold in l_threshold:
                    for area in l_area:
                        for filt in l_filter_name:
                            configurations.append([deblending, mincount,
                                                   threshold, area, filt])
        configurations_len = len(configurations)
        return configurations, configurations_len

    elif mode['type'] == 'scamp':
        l_crossid_radius = [10]  # [10] seconds
        l_pixscale_maxerr = [1.1]  # [1.2] scale-factor
        l_posangle_maxerr = [0.5]  # [0.5, 2.5] degrees
        l_position_maxerr = [0.04]

        configurations = []

        for crossid in l_crossid_radius:
            for pixscale in l_pixscale_maxerr:
                for posangle in l_posangle_maxerr:
                    for position in l_position_maxerr:
                        configurations.append([crossid, pixscale,
                                               posangle, position])

        configurations_len = len(configurations)

        return configurations, configurations_len


def create_scamp_dict(conf_num):
    """

    :param conf_num:
    :return:
    """

    scamp_list = []
    mode = {'type': 'scamp'}
    configurations, len_conf = create_configurations(mode)

    for conf in configurations:
        temp_list = [conf[0], conf[1], conf[2], conf[3]]
        scamp_list.append(temp_list)
    scamp_dict = {'crossid_radius': scamp_list[conf_num][0],
                  'pixscale_maxerr': scamp_list[conf_num][1],
                  'posangle_maxerr': scamp_list[conf_num][2],
                  'position_maxerr': scamp_list[conf_num][3]}

    return scamp_dict, len_conf


def setting_logger(prfs_d, logger_name):
    """ sets-up a logger object ready to be used

    TODO improve logger definition

    @return logger:
    """
    config.fileConfig(prfs_d['logger_config'])

    # TODO implement logger level setting
    """
    if argv[4] == '-INFO':
        logger = getLogger("main_process").setLevel(INFO)
    elif argv[4] == '-DEBUG':
        logger = getLogger("main_process").setLevel(DEBUG)
    else:
        raise Exception
    """
    logger = getLogger(logger_name)
    logger.info("Pipeline started")
    # logger.FileHandler('spam.log')

    return logger


def get_cats(dither):
    """

    :return:
    """
    prfs_d = extract_settings_elvis()
    cat_list = []

    files = listdir('{}/'.format(prfs_d['fits_dir']))
    for file_ in files:
        if file_[-4:] == '.cat':
            cat_list.append(file_)

    cats_out = []
    for file_ in cat_list:
        if file_[-5:-4] == str(dither):
            cats_out.append(file_)

    return cats_out


def get_cat(ccd):
    """ returns catalog from ccd name

    :param ccd:
    :return:
    """
    cats = []
    idx = 0

    for x_ in range(1, 7, 1):
        for y_ in range(1, 7, 1):
            for d_ in range(1, 5, 1):
                cat_name = 'x{}_y{}'.format(x_, y_, d_)
                cats.append([cat_name, d_, idx])

                idx += 1

    ccd_position = ccd[4:9]
    dither_n = ccd[11:12]

    for cat_ in cats:
        if str(ccd_position) == cat_[0] and int(dither_n) == cat_[1]:
            cat_n = cat_[2]

    return cat_n


def get_fits(dither):
    """

    :return:
    """
    prfs_d = extract_settings_elvis()
    fits_list = []

    files = listdir('{}/'.format(prfs_d['fits_dir']))
    for file_ in files:
        if file_[-5:] == '.fits':
            fits_list.append(file_)

    fits_out = []
    for file_ in fits_list:
        if file_[-6:-5] == str(dither):
            fits_out.append(file_)

    return fits_out


def get_fits_limits(fits_image):
    """ todo - to another new file?

    @param fits_image: fits image

    @return limits: a dict with ra/dec limits above_ra, below_ra,
                    above_dec_, below_dec
    """
    # logger.info('getting limits of {} image'.format(fits_image))

    data, header = fits.getdata(fits_image, header=True)
    w = WCS(fits_image)

    above_x, above_y = header['NAXIS1'], header['NAXIS2']
    above_ra, above_dec = w.all_pix2world(above_x, above_y, 0)

    below_ra, below_dec = w.all_pix2world(0, 0, 0)

    c = SkyCoord(ra=[above_ra, below_ra] * degree,
                 dec=[above_dec, below_dec] * degree)

    ra = [above_ra, below_ra]
    dec = [above_dec, below_dec]

    limits = {'below_ra': float(min(ra)), 'above_ra': float(max(ra)),
              'below_dec': float(min(dec)), 'above_dec': float(max(dec))}

    # check position
    # sometimes some values could be higher when are tagged as "lowest"
    return limits


def get_cats_elvis_d(dither):
    """

    :param mag_:
    :param dither:
    :return:
    """
    prfs_d = extract_settings_elvis()
    fits_list = []

    files = listdir('{}'.format(prfs_d['fits_dir']))
    for file_ in files:
        if file_[-4:] == '.cat':
            fits_list.append(file_)

    list_out = []
    for file_ in fits_list:
        if file_[-5:-4] == str(dither):
            list_out.append(file_)

    return list_out


def get_fits_d(mag_, dither):
    """

    :param mag_:
    :param dither:
    :return:
    """
    prfs_d = extract_settings_luca()
    fits_list = []

    files = listdir('{}/{}/CCDs/'.format(prfs_d['fits_dir'], mag_))
    for file_ in files:
        if file_[:1] == 'm' and file_[-5:] == '.fits':
            fits_list.append(file_)

    list_out = []
    for file_ in fits_list:
        if file_[-6:-5] == str(dither):
            list_out.append(file_)

    return list_out


def get_os():
    """ a function that gets the current operative system
    for now works in Debian, Fedora and Ubuntu shell (Microsoft version)

    @return os_system: a string which contains the operative system name
    """
    if 'fedora-19' in platform() or 'fedora-23' in platform():
        os_system = 'cab'
    elif 'centos' in platform():
        os_system = 'centos'
    else:
        raise Exception

    return os_system


def conf_map(config_, section):
    """

    @param config_:
    @param section:

    @return dict1:
    """
    dict1 = {}
    options = config_.options(section)
    for option in options:
        try:
            dict1[option] = config_.get(section, option)
            if dict1[option] == -1:
                print('skip: {}'.format(option))
        except KeyError:
            print('exception on {}'.format(option))
            dict1[option] = None
    return dict1


def extract_settings_luca():
    """ creates a dictionary with all the configuration parameters
        at this moment configuration file location is fixed at main directory

    @return prfs_d: a dictionary which contains all valuable data
    """
    cf = ConfigParser()
    cf.read(".settings.ini")

    prfs_d = {}
    os_version = get_os()

    if os_version == 'centos':
        prfs_d['version'] = conf_map(cf, "Version")['centos_version']
    elif os_version == 'cab':
        prfs_d['version'] = conf_map(cf, "Version")['cab_version']
    else:
        raise BadSettings('Operative system not chosen')

    if os_version == 'centos':
        prfs_d['home'] = conf_map(cf, "HomeDirs")['centos_home']
    elif os_version == 'cab':
        prfs_d['home'] = conf_map(cf, "HomeDirs")['cab_home']
    else:
        raise BadSettings('Operative system not chosen')

    prfs_d['fits_dir'] = conf_map(cf, "ImagesDirs")['fits_dir']
    prfs_d['fits_dir'] = '{}{}'.format(prfs_d['version'], prfs_d['fits_dir'])

    # todo - comment!
    prfs_d['output_cats'] = conf_map(cf, "CatsDirs")['output_cats']
    prfs_d['output_cats'] = prfs_d['version'] + prfs_d['output_cats']
    # todo - comment!
    prfs_d['references'] = conf_map(cf, "CatsDirs")['references']
    prfs_d['references'] = prfs_d['version'] + prfs_d['references']
    # todo - comment!
    prfs_d['filtered'] = conf_map(cf, "CatsDirs")['filtered']
    prfs_d['filtered'] = prfs_d['version'] + prfs_d['filtered']

    outputdirs_list = ['conf_scamp', 'conf_sex', 'params_sex', 'neural_sex',
                       'params_cat', 'logger_config']
    for conf_ in outputdirs_list:
        prfs_d[conf_] = conf_map(cf, "ConfigDirs")[conf_]
        prfs_d[conf_] = prfs_d['home'] + prfs_d[conf_]

    prfs_d['first_star'] = conf_map(cf, "CatsOrganization")['first_star']
    prfs_d['first_star'] = int(prfs_d['first_star'])
    prfs_d['first_galaxy'] = conf_map(cf, "CatsOrganization")['first_galaxy']
    prfs_d['first_galaxy'] = int(prfs_d['first_galaxy'])
    prfs_d['first_sso'] = conf_map(cf, "CatsOrganization")['first_sso']
    prfs_d['first_sso'] = int(prfs_d['first_sso'])

    prfs_d['detections'] = int(conf_map(cf, "Misc")['detections'])
    prfs_d['pm_low'] = float(conf_map(cf, "Misc")['pm_low'])
    prfs_d['pm_up'] = float(conf_map(cf, "Misc")['pm_up'])
    pms = conf_map(cf, "Misc")['pms']
    pms = pms.replace(",", " ")
    prfs_d['pms'] = [float(x) for x in pms.split()]
    mags = conf_map(cf, "Misc")['mags']
    mags = mags.replace(",", " ")
    prfs_d['mags'] = mags.split()
    prfs_d['r_fit'] = conf_map(cf, "Misc")['r_fit']
    prfs_d['cores_number'] = conf_map(cf, "Misc")['cores_number']
    if prfs_d['cores_number'] == '0':
        prfs_d['cores_number'] = int(str(cpu_count()))
        # TODO should leave free at least 20% of processors
    else:
        prfs_d['cores_number'] = int(prfs_d['cores_number'])
    prfs_d['tolerance'] = float(conf_map(cf, "Misc")['tolerance'])

    return prfs_d


def extract_settings_sc3():
    """ creates a dictionary with all the configuration parameters
        at this moment configuration file location is fixed at main directory

    @return prfs_d: a dictionary which contains all valuable data
    """
    cf = ConfigParser()
    cf.read(".settings_SC3.ini")

    prfs_d = {}
    os_version = get_os()

    if os_version == 'centos':
        prfs_d['version'] = conf_map(cf, "Version")['centos_version']
    elif os_version == 'cab':
        prfs_d['version'] = conf_map(cf, "Version")['cab_version']
    else:
        raise BadSettings('Operative system not chosen')

    if os_version == 'centos':
        prfs_d['home'] = conf_map(cf, "HomeDirs")['centos_home']
    elif os_version == 'cab':
        prfs_d['home'] = conf_map(cf, "HomeDirs")['cab_home']
    else:
        raise BadSettings('Operative system not chosen')

    prfs_d['fits_dir'] = conf_map(cf, "ImagesDirs")['fits_dir']
    prfs_d['fits_dir'] = '{}{}'.format(prfs_d['version'], prfs_d['fits_dir'])

    # todo - comment!
    prfs_d['output_cats'] = conf_map(cf, "CatsDirs")['output_cats']
    prfs_d['output_cats'] = prfs_d['version'] + prfs_d['output_cats']
    # todo - comment!
    prfs_d['references'] = conf_map(cf, "CatsDirs")['references']
    prfs_d['references'] = prfs_d['version'] + prfs_d['references']
    # todo - comment!
    prfs_d['filtered'] = conf_map(cf, "CatsDirs")['filtered']
    prfs_d['filtered'] = prfs_d['version'] + prfs_d['filtered']

    outputdirs_list = ['conf_scamp', 'conf_sex', 'params_sex', 'neural_sex',
                       'params_cat', 'logger_config']
    for conf_ in outputdirs_list:
        prfs_d[conf_] = conf_map(cf, "ConfigDirs")[conf_]
        prfs_d[conf_] = prfs_d['home'] + prfs_d[conf_]

    prfs_d['detections'] = int(conf_map(cf, "Misc")['detections'])
    prfs_d['pm_low'] = float(conf_map(cf, "Misc")['pm_low'])
    prfs_d['pm_up'] = float(conf_map(cf, "Misc")['pm_up'])
    prfs_d['r_fit'] = conf_map(cf, "Misc")['r_fit']
    prfs_d['cores_number'] = conf_map(cf, "Misc")['cores_number']
    if prfs_d['cores_number'] == '0':
        prfs_d['cores_number'] = int(str(cpu_count()))
        # TODO should leave free at least 20% of processors
    else:
        prfs_d['cores_number'] = int(prfs_d['cores_number'])
    prfs_d['tolerance'] = float(conf_map(cf, "Misc")['tolerance'])

    return prfs_d


def extract_settings_elvis():
    """ creates a dictionary with all the configuration parameters
        at this moment configuration file location is fixed at main directory

    @return prfs_d: a dictionary which contains all valuable data
    """
    cf = ConfigParser()
    cf.read(".settings_ELViS.ini")

    prfs_d = {}
    os_version = get_os()

    if os_version == 'centos':
        prfs_d['version'] = conf_map(cf, "Version")['centos_version']
    elif os_version == 'cab':
        prfs_d['version'] = conf_map(cf, "Version")['cab_version']
    else:
        raise BadSettings('Operative system not chosen')

    if os_version == 'centos':
        prfs_d['home'] = conf_map(cf, "HomeDirs")['centos_home']
    elif os_version == 'cab':
        prfs_d['home'] = conf_map(cf, "HomeDirs")['cab_home']
    else:
        raise BadSettings('Operative system not chosen')

    prfs_d['fits_dir'] = conf_map(cf, "ImagesDirs")['fits_dir']
    prfs_d['fits_dir'] = '{}{}'.format(prfs_d['version'], prfs_d['fits_dir'])
    prfs_d['fpas_dir'] = conf_map(cf, "ImagesDirs")['fpas_dir']
    prfs_d['fpas_dir'] = '{}{}'.format(prfs_d['version'], prfs_d['fpas_dir'])

    # todo - comment!
    prfs_d['output_cats'] = conf_map(cf, "CatsDirs")['output_cats']
    prfs_d['output_cats'] = prfs_d['version'] + prfs_d['output_cats']
    # todo - comment!
    prfs_d['references'] = conf_map(cf, "CatsDirs")['references']
    prfs_d['references'] = prfs_d['version'] + prfs_d['references']
    # todo - comment!
    prfs_d['filtered'] = conf_map(cf, "CatsDirs")['filtered']
    prfs_d['filtered'] = prfs_d['version'] + prfs_d['filtered']

    prfs_d['time_1'] = conf_map(cf, "ImagesTime")['time_1']  # 1st dither time
    prfs_d['time_2'] = conf_map(cf, "ImagesTime")['time_2']  # 2nd dither time
    prfs_d['time_3'] = conf_map(cf, "ImagesTime")['time_3']  # 3nd dither time
    prfs_d['time_4'] = conf_map(cf, "ImagesTime")['time_4']  # 4th dither time

    outputdirs_list = ['conf_scamp', 'conf_sex', 'params_sex', 'neural_sex',
                       'params_cat', 'logger_config']
    for conf_ in outputdirs_list:
        prfs_d[conf_] = conf_map(cf, "ConfigDirs")[conf_]
        prfs_d[conf_] = prfs_d['home'] + prfs_d[conf_]

    prfs_d['detections'] = int(conf_map(cf, "Misc")['detections'])
    prfs_d['pm_low'] = float(conf_map(cf, "Misc")['pm_low'])
    prfs_d['pm_up'] = float(conf_map(cf, "Misc")['pm_up'])
    prfs_d['pm_sn'] = float(conf_map(cf, "Misc")['pm_sn'])
    pms = conf_map(cf, "Misc")['pms']
    pms = pms.replace(",", " ")
    prfs_d['pms'] = [float(x) for x in pms.split()]
    prfs_d['r_fit'] = conf_map(cf, "Misc")['r_fit']
    prfs_d['cores_number'] = conf_map(cf, "Misc")['cores_number']
    if prfs_d['cores_number'] == '0':
        prfs_d['cores_number'] = int(str(cpu_count()))
        # TODO should leave free at least 20% of processors
    else:
        prfs_d['cores_number'] = int(prfs_d['cores_number'])
    prfs_d['tolerance'] = float(conf_map(cf, "Misc")['tolerance'])

    return prfs_d
