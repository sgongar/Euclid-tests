#!/usr/bin/python
# -*- coding: utf-8 -*-

"""Python script for time measurements

Todo:
    * Improve log messages

"""

from pandas import concat, read_csv, Series
from pipeline.misc import compare_floats

__author__ = "Samuel Gongora-Garcia"
__copyright__ = "Copyright 2017"
__credits__ = ["Samuel Gongora-Garcia"]
__license__ = "GPL"
__version__ = "0.1"
__maintainer__ = "Samuel Gongora-Garcia"
__email__ = "sgongora@cab.inta-csic.es"
__status__ = "Development"



def analyse_perfomance(logger, prefs_dict):
    """

    @param logger:
    @param prefs_dict:

    """
    list_analysis = []
    list_total_detections = []
    list_total_sources = []
    list_ok_sources = []
    list_no_sources = []
    porcentages = []
    list_f_dr = []
    list_f_pur = []
    list_f_com = []

    sources_cat = read_csv('regions.reg', names=['alpha', 'delta'],
                           delimiter=" ")
    confs_cat = read_csv('configurations_2.csv',
                         names=['deblending', 'mincount',
                                'threshold', 'minarea'], delimiter=",")

    for csv_num in range(0, 54, 1):
        cat_file = read_csv('reports/full_{}.csv'.format(csv_num), index_col=0)

        cats = {}
        # read the four catalogues from images
        for cat in range(0, 4, 1):
            cats['{}'.format(cat)] = cat_file[cat_file['source_cat'].isin([cat])]
            """
            cats['{}'.format(cat)].to_csv('reports/out/out_{}_{}'.format(csv_num,
                                                                         cat))
            """
        sources = 180
        total_detections = 0
        total_sources = 0
        ok_sources = 0
        no_sources = 0
        alpha_tolerance = 0.0006  # 1 arsec tolerance
        delta_tolerance = 0.0006  # 1 arsec tolerance

        # looping over catalogues read
        for cat in range(0, 4, 1):
            catalogue = cats['{}'.format(cat)]
            list_sources = catalogue['source_num'].tolist()
            total_detections = total_detections + len(list_sources)
            unique_sources = sorted(list(set(list_sources)))
            total_sources = total_sources + len(unique_sources)

            # check each source
            for source_num in range(catalogue['source_num'].size):
                # getting alpha, delta for each source sextracted
                source_alpha = catalogue['source_alpha'].iloc[source_num]
                source_delta = catalogue['source_delta'].iloc[source_num]

                source = False

                for orig_source in range(sources_cat['alpha'].size):
                    # getting original values
                    original_alpha = sources_cat['alpha'].iloc[orig_source]
                    original_delta = sources_cat['delta'].iloc[orig_source]
                    # comparing values
                    alpha_comp = compare_floats(original_alpha, source_alpha,
                                                alpha_tolerance)
                    delta_comp = compare_floats(original_delta, source_delta,
                                                delta_tolerance)

                    # checking if source is real
                    if alpha_comp and delta_comp:
                        ok_sources += 1
                        source = True
                if not source:
                    no_sources += 1

        ok_sources = float(ok_sources)
        no_sources = float(no_sources)
        sources = float(sources)

        f_dr = (ok_sources + no_sources) / sources
        f_pur = ok_sources / (ok_sources + no_sources)
        f_com = f_dr * f_pur

        porcentage = (float(no_sources) / 180.0) * 100

        list_analysis.append(csv_num)
        list_total_detections.append(total_detections)
        list_total_sources.append(total_sources)
        list_ok_sources.append(ok_sources)
        list_no_sources.append(no_sources)
        porcentages.append(porcentage)
        list_f_dr.append(f_dr)
        list_f_pur.append(f_pur)
        list_f_com.append(f_com)

    # final table with all analyses
    t_1 = Series(list_analysis, name='analysis')
    t_2 = confs_cat['deblending']
    t_3 = confs_cat['mincount']
    t_4 = confs_cat['threshold']
    t_5 = confs_cat['minarea']
    t_6 = Series(list_total_detections, name='total_detections')
    t_7 = Series(list_total_sources, name='unique_sources')
    t_8 = Series(list_ok_sources, name='ok_sources')
    t_9 = Series(list_no_sources, name='no_sources')
    t_10 = Series(porcentages, name='fake_sources')
    t_11 = Series(list_f_dr, name='detection_rate')
    t_12 = Series(list_f_pur, name='purity')
    t_13 = Series(list_f_com, name='completeness')

    test_table = concat([t_1, t_2, t_3, t_4, t_5, t_6,
                         t_7, t_8, t_9, t_10, t_11, t_12, t_13], axis=1)
    test_table.to_csv('test_table.csv')
