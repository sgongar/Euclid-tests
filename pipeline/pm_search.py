#!/usr/bin/python
# -*- coding: utf-8 -*-

from pandas import read_csv

source_n = input('Catalog number: ')

f_file = read_csv('/home/sgongora/Documents/CarpetaCompartida/filt_10_1.2_2.5_0.64_20-21_4.csv', index_col=0)

f_pm = f_file[f_file['SOURCE_NUMBER'].isin([source_n])]
f_pm = f_pm['PM'].iloc[0]

print f_pm

