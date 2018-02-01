#!/usr/bin/python
# -*- coding: utf-8 -*-

"""Plots pm and magnitude against A/B/elongation
Galaxies data

Versions:
- 0.1 Initial release

Todo:
    * Improve usability

"""
import matplotlib.pyplot as plt
import numpy as np
from pandas import read_csv


class PlotMagVSPMStars:

    def __init__(self):
        stars_df = read_csv('stars_df.csv', index_col=0)
        self.stars_d = stars_df.to_dict()

        self.data_d = {}
        self.manage_stars_dict()
        self.plot_figure()

    def manage_stars_dict(self):
        """

        :return:
        """
        # Stars
        stars_mag = []
        for key_ in self.stars_d['stars_mag'].keys():
            stars_mag.append(self.stars_d['stars_mag'][key_])
        stars_pm = []
        for key_ in self.stars_d['stars_pm'].keys():
            stars_pm.append(self.stars_d['stars_pm'][key_])

        self.data_d['x_stars'] = np.array(stars_mag)
        self.data_d['y_stars'] = np.array(stars_pm)

    def plot_figure(self):
        """

        :return:
        """
        fig = plt.figure(figsize=(16.53, 11.69), dpi=100)
        ax = fig.add_subplot(1, 1, 1)

        # ax.semilogy(self.data_d['x_galaxies'], self.data_d['y_galaxies'],
        #             'bs', ms=1)
        ax.plot(self.data_d['x_stars'], self.data_d['y_stars'], 'bs', ms=1)
        ax.set_ylim([0.001, 1])
        # ax.grid(b=True, which='major', ls='-', lw=2)
        # ax.grid(b=True, which='minor', ls='--', lw=1)

        plt.grid(True)
        plt.show()


class PlotMagVSPMGalaxies:

    def __init__(self):
        galaxies_df = read_csv('galaxies_df.csv', index_col=0)
        self.galaxies_d = galaxies_df.to_dict()

        self.data_d = {}
        self.manage_galaxies_dict()
        self.plot_figure()

    def manage_galaxies_dict(self):
        """

        :return:
        """
        # Galaxies
        galaxies_mag = []
        for key_ in self.galaxies_d['galaxies_mag'].keys():
            galaxies_mag.append(self.galaxies_d['galaxies_mag'][key_])
        galaxies_pm = []
        for key_ in self.galaxies_d['galaxies_pm'].keys():
            galaxies_pm.append(self.galaxies_d['galaxies_pm'][key_])

        self.data_d['x_galaxies'] = np.array(galaxies_mag)
        self.data_d['y_galaxies'] = np.array(galaxies_pm)

    def plot_figure(self):
        """

        :return:
        """
        fig = plt.figure(figsize=(16.53, 11.69), dpi=100)
        ax = fig.add_subplot(1, 1, 1)

        ax.semilogy(self.data_d['x_galaxies'], self.data_d['y_galaxies'],
                    'bs', ms=1)
        # ax.plot(self.data_d['x_galaxies'], self.data_d['y_galaxies'], 'bs',
        #         ms=1)
        ax.set_ylim([0.001, 1])
        # ax.grid(b=True, which='major', ls='-', lw=2)
        # ax.grid(b=True, which='minor', ls='--', lw=1)

        plt.grid(True)
        plt.show()


if __name__ == "__main__":
    test = PlotMagVSPMGalaxies()
