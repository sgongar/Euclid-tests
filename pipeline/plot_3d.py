#!/usr/bin/python
# -*- coding: utf-8 -*-

from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
import numpy as np
from pandas import read_csv


class Plot3D:

    def __init__(self):
        stars_df = read_csv('stars_df.csv', index_col=0)
        self.stars_d = stars_df.to_dict()
        galaxies_df = read_csv('galaxies_df.csv', index_col=0)
        self.galaxies_d = galaxies_df.to_dict()
        ssos_df = read_csv('ssos_df.csv', index_col=0)
        self.ssos_d = ssos_df.to_dict()

        self.data_d = {}
        self.manage_dict()
        self.plot_3d_figure()

    def manage_dict(self):
        """

        :return:
        """
        stars_mag = []
        for key_ in self.stars_d['stars_mag'].keys():
            stars_mag.append(self.stars_d['stars_mag'][key_])
        stars_pm = []
        for key_ in self.stars_d['stars_pm'].keys():
            stars_pm.append(self.stars_d['stars_pm'][key_])
        stars_a = []
        for key_ in self.stars_d['stars_a'].keys():
            stars_a.append(self.stars_d['stars_a'][key_])

        self.data_d['x_stars'] = np.array(stars_mag)
        self.data_d['y_stars'] = np.array(stars_pm)
        self.data_d['z_stars'] = np.array(stars_a)

        # Galaxies
        galaxies_mag = []
        for key_ in self.galaxies_d['galaxies_mag'].keys():
            galaxies_mag.append(self.galaxies_d['galaxies_mag'][key_])
        galaxies_pm = []
        for key_ in self.galaxies_d['galaxies_pm'].keys():
            galaxies_pm.append(self.galaxies_d['galaxies_pm'][key_])
        galaxies_a = []
        for key_ in self.galaxies_d['galaxies_a'].keys():
            galaxies_a.append(self.galaxies_d['galaxies_a'][key_])

        self.data_d['x_galaxies'] = np.array(galaxies_mag)
        self.data_d['y_galaxies'] = np.array(galaxies_pm)
        self.data_d['z_galaxies'] = np.array(galaxies_a)

        # SSOs
        ssos_mag = []
        for key_ in self.ssos_d['ssos_mag'].keys():
            ssos_mag.append(self.ssos_d['ssos_mag'][key_])
        ssos_pm = []
        for key_ in self.ssos_d['ssos_pm'].keys():
            ssos_pm.append(self.ssos_d['ssos_pm'][key_])
        ssos_a = []
        for key_ in self.ssos_d['ssos_b'].keys():
            ssos_a.append(self.ssos_d['ssos_b'][key_])

        self.data_d['x_ssos'] = np.array(ssos_mag)
        self.data_d['y_ssos'] = np.array(ssos_pm)
        self.data_d['z_ssos'] = np.array(ssos_a)

    def plot_3d_figure(self):
        """

        :return:
        """
        plt.figure()
        ax = plt.axes(projection='3d')
        ax.set_xlabel('X - Magnitude')
        ax.set_ylabel('Y - Proper Motion')
        ax.set_ylim3d([0, 5])
        ax.set_zlabel('Z - A size')
        """
        ax.scatter3D(self.data_d['x_stars'], self.data_d['y_stars'],
                     self.data_d['z_stars'], c='r')
        ax.scatter3D(self.data_d['x_galaxies'], self.data_d['y_galaxies'],
                     self.data_d['z_galaxies'], c='b')
        """
        ax.scatter3D(self.data_d['x_ssos'], self.data_d['y_ssos'],
                     self.data_d['z_ssos'], c='g')

        plt.show()


if __name__ == "__main__":
    test = Plot3D()