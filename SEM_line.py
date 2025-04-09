# -*- coding: utf-8 -*-
"""
Created on Tue Mar 26 13:17:13 2024

@author: Rui
"""
import numpy as np

import geopandas as gpd
from shapely import geometry
import matplotlib.pyplot as plt
import yaml

class SEM_Line:

    def __init__(self, cfg_file, **kwargs):
        
        with open(cfg_file, 'r') as file:
            # dictionary defined in configuration file
            config_dict = yaml.safe_load(file)        
    
        for key, value in config_dict.items():
            setattr(self, key, value)
        
        for key, value in kwargs.items():
            setattr(self, key, value)

        if self.filename == []:
            return
        try:
            self.original_line = gpd.read_file(self.filename)
        except IOError:
            print('An error occured trying to read the file  -> ', self.filename)

    def compute_mid_cells(self):
        self.xc = (self.x[1:] + self.x[:-1]) / 2
        self.yc = (self.y[1:] + self.y[:-1]) / 2

    def build_transects(self, transect_length):
        self.trans_dir = self.normals_from_cell_boundaries()
        self.x_trans, self.y_trans = self.xycoords_from_distance(transect_length)

        lines = []
        for xc, yc, x_trans, y_trans in zip(self.xc, self.yc, self.x_trans, self.y_trans):
            lines.append(geometry.LineString(((xc, yc), (x_trans, y_trans))))
        self.transects = gpd.GeoDataFrame(crs=self.original_line.crs, geometry=lines)

    def normals_from_cell_boundaries(self):
        dx = self.x[1:] - self.x[:-1]
        dy = self.y[1:] - self.y[:-1]
        return np.pi / 2 + np.arctan2(dy, dx)

    def normals_from_cell_centers(self):
        dx = self.xc[1:] - self.xc[:-1]
        dy = self.yc[1:] - self.yc[:-1]
        normals = np.pi / 2 + np.arctan2(dy, dx)

        return np.insert(normals, (0, -1), (normals[0], normals[-1]),
                         axis=0)  # first and last cells have the same azimuth

    def azimuth_from_normals_degress(self):  # azimuth from shore_normals
        normals = self.normals_from_cell_centers()
        return np.mod(90 - np.rad2deg(normals), 360)

    def cell_lenght(self):
        dx = self.x[1:] - self.x[:-1]
        dy = self.y[1:] - self.y[:-1]
        return np.hypot(dx, dy)

    def create_grid(self, dx,
                    transect_length):  # adpated from https://stackoverflow.com/questions/34906124/interpolating-every-x-distance-along-multiline-in-shapely/35025274#35025274
        geom = self.original_line.geometry
        self.num_cells = int(round(geom.length / dx))
        points = [geom.interpolate(float(n) / self.num_cells, normalized = True) for n in range(self.num_cells + 1)]

        self.x = np.array([p.x for p in points])
        self.y = np.array([p.y for p in points])

        self.compute_mid_cells()
        self.build_transects(transect_length)

    def xycoords_from_distance(self, dist):
        xc = self.xc + dist * np.cos(self.trans_dir)
        yc = self.yc + dist * np.sin(self.trans_dir)
        return xc, yc

    def plot_shp(self):
        xy = self.original_line.geometry[0].coords.xy

        plt.plot(xy[0], xy[1], color=self.shp_color, markersize=self.shp_markersize,
                 marker=self.shp_marker, markeredgecolor=self.shp_markeredgecolor, linestyle=self.shp_linestyle,
                 linewidth=self.shp_linewidth)

    def plot_x(self):
        plt.plot(self.x, self.y, color=self.x_color, markersize=self.x_markersize,
                 marker=self.x_marker, markeredgecolor=self.x_markeredgecolor, linestyle=self.x_linestyle,
                 linewidth=self.x_linewidth)

    def plot_xc(self):
        plt.plot(self.xc, self.yc, color=self.xc_color, markersize=self.xc_markersize,
                 marker=self.xc_marker, markeredgecolor=self.xc_markeredgecolor, linestyle=self.xc_linestyle,
                 linewidth=self.xc_linewidth)

    def plot(self):
        if self.shp_plot:
            self.plot_shp()
        if self.xc_plot:
            self.plot_xc()
        if self.x_plot:
            self.plot_x()

