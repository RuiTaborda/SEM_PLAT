# -*- coding: utf-8 -*-
"""
Created on Tue Mar 12 14:59:16 2024

@author: Rui
"""
from shapely import geometry
import matplotlib.pyplot as plt
import yaml
import numpy as np


class BeachProfile:
    def __init__(self, cfg_file, **kwargs):
        
        with open(cfg_file, 'r') as file:
            # dictionary defined in configuration file
            config_dict = yaml.safe_load(file)        
    
        for key, value in config_dict.items():
            setattr(self, key, value)
   
        for key, value in kwargs.items():
            setattr(self, key, value)

        self.build_rocky_profile()
        self.build_beach_profile()

    
    def build_beach_profile(self):
        self.y_shoreline = self.berm_with
        self.z_shoreline = self.z_coastline - self.berm_with * self.berm_slope


        # Dean profile starts where the slope matches the beach face slope
        self.distance_beach_face_slope = 8 * self.dean_parameter **3 / (27 * self.beach_face_slope ** 3)
        
        self.z_beach_face_toe = self.z_shoreline - self.distance_beach_face_slope * self.beach_face_slope
        coords = np.array([(0, self.z_coastline), (self.y_shoreline, self.z_shoreline), 
                  (self.y_shoreline + self.distance_beach_face_slope,  self.z_beach_face_toe)])
        
        coords_dean = self.build_dean_profile()
        coords = np.concatenate((coords, coords_dean))
        self.beach_polygon = geometry.Polygon(coords)
        
    def build_dean_profile(self):
        y_min = self.distance_beach_face_slope
        y_max = (self.closure_depth / self.dean_parameter) ** (3/2)
        y_dean = np.arange(y_min, y_max, self.y_resolution)
        z_dean =  self.dean_parameter * y_dean ** (2/3)
        
        y_coords = y_dean + self.y_shoreline + self.distance_beach_face_slope
        z_coords = self.z_shoreline - self.distance_beach_face_slope * self.beach_face_slope - z_dean 
        coordinates = [(y, z) for y, z in zip(y_coords, z_coords)]
        return coordinates
        

    def build_rocky_profile(self):
        y_min, y_max, z_min, z_max = self.get_bounds(self.domain_bounds)
        platform_z_ymax = self.platform_z_coastline - self.platform_slope * y_max
        
        coords = ((0, self.platform_z_coastline), (y_max, platform_z_ymax), (y_max, z_min), 
                (y_min, z_min), (y_min, self.cliff_height), (0, self.cliff_height))
        self.rocky_polygon = geometry.Polygon(coords)
                                                
    @staticmethod
    def get_bounds(bounds):
        y_min = bounds['z_min']
        y_max = bounds['y_max']
        z_min = bounds['z_min']
        z_max = bounds['z_max']
        return y_min, y_max, z_min, z_max        
            
    def plot(self):
        
        fig, ax = plt.subplots()
        y_min, y_max, z_min, z_max = self.get_bounds(self.view_bounds)
        
        if self.rocky_polygon.area >0 :
            x, y = self.rocky_polygon.exterior.xy
            ax.fill(x, y, color = self.platform_color) 
      
        if self.beach_polygon.area >0 :
            x, y = self.beach_polygon.exterior.xy
            ax.fill(x, y, color = self.sand_color) 
        
      
        ax.set_xlim(y_min, y_max)
        ax.set_ylim(z_min, z_max)


profile1 = BeachProfile("profile.yaml", y_cl = 0)
profile1.plot()

