# -*- coding: utf-8 -*-
"""
Created on Tue Mar 12 14:59:16 2024

@author: Rui
"""
from shapely import geometry
from shapely.geometry import Point
import matplotlib.pyplot as plt
import numpy as np
import yaml
from copy import deepcopy


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


    def build_beach_profile_from_area(self, area):
        c = Point(self.x, 0, self.c_z)
        mb =self.bf_m #beach face slope
        mp = self.p_m # platform slope
        pz = self.c_p_z 
        h  = c.z - pz 
        l =  h / (mb - mp)
        critical_area = h * l / 2
        if area <= critical_area:   # beach face only
            t_y = np.sqrt(area * 2 / ((mb - mp)))
            t_z = pz - mp * t_y 
            c_z = (mb - mp) * t_y + pz
            s_y = 0
            s_z = c_z
            
        else:
            if mp == 0:
               s_y = -(-2*area*mb + c.z**2 - 2*c.z*pz + pz**2) / (2*c.z* mb - 2*mb*pz)
               
               s_z = c.z
               t_y = (s_z - pz) / mb + s_y
               t_z = pz
               print(s_y, s_z, t_y, t_z)
               
            else:
                b = ((2*np.sqrt(2*area * mp + c.z**2 - 2*c.z*pz + pz**2) /np.sqrt(mb**2-mb*mp)) - 2*c.z/mb + 2*pz/(mb-mp))/(2*(1/(mb-mp) - 1/mb))
                s_y = (b - c.z) / mb
                s_z = c.z
                t_y = (pz - b) / (mp - mb)
                t_z =  -mp * t_y + pz

        
        c = Point(self.x, 0, c.z)
        s = Point(self.x, s_y, s_z)
        t = Point(self.x, t_y, t_z)
            
            
        self.c = c
        self.s = s
        self.t = t
        coords = ((c.y, c.z), (s.y, s.z), (t.y, t.z), (c.y, self.c_p_z), (c.y, c.z))
        self.beach_polygon = geometry.Polygon(coords)
             
            
    
    def build_beach_profile(self):
        c = Point(self.x, 0, self.c_z)
        s = Point(self.x, c.y + self.b_w, c.z)
       
        b = s.z + self.bf_m * s.y
        t_y = (self.c_p_z - b) / (self.p_m - self.bf_m)
        t_z = -self.p_m * t_y + self.c_p_z
        t = Point(self.x, t_y, t_z)
        coords = ((c.y, c.z), (s.y, s.z), (t.y, t.z), (c.y, self.c_p_z), (c.y, c.z))
        self.beach_polygon = geometry.Polygon(coords)
        self.c = c
        self.s = s
        self.t = t
        
    def build_rocky_profile(self):
        y_min, y_max, z_min, z_max = self.get_bounds(self.domain_bounds)
        c_p = Point(self.x, 0, self.c_p_z)
        ymax_p = Point(self.x, y_max, c_p.z - self.p_m * y_max)
        coords = ((c_p.y, c_p.z), (ymax_p.y, ymax_p.z), (y_max, z_min), 
                (y_min, z_min), (y_min, self.cliff_elevation), (0, self.cliff_elevation))
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


p = BeachProfile("beach_profile.yaml")
print(p.beach_polygon.area)
p.plot()
p2 = deepcopy(p)
p2.build_beach_profile_from_area(1000)
p2.plot()
