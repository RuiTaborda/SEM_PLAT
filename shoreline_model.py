# -*- coding: utf-8 -*-
"""
Created on Tue Mar 26 13:17:13 2024

@author: Rui
"""
import numpy as np

import geopandas as gpd
from shapely import geometry
from shapely.geometry import Point
import matplotlib.pyplot as plt
import yaml
import pandas as pd
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


    def update_beach_profile_from_area(self, area):
        if area < 1e-2:
            area = 0
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
               
               
            else:
                b = ((2*np.sqrt(2*area * mp + c.z**2 - 2*c.z*pz + pz**2) /np.sqrt(mb**2-mb*mp)) - 2*c.z/mb + 2*pz/(mb-mp))/(2*(1/(mb-mp) - 1/mb))
                s_y = (b - c.z) / mb
                s_z = c.z
                t_y = (pz - b) / (mp - mb)
                t_z =  -mp * t_y + pz

        print(s_y, s_z, t_y, t_z)  
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



class Line:
    def __init__(self, cfg_file, **kwargs):
        
        with open(cfg_file, 'r') as file:
            # dictionary defined in configuration file
            config_dict = yaml.safe_load(file)        
    
        for key, value in config_dict.items():
            setattr(self, key, value)
        
        for key, value in kwargs.items():
            setattr(self, key, value)

        if  not hasattr(self, 'filename'):
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

    def create_grid(self, dx, transect_length):  # adpated from https://stackoverflow.com/questions/34906124/interpolating-every-x-distance-along-multiline-in-shapely/35025274#35025274
        geom = self.original_line.geometry
        self.num_cells = int(round(geom.length / dx).iloc[0])
        points = [geom.interpolate(float(n) / self.num_cells, normalized = True) for n in range(self.num_cells + 1)]

        self.x = np.squeeze(np.array([p.x for p in points]))
        self.y = np.squeeze(np.array([p.y for p in points]))

        self.compute_mid_cells()
        self.build_transects(transect_length)

    def xycoords_from_distance(self, dist):
        xc = self.xc + dist * np.cos(self.trans_dir)
        yc = self.yc + dist * np.sin(self.trans_dir)
        return xc, yc

    def plot_shp(self, ax = None):
        xy = self.original_line.geometry[0].coords.xy

        if ax is None:
            fig, ax = plt.subplots()
            
        
        ax.plot(xy[0], xy[1], color=self.shp_color, markersize=self.shp_markersize,
                 marker=self.shp_marker, markeredgecolor=self.shp_markeredgecolor, linestyle=self.shp_linestyle,
                 linewidth=self.shp_linewidth)
        ax.axis('equal')
        return ax

    def plot_x(self, ax = None):
        if ax is None:
            fig, ax = plt.subplots()
        
        ax.plot(self.x, self.y, color=self.x_color, markersize=self.x_markersize,
                 marker=self.x_marker, markeredgecolor=self.x_markeredgecolor, linestyle=self.x_linestyle,
                 linewidth=self.x_linewidth)
        ax.axis('equal')
        return ax


    def plot_xc(self, ax = None):
        if ax is None:
            fig, ax = plt.subplots()

        ax.plot(self.xc, self.yc, color=self.xc_color, markersize=self.xc_markersize,
                 marker=self.xc_marker, markeredgecolor=self.xc_markeredgecolor, linestyle=self.xc_linestyle,
                 linewidth=self.xc_linewidth)
        ax.axis('equal')
        return ax


    def plot(self, ax = None):
        
        if self.shp_plot:
            ax = self.plot_shp(ax = ax)
        if self.xc_plot:
            ax = self.plot_xc(ax = ax)
        if self.x_plot:
            ax = self.plot_x(ax = ax)
            
        return ax
        

class Grid:
    def __init__(self, cfg_file, **kwargs):
        
        with open(cfg_file, 'r') as file:
            # dictionary defined in configuration file
            config_dict = yaml.safe_load(file)        
    
        for key, value in config_dict.items():
            setattr(self, key, value)
        
        for key, value in kwargs.items():
            setattr(self, key, value)
            
        self.baseline = Line(self.sem_line_style, **self.opt_baseline)
        self.shoreline = Line(self.sem_line_style, **self.opt_shoreline)
        self.coastline = Line(self.sem_line_style, **self.opt_coastline)
        self.beach_toe = Line(self.sem_line_style, **self.opt_beach_toe)
        
        self.nearshore_contour = Line(self.sem_line_style, **self.opt_nearshore_contour)
        self.create_computational_grid()
        
        if self.validation_filename:
            self.validation_lines = gpd.read_file(self.validation_filename)
            xc, yc = self.transect_intersection(self.validation_lines, self.baseline.transects, multiple_lines = True)
            dist = np.hypot(self.baseline.xc.T - xc, self.baseline.yc.T - yc)
            self.validation_dist = dist - dist.iloc[0]
        
    def create_computational_grid(self):
        self.baseline.create_grid(self.dx, self.transect_length)
        
        # build coastline - only defined at cell centers. 
        self.coastline.xc, self.coastline.yc = self.transect_intersection(self.coastline.original_line, self.baseline.transects)
        self.shoreline.xc, self.shoreline.yc = self.transect_intersection(self.shoreline.original_line, self.baseline.transects)
        
        #self.nearshore_contour.xc, self.nearshore_contour.yc = self.transect_intersection(self.nearshore_contour.original_line, self.baseline.transects)
        self.b_w = self.dist(self.shoreline, self.coastline)
        
        n = self.baseline.num_cells
        self.profiles = [BeachProfile(self.beach_profile_yaml, b_w = self.b_w[i]) for i in range(n)]
        self.volume = [profile.beach_polygon.area for profile in self.profiles]
 
        
        self.beach_toe.dist = [profile.t.y for profile in self.profiles]
        self.beach_toe.xc, self.beach_toe.yc = self.baseline.xycoords_from_distance(self.beach_toe.dist + self.dist(self.baseline, self.coastline))
        
 
    def read_validation_lines(self, filename):
        self.validation_lines = gpd.read_file(filename)
        xv, yv = self.transect_intersection(self.validation_lines, self.baseline.transects)
        

    def Q_pot(self, wave):
        local_angle = self.beach_toe.azimuth_from_normals_degress() - wave.Dir
        
        if self.Q_method == 'AM':
            
            k2 = np.sqrt(self.g*self.gama/(2* np.pi))**(1/5) * self.k1
            return k2*wave.Hs**(12./5) * wave.Tp**(1./5)* np.cos(np.radians(local_angle))**(6/5)*np.sin(np.radians(local_angle))
        
        else:
            b = self.K * self.rho * np.sqrt(self.g / self.K) / (16 * (self.rho_sed - self.rho) * (1 - self.porosity))
            return  b * wave.Hs**(5/2) * np.sin(np.radians(2 * local_angle))

    def next_step(self, wave):
        
        q = self.Q_pot(wave)
        
        if self.boundary['left'] == 'closed':
            q[0] =0
        if self.boundary['right'] == 'closed':
             q[-1] = 0
        
        if self.boundary['left'] == 'gated':
            q[0] = min(q[0], 0)
   
            
        if self.boundary['right'] == 'gated':
            q[-1] = max(q[-1], 0)
        
        
        q_out_left_pot = q[:-1]
        q_out_right_pot = q[1:]
    
        div_q = + q_out_right_pot - q_out_left_pot 
    
        vol = self.volume    
        
        print(np.array(vol).sum(), '\n\n')   
    
        ncells = self.baseline.num_cells
        
        t = 0
    
        while t < self.dt_next_step:
            vol_old = vol.copy()
            
            
            q_out_left_pot = q[:-1]
            q_out_right_pot = q[1:]

            q_net = q.copy()    
    
            div_q = -q_out_left_pot + q_out_right_pot
                    
            for i in range( ncells):
              if (vol_old[i] < -q_out_left_pot[i]): 
                    q_net[i] = -vol_old[i]
                    
              if (vol_old[i] < q_out_right_pot[i]): 
                    q_net[i+1] = vol_old[i]
    
                  #check for divergent cells   
              if (q_out_left_pot[i] < 0) & (q_out_right_pot[i] > 0):
                    if vol_old[i] < div_q[i]:
                        q_net[i+1] = vol_old[i] * q_out_right_pot[i]/div_q[i]
                        q_net[i] = vol_old[i] * q_out_left_pot[i]/div_q[i]
    
            q_out_left_pot = q_net[:-1] 
            q_out_right_pot = q_net[1:]
    
            div_q = -q_out_left_pot + q_out_right_pot
            if self.dt == 0: #adaptative time step
                q_max = max(np.abs(q_out_left_pot).max(), np.abs(q_out_right_pot).max())
                dt = self.dx**2 / (4*np.abs(q_max))
            else:
                dt = self.dt
            
            t += dt
    
            dq = div_q*dt/self.dx
            
            vol = vol_old - dq 
            
        
        #print('t = ', t , dt)
            
        self.volume = vol
        
        [self.profiles[i].update_beach_profile_from_area(vol[i]) for i in range(self.baseline.num_cells)]

        self.shoreline.dist = [profile.s.y for profile in self.profiles]
        self.shoreline.xc, self.shoreline.yc = self.baseline.xycoords_from_distance(self.shoreline.dist + self.dist(self.baseline, self.coastline))
 
        self.beach_toe.dist = [profile.t.y for profile in self.profiles]
        self.beach_toe.xc, self.beach_toe.yc = self.baseline.xycoords_from_distance(self.beach_toe.dist + self.dist(self.baseline, self.coastline))
 
        
            
        
         
    
    



    
    

    @staticmethod
    def transect_intersection(line, transects, multiple_lines = False, date_field = 'DataFinal'):
        points = [line.intersection(transect) for transect in transects.geometry]
        
        xc = np.array([])
        yc = np.array([])
        for point in points:
            if point[0].is_empty:
                xc = np.append(xc, np.nan)
                yc = np.append(yc, np.nan)
            else:
                xc = np.append(xc, np.array(point.x))
                yc = np.append(yc, np.array(point.y))
        
        if multiple_lines:
            dates = line[date_field]
        
            n_transects = xc.shape[0]
            xc = pd.DataFrame(xc.T, index=dates, columns=[f'T_{i}' for i in range(1, n_transects + 1)])
            yc = pd.DataFrame(yc.T, index=dates, columns=[f'T_{i}' for i in range(1, n_transects + 1)])
            return xc, yc

        else:
        
            return xc, yc

    @staticmethod
    def dist(line1, line2):
        return np.hypot(line1.xc - line2.xc, line1.yc - line2.yc)


 
        


    def plot(self, ax = None):
       
        if ax == None:
            ax = self.baseline.plot()
        else:
            self.baseline.plot(ax = ax)
            
        if self.transect_plot:
          self.baseline.transects.plot(ax=ax, color=self.transect_color, linestyle='--', linewidth=0.5)
 
        self.coastline.plot(ax = ax)
        self.shoreline.plot(ax = ax)
        self.beach_toe.plot(ax = ax)
    
        #self.nearshore_contour.plot(ax = ax)
        
        
        if self.beachface_plot:
            x = np.concatenate([self.beach_toe.xc.flatten(), self.shoreline.xc[::-1].flatten()])
            y = np.concatenate([self.beach_toe.yc.flatten(), self.shoreline.yc[::-1].flatten()])
            plt.fill(x, y, color='gold')  # color names from https://matplotlib.org/examples/color/named_colors.html
    
        if self.beachberm_plot:
            x = np.concatenate([self.coastline.xc.flatten(), self.shoreline.xc[::-1].flatten()])
            y = np.concatenate([self.coastline.yc.flatten(), self.shoreline.yc[::-1].flatten()])
            plt.fill(x, y, color='sandybrown')
    
    
    
        if self.cell_annotation:
            for i, x, y in zip(range(self.baseline.xc.size), self.baseline.xc, self.baseline.yc):
                ax.annotate('%s' % i, xy=(x, y), xytext=(5, 0), textcoords='offset points', fontsize=8)

        plt.show()
        return ax
    
      
    
        
