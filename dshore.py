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



class Shore:
    def __init__(self, cfg_file, **kwargs):
        
 
        with open(cfg_file, 'r') as file:
            # dictionary defined in configuration file
            config_dict = yaml.safe_load(file)        
        
        for key, value in config_dict.items():
                setattr(self, key, value)
     
 
        for key, value in kwargs.items():
            setattr(self, key, value)

        self.shoreline = SEM_Line(**self.opt_shoreline)
        self.baseline = SEM_Line(**self.opt_baseline)
        self.coastline = SEM_Line(**self.opt_coastline)
        self.beach_toe = SEM_Line(**self.opt_beach_toe)
        self.create_computational_grid()

    def create_computational_grid(self):
        self.baseline.create_grid(self.dx, self.transect_length)
        self.coastline.xc, self.coastline.yc = self.transect_intersection(self.coastline.original_line,
                                                                          self.baseline.transects)
        self.beach_toe.xc, self.beach_toe.yc = self.transect_intersection(self.beach_toe.original_line,
                                                                          self.baseline.transects)

        b_w = self.dist(self.beach_toe, self.coastline)
        self.beach_width = b_w.flatten()

        self.beach_profile = [beachpy.BeachProfile(x_beachface_toe=beach_width) for beach_width in self.beach_width]
        self.compute_shoreline_from_profiles()

    @staticmethod
    def transect_intersection(line, transects):
        points = [line.intersection(transect) for transect in transects.geometry]
        xc = np.array([(point.x) for point in points])
        yc = np.array([(point.y) for point in points])
        return xc, yc

    @staticmethod
    def dist(line1, line2):
        return np.hypot(line1.xc - line2.xc, line1.yc - line2.yc)

    def compute_shoreline_from_profiles(self):
        dist_shoreline_to_coastline = np.array([profile.p_shoreline.x for profile in self.beach_profile])
        dist_shoreline_to_coastline.resize(dist_shoreline_to_coastline.size, 1)  # put in a compatible format
        dist_coastline_to_baseline = self.dist(self.coastline, self.baseline)
        self.shoreline.xc, self.shoreline.yc = self.baseline.xycoords_from_distance(
            dist_shoreline_to_coastline + dist_coastline_to_baseline)

    def compute_beachface_toe_from_profiles(self):
        dist_beachface_toe_to_coastline = np.array([profile.x_beachface_toe for profile in self.beach_profile])
        dist_beachface_toe_to_coastline.resize(dist_beachface_toe_to_coastline.size, 1)  # put in a compatible format
        dist_coastline_to_baseline = self.dist(self.coastline, self.baseline)
        self.beach_toe.xc, self.beach_toe.yc = self.baseline.xycoords_from_distance(
            dist_beachface_toe_to_coastline + dist_coastline_to_baseline)

    def propagate_waves(self, offshore_wave):

        nearshore_angle = self.baseline.azimuth_from_normals_degress()
        shoreline_angle = self.beach_toe.azimuth_from_normals_degress()

        self.offshore_waves = []
        self.nearshore_waves = []
        self.breaking_waves = []

        for transect_boundary in range(nearshore_angle.size):
            self.offshore_waves.append(deepcopy(offshore_wave))
            self.offshore_waves[transect_boundary].set_dir_bottom(nearshore_angle[transect_boundary])

            self.nearshore_waves.append(deepcopy(self.offshore_waves[transect_boundary]))
            self.nearshore_waves[transect_boundary].shoal(self.nearshore_depth)

            self.breaking_waves.append(deepcopy(self.nearshore_waves[transect_boundary]))
            self.breaking_waves[transect_boundary].set_dir_bottom(shoreline_angle[transect_boundary])
            self.breaking_waves[transect_boundary].break_it()

    def Q_potential(self):
        Q = np.array([self.dt * 0.233 * self.K * wave.H ** (5 / 2) * math.sin(np.deg2rad(wave.alpha()) * 2) for wave in
                      self.breaking_waves])  # CEM Eq III-2-7b (Rosati et al., 2002)

        if self.left_boundary == 'closed':
            Q[0] = 0
        elif self.left_boundary == 'open':
            Q[0] = Q[1]
        else:
            exit('SEM_Grid error: left boundary not defined')

        if self.right_boundary == 'closed':
            Q[-1] = 0
        elif self.right_boundary == 'open':
            Q[-2] = Q[-1]
        else:
            exit('SEM_Grid error: right boundary not defined')

        return Q

    def profile_volume(self):
        return np.array([profile.volume() for profile in self.beach_profile])

    def beach_volume(self):
        dx = self.baseline.cell_lenght().flatten()
        vol = self.profile_volume()
        return np.sum(dx * vol)

    def Q_net(self):
        return self.q_net(self.Q_potential(), self.profile_volume())

    def nextstep(self, offshore_wave):
        self.propagate_waves(offshore_wave)
        Q = self.Q_net()
        dv_cell = np.diff(Q)
        dv_profile = dv_cell / self.baseline.cell_lenght().flatten()
        self.update_profiles(-dv_profile)

    def update_profiles(self, dv):
        for dvol, profile in zip(dv, self.beach_profile):
            profile.update_volume(dvol)
        self.compute_beachface_toe_from_profiles()

    @staticmethod
    def divergent_cells(q_pot):
        return (np.sign(q_pot[1:]) - np.sign(q_pot[:-1])) == 2

    @staticmethod
    def recalculate_q(q_old, v, vnew):
        q_left = np.sign(q_old[:-1]) == -1
        q_right = np.sign(q_old[1:]) == 1

        vnew[vnew > 0] = 0

        v_correction = vnew.copy()
        v_correction[v_correction > 0] = 0

        div_cells = SEM_Grid.divergent_cells(q_old)
        v_correction[div_cells] = vnew[div_cells] / 2

        q_correction = q_old.copy()
        q_correction_left = q_correction[:-1].copy()
        q_correction_right = q_correction[1:].copy()

        q_correction_left[q_left] -= v_correction[q_left]
        q_correction_left[~q_left] = 0

        q_correction_right[q_right] += v_correction[q_right]
        q_correction_right[~q_right] = 0

        q_corrected = np.concatenate(([q_old[0]], q_correction_right)) + np.concatenate(
            (q_correction_left, [q_old[-1]]))

        q_corrected[0] = q_old[0]
        q_corrected[-1] = q_old[-1]

        return q_corrected

    @staticmethod
    def q_net(q_pot, v):

        div_cells = SEM_Grid.divergent_cells(q_pot)

        # at divergent cells potential transport if assumed to be equal at both sides 
        # this constrain should be only applicable at very specific locations  
        divergent_q = abs(np.array([q_pot[1:][div_cells], q_pot[:-1][div_cells]])).min(0)
        q_pot[:-1][div_cells] = -divergent_q
        q_pot[1:][div_cells] = divergent_q

        for x in range(SEM_Grid.max_iter_qpot2qnet):
            dq = -np.diff(q_pot)
            vnew = v + dq

            if np.sum(vnew < 0):
                q_pot = SEM_Grid.recalculate_q(q_pot, v, vnew)

            else:
                break

        return q_pot

    def plot(self):

        self.compute_shoreline_from_profiles()

        self.baseline.plot()
        self.shoreline.plot()
        self.coastline.plot()
        self.beach_toe.plot()

        if self.beachface_plot:
            x = np.concatenate([self.beach_toe.xc.flatten(), self.shoreline.xc[::-1].flatten()])
            y = np.concatenate([self.beach_toe.yc.flatten(), self.shoreline.yc[::-1].flatten()])
            plt.fill(x, y, color='gold')  # color names from https://matplotlib.org/examples/color/named_colors.html

        if self.beachberm_plot:
            x = np.concatenate([self.coastline.xc.flatten(), self.shoreline.xc[::-1].flatten()])
            y = np.concatenate([self.coastline.yc.flatten(), self.shoreline.yc[::-1].flatten()])
            plt.fill(x, y, color='sandybrown')

        ax = plt.gca()

        if self.transect_plot:
            self.baseline.transects.plot(ax=ax, color=self.transect_color, linestyle='--', linewidth=0.5)

        if self.cell_annotation:
            for i, x, y in zip(range(self.baseline.xc.size), self.baseline.xc, self.baseline.yc):
                ax.annotate('%s' % i, xy=(x, y), xytext=(5, 0), textcoords='offset points', fontsize=8)

        plt.axis('equal')

    def plot_waves(self):
        df = pandas.DataFrame({'H': [wave.H for wave in self.breaking_waves],
                               'alpha': [wave.alpha() for wave in self.breaking_waves]})
        fig, ax = plt.subplots()
        plt.plot(df['H'], '--r', label=r' $H_b$')
        plt.legend(loc='upper left')
        plt.xlabel('transect number')

        ax.twinx()
        plt.plot(df['alpha'], '--b', label=r'$\alpha_b$')
        plt.legend(loc='upper right')

    def plot_Q(self):

        fig, ax = plt.subplots()
        plt.plot(self.Q_net(), '-k', label=r' $Q_net$')
        plt.legend(loc='upper left')
        plt.xlabel('transect number')

        ax.twinx()
        plt.plot(self.Q_potential(), '--r', label=r' $Q_potential$')
        plt.legend(loc='upper right')
     
        
        
        
        
        
        
        

n = 50
q = np.array([0, -0.41, -0.85, 0.01 , 0.25, 0])

#array = np.random.uniform(-1, 1, size=n-1)
#q = np.insert(array, [0, len(array)], [0, 0])

q =  np.array([ 0.        , 1.70, -0.4,  0.7, -0.9,        0.        ])

q = np.random.uniform(-1, 1, n-1)

q = np.append([0], q)
q = np.append(q, [0])

#vol = np.array([0, 0, 10, 0, 1.])

vol = np.ones((n))

q_out_left_pot = q[:-1]
q_out_right_pot = q[1:]

div_q = + q_out_right_pot - q_out_left_pot 
      

ncells = n

time = 2000
for t in range(time):
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
    
    vol = vol_old - div_q 
    print(vol.sum(), '\n\n')    




