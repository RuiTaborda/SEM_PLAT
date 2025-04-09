# -*- coding: utf-8 -*-
"""
Created on Fri Apr 26 18:14:50 2024

@author: Rui
"""

import shoreline_model as SEM
import wavetimeseries as WTS
import numpy as np
import imageio


create_movie = True

waves = WTS.WaveTimeSeries(filename = 'PT_1979_2023_sea_and_swell.nc', datafile_type = 'era5')

grid = SEM.Grid('sem_grid_teo.yaml')

grid.plot()

w=waves.wave_data.iloc[0]

w.Dir = 20
w.Hs = 2
w.Tp = 10


yc = grid.shoreline.yc
v1 = grid.volume

ax = grid.plot()

frames = []
for i in range(3000):
    grid.next_step(w)
    if create_movie and i % 20 == 0:
            ax.cla()
            grid.plot(ax = ax)
            fig = ax.figure

            fig.canvas.draw()
            image = np.frombuffer(fig.canvas.tostring_rgb(), dtype='uint8')
            image = image.reshape(fig.canvas.get_width_height()[::-1] + (3,))
    
            # Append the color image to the frames
            frames.append(image)
    
    print(i)

imageio.mimsave('teo2.gif', frames, fps=2)

for i in range(3000):
    grid.next_step(w)
    
dv = grid.volume - v1

grid.plot()