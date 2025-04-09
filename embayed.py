# -*- coding: utf-8 -*-
"""
Created on Fri Apr 26 18:14:50 2024

@author: Rui
"""

import shoreline_model as SEM
import wavetimeseries as WTS
import matplotlib.pyplot as plt
import numpy as np
import imageio

import wavepy as wv
#waves = WTS.WaveTimeSeries(filename = 'PT_1979_2023_sea_and_swell.nc', datafile_type = 'era5')

grid = SEM.Grid('sem_grid_embayed.yaml')

grid.plot()
plt.savefig('figinitial.svg')

create_movie = True

#w=waves.wave_data.iloc[0]
w = wv.Wave(Dir = 330, Hs = 1, Tp = 8)



yc = grid.shoreline.yc
ax = grid.plot()

frames = []
for i in range(3000):
    grid.next_step(w)
    if create_movie and i % 50 == 0:
            ax.cla()
            grid.plot(ax = ax)
            fig = ax.figure

            fig.canvas.draw()
            image = np.frombuffer(fig.canvas.tostring_rgb(), dtype='uint8')
            image = image.reshape(fig.canvas.get_width_height()[::-1] + (3,))
    
            # Append the color image to the frames
            frames.append(image)
    
    print(i)

imageio.mimsave('coast.gif', frames, fps=2)

plt.savefig('fig.svg')


w.Dir = 300
w.Hs = 1
for i in range(150):
    grid.next_step(w)
grid.plot()
plt.savefig('fig270.svg')

w.Dir = 290
w.Hs = 1

for i in range(150):
    grid.next_step(w)

grid.plot()
plt.savefig('fig300.svg')

w.Dir = 280
w.Hs = 1

for i in range(150):
    grid.next_step(w)

grid.plot()
plt.savefig('fig300.svg')



w.Dir = 320
w.Hs = 1

yc = grid.shoreline.yc

for i in range(200):
    grid.next_step(w)
    print(i)

grid.plot()
plt.savefig('fig320.svg')

