# -*- coding: utf-8 -*-
"""
Created on Fri Apr 26 18:14:50 2024

@author: Rui
"""

import shoreline_model as SEM
import wavetimeseries as WTS
import numpy as np

waves = WTS.WaveTimeSeries(filename = 'PT_1979_2023_sea_and_swell.nc', datafile_type = 'era5')

grid = SEM.Grid('sem_grid_aguda.yaml')

grid.plot()

w=waves.wave_data.iloc[0]

w.Dir = 300
w.Hs = 0.5

yc = grid.shoreline.yc
v1 = grid.volume

for i in range(100):
    grid.next_step(w)
dv = grid.volume - v1

grid.plot()
