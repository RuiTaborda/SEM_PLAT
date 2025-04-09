# -*- coding: utf-8 -*-
"""
Created on Fri Apr 26 18:14:50 2024

@author: Rui
"""

import shoreline_model as SEM
import wavetimeseries as WTS

#waves = WTS.WaveTimeSeries(filename = 'PT_1979_2023_sea_and_swell.nc', datafile_type = 'era5')

grid = SEM.Grid('sem_grid.yaml')



grid.plot()




