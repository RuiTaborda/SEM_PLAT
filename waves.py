# -*- coding: utf-8 -*-
"""
Created on Tue Apr 30 21:46:10 2024

@author: rtabo
"""

import wavetimeseries as wts
import matplotlib.cm as cm
import matplotlib.pyplot as plt

era5 = wts.WaveTimeSeries(filename = 'C:/Users/rtabo/OneDrive - Universidade de Lisboa/mytools/wave/WaveTools/PT_1979_2021_sea_and_swell.nc', datafile_type = 'era5', lat = 39, long = -9.5, label_style = 'default')
