# -*- coding: utf-8 -*-
"""
Created on Tue Jun 11 09:34:13 2024

@author: Rui
"""

import math
import sympy as sp
import numpy as np

# Given values (example)
mb = 0.1
mp = 0.01
pz = -1
cz = 4.0
A = 195



b1 = ((2*np.sqrt(2*A * mp + cz**2 - 2*cz*pz + pz**2) /np.sqrt(mb**2-mb*mp)) - 2*cz/mb + 2*pz/(mb-mp))/(2*(1/(mb-mp) - 1/mb))
b2 = (-2*np.sqrt(2*A * mp + cz**2 - 2*cz*pz + pz**2) /np.sqrt(mb**2-mb*mp) - 2*cz/mb + 2*pz/(mb-mp))/(2*(1/(mb-mp) - 1/mb))
                
sy = (b1 - cz) / mb
ty = (b1+ pz) / (mb - mp)  + sy
