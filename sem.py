# -*- coding: utf-8 -*-
"""
Created on Tue Mar 26 13:17:13 2024

@author: Rui
"""
import numpy as np



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




