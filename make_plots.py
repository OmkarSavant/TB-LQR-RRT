# -*- coding: utf-8 -*-
"""
Created on Tue May 12 15:24:15 2020

@author: samwe
"""
import matplotlib.pyplot as plt
import time
from mpl_toolkits import mplot3d
import numpy as np


# %% plot tree

fig = plt.figure()
ax = plt.axes(projection='3d')

#theta
xdata = [node.state_time[0] for node in rrt.Tree]

#theta_dot
ydata = [node.state_time[1] for node in rrt.Tree]

#time
zdata = [node.state_time[2] for node in rrt.Tree]

ax.scatter3D(xdata,ydata,zdata, color = 'blue',alpha = 0.3)
ax.set_xlabel('theta (rad)')
ax.set_ylabel('theta_dot (rad/s)')
ax.set_zlabel('time (s)')

#plot nodes on path
x_path = [node.state_time[0] for node in path]
y_path = [node.state_time[1] for node in path]
z_path = [node.state_time[2] for node in path]

ax.plot(x_path,y_path,z_path, linewidth = 5.0, color = 'red', alpha=1.0)

ax.scatter3D(rrt.goal[0], rrt.goal[1], rrt.goal[2], s = 150.0, color = 'green', alpha = 1.0)
ax.scatter3D(rrt.start[0], rrt.start[1], rrt.start[2], s = 100.0, color = 'black', alpha = 1.0)

# %% plot

control = [node.u for node in path]
plt.plot(np.linspace(0, len(control)*rrt.dt, len(control)), control, linewidth=3.0)
plt.xlabel('Timestep')
plt.ylabel('Control Input')
plt.grid()
plt.legend(["Total Cost: " + str(cost)], loc="upper")