# LWR solver
# Upwind discretisation

import numpy as np

nx = 1000+1  # number of nodes
x0 = -200   # 'left' x value
xm = 200  # 'right' x value
dx = (xm-x0)/nx  # x increment
tf = 5  # final time
nt = 500  # number of time steps
dt = tf / nt  # time step
uf = 25  # free velocity
rho_m = 0.04  # maximum (jam) density

# create empty rho arrays
rho_c = np.zeros(nx+1)  # 'current' rho array
rho_n = rho_c              # 'next' rho array

# dimension arrays
x = np.linspace(x0, xm, nx+1)
time = np.linspace(dt, tf, num=nt)

# Final arrays
density = flow = velocity = np.zeros((nx+1, nt+1))

for i in range(0, nx+1):  # initialise rho distribution
    rho_c[i] = (1+abs(x[i])/x[i])/100+0.01

u = uf*(1-rho_c/rho_m)  # Greenshields velocity
f = rho_c*u             # Flow-density-speed relation

# Store initial values in large arrays
density[:, 0] = rho_c
velocity[:, 0] = u
flow[:, 0] = f

for t in time:  # time loop
    time_index = int(round(t/dt))

    for i in range(0, nx+1):  # space loop
        rho_n[i] = rho_c[i] - (dt/dx)*(f[i]-f[i-1])  # Upwind update scheme

    u = uf * (1 - rho_n / rho_m)  # Greenshields velocity
    f = rho_n * u  # Flow-density-speed relation

    rho_c = rho_n  # save the new solution as 'current' for the next time step

    density[:, time_index] = rho_c
    velocity[:, time_index] = u
    flow[:, time_index] = f

# Plot results

# import matplotlib.pyplot as plt
# from matplotlib.animation import FuncAnimation
# plt.style.use('seaborn-pastel')
#
# fig = plt.figure()
# ax = plt.axes(xlim=(-200, 200), ylim=(0, 0.05))
# line, = ax.plot([], [], lw=2)
#
# def init():
#     line.set_data([], [])
#     return line,
# def animate(i):
#     y = density[:,i]
#     line.set_data(x,y)
#     return line,
#
# anim = FuncAnimation(fig, animate, init_func=init, frames=500, interval=20, blit=True)
# anim.save('density.gif', writer='imagemagick')

# for i in range(0,nt+1):
#     xvals = x
#     yvals = density[:, i]
#     plt.plot(xvals,yvals)
#     plt.show()