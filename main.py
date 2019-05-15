# LWR solver
# Upwind discretisation

import numpy as np
import matplotlib.pyplot as plt

nx = 10  # number of nodes
x0 = 0   # 'left' x value
xm = 10  # 'right' x value
dx = (xm-x0)/nx  # x increment
tf = 10  # final time
nt = 10  # number of time steps
dt = tf / nt  # time step
uf = 10  # free velocity
rho_m = 3  # maximum (jam) density

# create empty rho arrays
rho_c = np.zeros(nx+1)  # 'current' rho array
rho_n = rho_c              # 'next' rho array

# Final arrays
density = flow = velocity = np.zeros((nx+1, nt-1))
print(density.shape)

x = np.linspace(x0, xm, nx+1)  # create space grid

for i in range(0, nx+1):  # initialise rho distribution
    rho_c[i] = 2

u = uf*(1-rho_c/rho_m)  # Greenshields velocity
f = rho_c*u             # Flow-density-speed relation

# Store initial values in large arrays
density[:, 0] = rho_c
velocity[:, 0] = u
flow[:, 0] = f

for t in np.linspace(1, tf, num=nt):  # time loop
    time_index = int(round(t/dt))

    for i in range(0, nx+1):  # space loop
        rho_n[i] = rho_c[i] - (dt/dx)*(f[i]-f[i-1])  # Upwind update scheme

    u = uf * (1 - rho_n / rho_m)  # Greenshields velocity
    f = rho_n * u  # Flow-density-speed relation

    rho_c = rho_n  # save the new solution as 'current' for the next time step

    density[:, time_index] = rho_c
    velocity[:, time_index] = u
    flow[:, time_index] = f
