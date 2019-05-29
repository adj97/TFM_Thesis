# LWR solver
# Upwind discretisation

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation

# DEFINE GEOMETRY

length = 5  # km
vmax = 90  # km/h
def demand(rho): return (90*rho)*(rho <= 30) + 2700*(rho > 30)
def supply(rho): return 2700*(rho <= 30) + (15*(30-rho)+2700)*(rho > 30)


# INITIAL DENSITY PROFILE

initial_density = 50  # veh/km

# SUPPY AND DEMAND

demand_upstream = 1000  # veh/hr
supply_downstream = 1500  # veh/hr

# SPATIAL STEP AND FINAL TIME

dx = 0.2  # km
T = 0.5  # hr

# GODUNOV SOLVER

# total length
total_length = length

# cell index holder
X = np.zeros((2, int(total_length/dx)))
X[0, ] = np.arange(dx/2, total_length+dx/2, dx)
X[1, ] = np.ones(int(total_length/dx), dtype=int)

# CFL safety factor
CFL = 1.5
dt = dx/(CFL*vmax)

# initialise
n_t = len(np.arange(0, T, dt))
n_x = len(X[0, ])
rho = np.zeros((n_t, n_x))
rho[0, ] = 1250*np.ones((1, n_x))

# time loop
for i in range(1, n_t-1):

    # for
    j = 0
    inflow = supply_downstream  # min(demand_upstream, supply)
    outflow = demand_upstream  # min(demand, supply_downstream)
    rho[i, j] = rho[i-1, j]+(dt/dx)*(inflow - outflow)
    inflow = outflow

    for j in range(1, n_x-2):
        outflow = demand_upstream
        rho[i, j] = rho[i - 1, j] + (dt / dx) * (inflow - outflow)
        inflow = outflow

    # for
    j = n_x-1
    outflow = demand_upstream
    rho[i, j] = rho[i - 1, j] + (dt / dx) * (inflow - outflow)

# Plot results

plt.style.use('seaborn-pastel')

fig = plt.figure()
ax = plt.axes(xlim=(0, length), ylim=(0, 2000))
line, = ax.plot([], [], lw=2)

def init():
    line.set_data([], [])
    return line,
def animate(i):
    y = rho[i, ]
    line.set_data(X[0, ],y)
    return line,

anim = FuncAnimation(fig, animate, init_func=init, frames=n_t, interval=20, blit=True)
anim.save('density.gif', writer='imagemagick')
