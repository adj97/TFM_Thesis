# LWR solver

import numpy as np

# DEFINE single road GEOMETRY

length = 5  # km
vmax = 90  # km/h
def demand(rho): return (90*rho)*(rho <= 30) + 2700*(rho > 30)
def supply(rho): return 2700*(rho <= 30) + (15*(30-rho)+2700)*(rho > 30)


# SUPPY AND DEMAND

demand_upstream = 1200  # veh/hr
supply_downstream = 2000  # veh/hr

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

for i in range(0, n_x):
    if X[0, i] <= 0.5:
        rho[0, i] = 20
    else:
        rho[0, i] = 100

# time loop
for i in range(1, n_t-1):

    # for
    j = 0
    inflow = min(demand_upstream, supply(rho[i, j]))
    outflow = min(demand(rho[i, j]), supply_downstream)
    rho[i, j] = rho[i-1, j]+(dt/dx)*(inflow - outflow)
    inflow = outflow

    for j in range(1, n_x-2):
        outflow = min(demand(rho[i, j]), supply_downstream)
        rho[i, j] = rho[i - 1, j] + (dt / dx) * (inflow - outflow)
        inflow = outflow

    # for
    j = n_x-1
    outflow = min(demand(rho[i, j]), supply_downstream)
    rho[i, j] = rho[i - 1, j] + (dt / dx) * (inflow - outflow)

np.savetxt('density.txt',rho)
