# LWR solver

import numpy as np
import os
from tqdm import tqdm

# DEFINE road characteristics

geometry = {}

# first road
geometry[1] = {}
length = 5  # km
vmax = 90  # km/h
def demand(rho): return (90*rho)*(rho <= 30) + 2700*(rho > 30)
def supply(rho): return 2700*(rho <= 30) + (15*(30-rho)+2700)*(rho > 30)
geometry[1]['length'] = length
geometry[1]['vmax'] = vmax
geometry[1]['demand'] = demand
geometry[1]['supply'] = supply

# second road
geometry[2] = {}
length = 5  # km
vmax = 90  # km/h
def demand(rho): return (90*rho)*(rho <= 30) + 2700*(rho > 30)
def supply(rho): return 2700*(rho <= 30) + (15*(30-rho)+2700)*(rho > 30)
geometry[2]['length'] = length
geometry[2]['vmax'] = vmax
geometry[2]['demand'] = demand
geometry[2]['supply'] = supply

# number of road sections
nb_link = len(geometry)

# NETWORK SUPPLY AND DEMAND

def demand_upstream(t): return 1200  # veh/hr
def supply_downstream(t): return 2000  # veh/hr

# SPATIAL STEP AND FINAL TIME
dx = 0.2  # km
T = 0.5  # hr

# Maximum network speed
V_max = geometry[1]['vmax']-1e-3
for road in geometry:
    V_max = max(V_max, geometry[road]['vmax'])

# CFL safety factor
CFL = 1.5
dt = dx/(CFL*V_max)

# total length
total_length = 0
for road in geometry:
    total_length += geometry[road]['length']

# space grid
x = np.arange(dx/2, total_length+dx/2, dx)

# Define initial density profile
Rho_0 = np.zeros(len(x))
for i in range(0, len(x)):
    if x[i] <= 0.5 + 1e-14:
        Rho_0[i] = 20
    elif x[i] <= 5 + 1e-14:
        Rho_0[i] = 100
    elif x[i] <= 5.5 + 1e-14:
        Rho_0[i] = 20
    else:
        Rho_0[i] = 100

# Density structure
Density = {}
for road in geometry:
    start = 0 if road == 1 else int(geometry[road-1]['length']/dx)
    end = start+int(geometry[road]['length']/dx)
    Density[road] = Rho_0[start:end]


# Set initial density as the first row of Rho array
n_t = len(np.arange(0, T, dt))
n_x = int(total_length/dx)
Rho = np.zeros((n_t, n_x))
Rho[0, ] = Rho_0

# GODUNOV SOLVER
def godunov(geometry, Rho_0, f_demand_upstream, f_supply_downstream, dx, T):

    # road length
    length = geometry['length']

    # number of x points
    n_x = int(length / dx)

    rho = np.zeros(len(Rho_0))

    # spatial range index j : 'loop'

    # first cell
    j = 0
    supply = geometry['supply']
    demand = geometry['demand']
    temp_demand_upstream = f_demand_upstream if isinstance(f_demand_upstream, np.float64) else f_demand_upstream((i+1)*dt)
    inflow = min(temp_demand_upstream, supply(Rho_0[j]))
    outflow = min(demand(Rho_0[j]), supply(Rho_0[j+1]))
    rho[j] = Rho_0[j]+(dt/dx)*(inflow - outflow)
    inflow = outflow

    # internal cells
    for j in range(1, n_x-1):
        supply = geometry['supply']
        demand = geometry['demand']
        outflow = min(demand(Rho_0[j]), supply(Rho_0[j+1]))
        rho[j] = Rho_0[j] + (dt/dx)*(inflow - outflow)
        inflow = outflow

    # end cell
    j = n_x-1
    demand = geometry['demand']
    temp_supply_downstream = f_supply_downstream if isinstance(f_supply_downstream, np.float64) else f_supply_downstream((i+1)*dt)
    outflow = min(demand(Rho_0[j]), temp_supply_downstream)
    rho[j] = Rho_0[j] + (dt/dx)*(inflow - outflow)

    return rho

# Junction function
def junction(geometry,A,rho_0):
    f_demand = geometry[1]['demand']
    f_supply = geometry[2]['supply']

    f_D = demand(rho_0[0])
    f_S = supply(rho_0[1])

    flows = np.zeros(len(rho_0))

    flows[0] = min(f_D, f_S)
    flows[1] = A*flows[0]

    return flows

# Time loop
i = 0  # time iteration index
for t in tqdm(np.arange(dt, T, dt)):

    # Junction function
    rho_0 = [Rho[i, int(geometry[1]['length']/dx)-1],
             Rho[i, int(geometry[1]['length']/dx)]]

    A = 1 # traffic distribution matrix
    flows = junction(geometry, A, rho_0)
    road_outflow = flows[0]
    road_inflow = flows[1]

    # Network global supply/demand definitions
    def net_glob_demand (road, t):
        n_g_d = demand_upstream(t)*(road == 1) + road_inflow*(road == 2)
        return n_g_d
    def net_glob_supply (road, t):
        n_g_s = road_outflow*(road == 1) + supply_downstream(t)*(road == 2)
        return n_g_s


    for road in geometry:
        # Update 'initial' density as the previous time step solution
        start = 0 if road == 1 else int(geometry[road - 1]['length'] / dx)
        end = start + int(geometry[road]['length'] / dx)
        rho_0 = Rho[i, start:end]

        # Run Godunov for each road over a single time step
        Rho[i+1, start:end] = godunov(geometry[road],
                                      rho_0,
                                      net_glob_demand(road, t),
                                      net_glob_supply(road, t),
                                      dx,
                                      dt)

        # First iteration complete
        # print(Rho[0, start:end])
        # print(godunov(geometry[road], rho_0, demand_upstream, supply_downstream, dx, dt))
        # exit()


    i = i+1

os.remove('density.txt')
np.savetxt('density.txt', Rho)
