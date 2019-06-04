#
#   Andrew Dixon
#
#   Created   03/05/2019
#   Modified  03/06/2019
#
#   Thesis project for Cranfield University
#   MSc in Computational Fluid Dynamics
#
#   Description :
#

import os                    # standard
import numpy as np           # numerical programming
from tqdm import tqdm        # progressbar in time loop

# Define Network and Road Characteristics
geometry = {}

# Road 1
geometry[1] = {}
length = 5  # km
vmax = 90  # km/h
def demand(rho): return (90*rho)*(rho <= 30) + 2700*(rho > 30)
def supply(rho): return 2700*(rho <= 30) + (15*(30-rho)+2700)*(rho > 30)
geometry[1]['length'] = length
geometry[1]['vmax'] = vmax
geometry[1]['demand'] = demand
geometry[1]['supply'] = supply

# Road 2
geometry[2] = {}
length = 5  # km
vmax = 90  # km/h
def demand(rho): return (90*rho)*(rho <= 30) + 2700*(rho > 30)
def supply(rho): return 2700*(rho <= 30) + (15*(30-rho)+2700)*(rho > 30)
geometry[2]['length'] = length
geometry[2]['vmax'] = vmax
geometry[2]['demand'] = demand
geometry[2]['supply'] = supply

# Road 3
geometry[3] = {}
length = 5  # km
vmax = 90  # km/h
def demand(rho): return (90*rho)*(rho <= 30) + 2700*(rho > 30)
def supply(rho): return 2700*(rho <= 30) + (15*(30-rho)+2700)*(rho > 30)
geometry[3]['length'] = length
geometry[3]['vmax'] = vmax
geometry[3]['demand'] = demand
geometry[3]['supply'] = supply

# Number of road sections
nb_link = len(geometry)

# Network inlet demand and outlet supply
# What if the network has more than one inlet or outlet
def demand_upstream(t): return 1200  # veh/hr
def supply_downstream(t): return 2000  # veh/hr

# Spatial resolution
dx = 0.2  # km

# Final time
T = 0.5  # hr

# Global maximum network speed
V_max = geometry[1]['vmax']-1e-3
for road in geometry:
    V_max = max(V_max, geometry[road]['vmax'])

# Courant-Friedrichs-Lewy (CFL) safety factor
CFL = 1.5

# Time resolution from CFL constraint
dt = dx/(CFL*V_max)

# Total network length
total_length = 0
for road in geometry:
    total_length += geometry[road]['length']

# Spatial grid
x = np.arange(dx/2, total_length+dx/2, dx)

# Define initial global density profile
Rho_0 = np.zeros(len(x))
for i in range(0, len(x)):
    if x[i] <= 0.5 + 1e-14:
        Rho_0[i] = 100
    #elif x[i] <= 5 + 1e-14:
    #    Rho_0[i] = 100
    #elif x[i] <= 5.5 + 1e-14:
    #    Rho_0[i] = 20
    #elif x[i] <= 10 + 1e-14:
    #    Rho_0[i] = 100
    #elif x[i] <= 10.5 + 1e-14:
    #    Rho_0[i] = 20
    else:
        Rho_0[i] = 20

# Density structure
# In the future this will be the output density data structure
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

# Godunov solver
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
    f_demand_1 = geometry[1]['demand']
    f_supply_2 = geometry[2]['supply']
    f_supply_3 = geometry[3]['supply']

    f_D1 = f_demand_1(rho_0[0])
    f_S2 = f_supply_2(rho_0[1])
    f_S3 = f_supply_3(rho_0[2])

    flows = np.zeros(len(rho_0))

    flows[0] = min(f_D1, min(f_S2/A[0], f_S3/A[1]))
    flows[1] = A[0]*flows[0]
    flows[2] = A[1]*flows[0]

    return flows

# Time loop
# for t in tqdm(np.arange(dt, T, dt)): progress bar
for t in np.arange(dt, T, dt):

    # Loop for each junction
    # Each junction will have its own traffic distribution matrix

    # Iteration index from time t
    i = int(round(t/dt)-1)

    # Get this one junction boundary densities
    # The Density array will be useful here for the general junction
    rho_0 = [Rho[i, int(geometry[1]['length']/dx)-1],
             Rho[i, int(geometry[1]['length']/dx)],
             Rho[i, int((geometry[1]['length']+geometry[2]['length'])/dx)]]

    # Define traffic distribution matrix
    A = np.array([0.8, 0.2])

    # Call junction to get input/output flows
    flows = junction(geometry, A, rho_0)
    road_outflow = flows[0]
    road_inflow_2 = flows[1]
    road_inflow_3 = flows[2]

    # Network global supply/demand definitions
    def net_glob_demand(road, t):
        n_g_d = demand_upstream(t)*(road == 1) + road_inflow_2*(road == 2) + road_inflow_3*(road == 3)
        return n_g_d
    def net_glob_supply(road, t):
        n_g_s = road_outflow*(road == 1) + supply_downstream(t)*(road == 2) + supply_downstream(t)*(road == 3)
        return n_g_s

    # Each road seperately
    for road in geometry:

        # Update the new 'initial' density as the previous-time-step solution
        if road == 1:
            start = 0
        else:
            start = 0
            for road_index in range(1, road):
                start += int(geometry[road]['length'] / dx)

        end = start + int(geometry[road]['length'] / dx)
        rho_0 = Rho[i, start:end]

        # Run Godunov for each road over a single time step
        Rho[i+1, start:end] = godunov(geometry[road],
                                      rho_0,
                                      net_glob_demand(road, t),
                                      net_glob_supply(road, t),
                                      dx,
                                      dt)


# Update the saved density profile in pwd
os.remove('density.txt')
np.savetxt('density.txt', Rho)
