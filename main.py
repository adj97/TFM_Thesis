#
#   Andrew Dixon
#
#   Created   03/05/2019
#   Modified  05/06/2019
#
#   Thesis project for Cranfield University
#   MSc in Computational Fluid Dynamics
#
#   ----------------------------------------------------------------
#     demand -->    --    --    --    --    --    --    supply -->
#   ----------------------------------------------------------------
#
#   demand of road i - the flow (veh/hr) demanded by road i from road i-1
#   supply of road i - the flow (veh/hr) supplied to road i+1 from road i
#
#   network demand/supply is inlet/outlet flow in/out of the sources/sinks
#

import os                    # standard
import numpy as np           # numerical programming
from tqdm import tqdm        # progressbar in time loop
import time                  # recording execution time

# Start time
time0 = time.time()

##### PLAYGROUND #####
playground = False
if playground:
    some_dict={'nest':{'length':5,'speed':4},
               'type':'bank'}
    some_vector = [1,2,3,4,5,6,6]
    print(type(some_dict))
    some_dict = some_vector
    print(type(some_dict))
    exit()
######################

#
# Define Road Characteristics
# Each road must have:
#   length - road segment length (in km)
#   vmax   - max speed (in km/h)
#   source - binary selection if source(1)/not(0)
#   sink   - binary selection if sink(1)/not(0)
#   demand - demand function of density
#   supply - supply function of density
#
# If a road is a source or sink:
#   change source(1)/sink(1) to the function
#   demand_upstream (source) or supply_downstream (sink)
#   a function of time and returns a flow (in veh/hr)
#
# Define a new road by network[road] with road = 1, 2, 3, ..., n
# the network object is a nested dictionary
# network = {road_1: {'length_1': ___, 'vmax_1': ___, 'source_1': ___, 'sink_1': ___, 'demand_1': ___, 'supply_1': ___},
#             ...
#            road_n: {'length_n': ___, 'vmax_n': ___, 'source_n': ___, 'sink_n': ___, 'demand_n': ___, 'supply_n': ___}
#           }
#
network = {}

# Road 1
network[1] = {'length': 5, 'vmax': 90, 'source': 1, 'sink': 0}
def demand(rho): return (90*rho)*(rho <= 30) + 2700*(rho > 30)
def supply(rho): return 2700*(rho <= 30) + (15*(30-rho)+2700)*(rho > 30)
network[1]['demand'] = demand
network[1]['supply'] = supply
def demand_upstream(t): return 0
network[1]['source'] = demand_upstream

# Road 2
network[2] = {'length': 5, 'vmax': 90, 'source': 0, 'sink': 1}
def demand(rho): return (90*rho)*(rho <= 30) + 2700*(rho > 30)
def supply(rho): return 2700*(rho <= 30) + (15*(30-rho)+2700)*(rho > 30)
network[2]['demand'] = demand
network[2]['supply'] = supply
def supply_downstream(t): return 10000
network[2]['sink'] = supply_downstream

# Road 3
network[3] = {'length': 5, 'vmax': 90, 'source': 0, 'sink': 0}
def demand(rho): return (90*rho)*(rho <= 30) + 2700*(rho > 30)
def supply(rho): return 2700*(rho <= 30) + (15*(30-rho)+2700)*(rho > 30)
network[3]['demand'] = demand
network[3]['supply'] = supply

# Road 4
network[4] = {'length': 5, 'vmax': 90, 'source': 0, 'sink': 1}
def demand(rho): return (90*rho)*(rho <= 30) + 2700*(rho > 30)
def supply(rho): return 2700*(rho <= 30) + (15*(30-rho)+2700)*(rho > 30)
network[4]['demand'] = demand
network[4]['supply'] = supply
def supply_downstream(t): return 10000
network[4]['sink'] = supply_downstream

# Road 5
network[5] = {'length': 5, 'vmax': 90, 'source': 0, 'sink': 1}
def demand(rho): return (90*rho)*(rho <= 30) + 2700*(rho > 30)
def supply(rho): return 2700*(rho <= 30) + (15*(30-rho)+2700)*(rho > 30)
network[5]['demand'] = demand
network[5]['supply'] = supply
def supply_downstream(t): return 10000
network[5]['sink'] = supply_downstream

# Define Junction Characteristics
# To unambiguously describe a traffic network, a description of the junctions is required
# Each junction i in (1, n) must be prescribed:
#   in  - a list of the road indexes which feed traffic IN to junction i
#   out - a lsit of the road indexes which lead OUT from junction i
#   tdm - a unique traffic distribution matrix
#
# The Traffic Distribution Matrix, A
#   The general junction with m incoming roads and n outgoing roads
#   A has m rows, n columns
#   The elements A=(a)_mn describe the proportion of traffic leaving incoming road m and travelling on outgoing road n
#   The sum of a_mn over index m is 1, as all traffic must leave its current road and also be conserved
#
junction_info = {}

# Junction 1
junction_info[1] = {'in': [1], 'out': [2, 3], 'tdm': np.array([[0.5, 0.5]])}
junction_info[2] = {'in': [3], 'out': [4, 5], 'tdm': np.array([[0.5, 0.5]])}

# Initialise for road loop
global_flows = {}  # To store each road's junction in/outflow
#
# global_flows is a nested dictionary which stores the supply/demand (or in/outflow) of each road
# non-source/sink roads have supply/demand which is updated at each time step according to the junction function
# source or sink roads have a fixed demand or supply function respectively
# global_flows = {road_1: {'supply': ___, 'demand': ___},
#                  ...
#                 road_n: {'supply': ___, 'demand': ___}
#                }
#
V_max = network[1]['vmax']-1e-3  # Global maximum network speed
total_length = 0  # Total network length
sources = []  # list of source-road indexes
sinks = []  # list of sink-road indexes
# Road loop
for road in network:

    global_flows[road] = {}  # Each road is a dictionary with demand and supply

    V_max = max(V_max, network[road]['vmax'])  # compare all speeds

    total_length += network[road]['length']  # add up all road lengths

    if network[road]['source'] != 0:
        sources.append(road)
    elif network[road]['sink'] != 0:
        sinks.append(road)

# Geometry and junctions (map) defined
time1 = time.time()

# Number of road sections
nb_link = len(network)

# Spatial resolution
dx = 0.2  # km

# Final time
T = 0.5  # hr

# Courant-Friedrichs-Lewy (CFL) safety factor constraint
CFL = 1.5
dt = dx/(CFL*V_max)  # Time resolution from CFL constraint

# Spatial grid
x = np.arange(dx/2, total_length+dx/2, dx)

# Define initial global density profile
Rho_0 = np.zeros(len(x))
for i in range(0, len(x)):
    if x[i] <= 1:
        Rho_0[i] = 100

# Set initial density as the first row of Rho array
n_t = len(np.arange(0, T, dt))
n_x = int(total_length/dx)
Rho = np.zeros((n_t, n_x))
Rho[0, ] = Rho_0

# Godunov solver function
def godunov(network, Rho_0, f_demand_upstream, f_supply_downstream, dx, T):

    # road length
    length = network['length']

    # number of x points
    n_x = int(length / dx)

    rho = np.zeros(len(Rho_0))

    # Get current road supply and demand
    supply = network['supply']
    demand = network['demand']

    # types list for callable function or non-callable number
    allowed_types = (np.float64, int, np.int64)

    # first cell
    j = 0
    temp_demand_upstream = f_demand_upstream if isinstance(f_demand_upstream, allowed_types) else f_demand_upstream((i+1)*dt)
    inflow = min(temp_demand_upstream, supply(Rho_0[j]))
    outflow = min(demand(Rho_0[j]), supply(Rho_0[j+1]))
    rho[j] = Rho_0[j]+(dt/dx)*(inflow - outflow)
    inflow = outflow

    # internal cells
    for j in range(1, n_x-1):
        outflow = min(demand(Rho_0[j]), supply(Rho_0[j+1]))
        rho[j] = Rho_0[j] + (dt/dx)*(inflow - outflow)
        inflow = outflow

    # end cell
    j = n_x-1
    temp_supply_downstream = f_supply_downstream if isinstance(f_supply_downstream, allowed_types) else f_supply_downstream((i+1)*dt)
    outflow = min(demand(Rho_0[j]), temp_supply_downstream)
    rho[j] = Rho_0[j] + (dt/dx)*(inflow - outflow)

    return rho

# Junction flow function
def junction(network,A,rho_0, junction_number):

    # Get number of in and out roads
    n_in = len(junction_info[junction_number]['in'])
    n_out = len(junction_info[junction_number]['out'])

    # Initialise demand(in)/supply(out) array
    sup_dem = np.zeros(len(rho_0))

    # Get demand and supply functions
    # Evaluate at the appropriate boundary density
    # Road IN - demand f'n
    for df_i in range(0, n_in):
        road_identifier = junction_info[junction_number]['in'][df_i]
        sup_dem[df_i] = network[road_identifier]['demand'](rho_0[df_i])

    # Road OUT - supply f'n
    for sf_i in range(0, n_out):
        road_identifier = junction_info[junction_number]['out'][sf_i]
        sup_dem[n_in + sf_i] = network[road_identifier]['supply'](rho_0[n_in + sf_i])

    # Initialise empty flows output array
    flows = np.zeros(len(rho_0))

    # Fill with values
    for flow_i in range(0, len(flows)):
        if flow_i <= n_in-1:
            # These are the supply of in-roads

            # Create supply min list
            sup_min_proportion = sup_dem[flow_i+n_in]/A[flow_i, 0]
            for out_i in range(1, n_out):
                sup_min_proportion = min(sup_min_proportion, sup_dem[out_i+n_in]/A[flow_i, out_i])

            # Calculate in-road new boundary flow
            flows[flow_i] = min(sup_dem[flow_i], sup_min_proportion)
        else:
            # These are the demands of out-roads

            # Construct the sum q_out(j)=sum_j(A(i,j)*q_in(i))
            flow_out = 0
            for in_i in range(0, n_in): # for all i in roads
                flow_out += A[in_i, flow_i-n_in]*flows[in_i]

            # Calculate the out-road new boundary flow
            flows[flow_i] = flow_out

    return flows

# Functions and arrays created
time2 = time.time()

# Time loop
for t in tqdm(np.arange(dt, T, dt)):  # loop with progress bar
# for t in np.arange(dt, T, dt):  # loop without progress bar

    # Iteration index from time t
    i = int(round(t/dt)-1)

    # Loop for each junction
    for junction_index in junction_info:

        # Extract from the junction information dict
        A = junction_info[junction_index]['tdm']
        roads_in = junction_info[junction_index]['in']
        roads_out = junction_info[junction_index]['out']

        # Initialise an empty array for boundary densities
        rho_0 = np.zeros(len(roads_in)+len(roads_out))
        rho_0_index_values = np.zeros(len(roads_in)+len(roads_out))

        # Fill these boundary densities
        #   with the end of in-roads and the start of out-roads
        for rho_0_index in range(0, len(rho_0)):

            # Get road value
            road_identifier = (roads_in+roads_out)[rho_0_index]

            # Get start and end index for this road
            start = 0
            for road_index in range(1, road_identifier):
                start += int(network[road_index]['length'] / dx)

            end = start + int(network[road_identifier]['length'] / dx) - 1

            if rho_0_index <= len(roads_in)-1:
                # first len(roads_in) rho_0 elements are for in-road ends
                rho_0[rho_0_index] = Rho[i, end]
                rho_0_index_values[rho_0_index]=end
                #print('In-Road', road_identifier, 'element', end)
            else:
                # later rho_0 elements are out-road starts
                rho_0[rho_0_index] = Rho[i, start]
                rho_0_index_values[rho_0_index] = start
                #print('Out-Road', road_identifier, 'element', start)

        # Call junction to get input/output flows
        local_flows = junction(network, A, rho_0, junction_index)

        # Store each road's new supply/demand
        for local_flow_i in range(0, len(local_flows)):
            local_flow_val = local_flows[local_flow_i]
            if local_flow_i < len(roads_in):
                road_identifier = junction_info[junction_index]['in'][local_flow_i]
                global_flows[road_identifier]['supply'] = local_flow_val
            else:
                road_identifier = junction_info[junction_index]['out'][local_flow_i-len(roads_in)]
                global_flows[road_identifier]['demand'] = local_flow_val

    # Fill the source/sink supply/demand values in global_flows
    for source in sources:
        global_flows[source]['demand'] = network[source]['source']
    for sink in sinks:
        global_flows[sink]['supply'] = network[sink]['sink']

    # Network global supply/demand definitions
    def net_glob_demand(road, t):

        if road in sources:
            n_g_d = global_flows[road]['demand'](t)
        else:
            n_g_d = global_flows[road]['demand']

        return n_g_d

    def net_glob_supply(road, t):

        if road in sinks:
            n_g_s = global_flows[road]['supply'](t)
        else:
            n_g_s = global_flows[road]['supply']

        return n_g_s

    # Each road separately
    for road in network:

        # Update the new 'initial' density as the previous-time-step solution
        start = 0
        for road_index in range(1, road):
            start += int(network[road_index]['length'] / dx)

        end = start + int(network[road]['length'] / dx)
        rho_0 = Rho[i, start:end]

        # Run Godunov for each road over a single time step
        Rho[i+1, start:end] = godunov(network[road],
                                      rho_0,
                                      net_glob_demand(road, t),
                                      net_glob_supply(road, t),
                                      dx,
                                      dt)

# Time loop completed
time3 = time.time()

# Update the saved density profile in pwd
os.remove('density.txt')
np.savetxt('density.txt', Rho)

# Results saved - program complete
time4 = time.time()

# Print program times
do_print = False
if do_print:
    print('\n')
    print('  Program time statistics')
    print('  -----------------------------------')
    print('  Total program time : {0:.4e} [s]'.format(time4-time0))
    print('  Defining map       : {0:.4e} [s]'.format(time1-time0))
    print('  Initialising       : {0:.4e} [s]'.format(time2-time1))
    print('  Time loop          : {0:.4e} [s]'.format(time3-time2))
    print('  Save results       : {0:.4e} [s]'.format(time4-time3))

print('Done')
exit()
