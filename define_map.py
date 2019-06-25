import numpy as np           # numerical programming

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

network = {}

# Road Template
# network[1] = {'length': 5, 'vmax': 90, 'source': 1, 'sink': 0}
# def demand(rho): return (90*rho)*(rho <= 30) + 2700*(rho > 30)
# def supply(rho): return 2700*(rho <= 30) + (15*(30-rho)+2700)*(rho > 30)
# network[1]['demand'] = demand
# network[1]['supply'] = supply
# def demand_upstream(t): return 0
# network[1]['source'] = demand_upstream

# SOURCE ROADS

# road 1
network[1] = {'length': 5, 'vmax': 90, 'source': 1, 'sink': 0}
def demand(rho): return (90*rho)*(rho <= 30) + 2700*(rho > 30)
def supply(rho): return 2700*(rho <= 30) + (15*(30-rho)+2700)*(rho > 30)
network[1]['demand'] = demand
network[1]['supply'] = supply
def demand_upstream(t): return 800 + 200*np.sin(80*(t+0.05)) - 50*np.sin(160*(t-0.05))
network[1]['source'] = demand_upstream

# road 2
network[2] = {'length': 5, 'vmax': 90, 'source': 1, 'sink': 0}
def demand(rho): return (90*rho)*(rho <= 30) + 2700*(rho > 30)
def supply(rho): return 2700*(rho <= 30) + (15*(30-rho)+2700)*(rho > 30)
network[2]['demand'] = demand
network[2]['supply'] = supply
def demand_upstream(t): return 600 + 200*np.sin(80*t) - 50*np.sin(160*(t+0.1))
network[2]['source'] = demand_upstream

# INNER ROADS

# road 3
network[3] = {'length': 5, 'vmax': 90, 'source': 0, 'sink': 0}
def demand(rho): return (90*rho)*(rho <= 30) + 2700*(rho > 30)
def supply(rho): return 2700*(rho <= 30) + (15*(30-rho)+2700)*(rho > 30)
network[3]['demand'] = demand
network[3]['supply'] = supply

# Identical inner roads
for i in range(4, 14):
    network[i] = network[3]

# SINK ROADS

# road 14
network[14] = {'length': 5, 'vmax': 90, 'source': 0, 'sink': 1}
def demand(rho): return (90*rho)*(rho <= 30) + 2700*(rho > 30)
def supply(rho): return 2700*(rho <= 30) + (15*(30-rho)+2700)*(rho > 30)
network[14]['demand'] = demand
network[14]['supply'] = supply
def supply_downstream(t): return 100000
network[14]['sink'] = supply_downstream

# road 15
network[15] = {'length': 5, 'vmax': 90, 'source': 0, 'sink': 1}
def demand(rho): return (90*rho)*(rho <= 30) + 2700*(rho > 30)
def supply(rho): return 2700*(rho <= 30) + (15*(30-rho)+2700)*(rho > 30)
network[15]['demand'] = demand
network[15]['supply'] = supply
def supply_downstream(t): return 100000
network[15]['sink'] = supply_downstream

# road 16
network[16] = {'length': 5, 'vmax': 90, 'source': 0, 'sink': 1}
def demand(rho): return (90*rho)*(rho <= 30) + 2700*(rho > 30)
def supply(rho): return 2700*(rho <= 30) + (15*(30-rho)+2700)*(rho > 30)
network[16]['demand'] = demand
network[16]['supply'] = supply
def supply_downstream(t): return 100000
network[16]['sink'] = supply_downstream

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
#   The numpy syntax for arrays in 1D is [a, b, c, ...] so for 2D arrays each of a, b, c, .. are replaced by arrays:
#       [ [a11, a12, a13, ...], [a21, a22, a23, ...], [a31, a32, a33, ...], ... ]

junction_info = {}

# Junction template
# junction_info[1] = {'in': [ list of roads in by integer],
#                     'out': [ list of roads out by integer],
#                     'tdm': np.array([[first row of TDM elements], [second row of TDM elements], ...])}

# Junctions
junction_info[1] = {'in': [1], 'out': [3, 4], 'tdm': np.array([[0.8, 0.2]])}
junction_info[2] = {'in': [2], 'out': [5], 'tdm': np.array([[1]])}
junction_info[3] = {'in': [4, 5], 'out': [6, 8, 9], 'tdm': np.array([[0.5, 0.4, 0.1], [0.3, 0.4, 0.3]])}
junction_info[4] = {'in': [3, 6], 'out': [7], 'tdm': np.array([[1], [1]])}
junction_info[5] = {'in': [8], 'out': [10, 11], 'tdm': np.array([[0.4, 0.6]])}
junction_info[6] = {'in': [9], 'out': [12, 13], 'tdm': np.array([[0.234, 0.766]])}
junction_info[7] = {'in': [7, 10], 'out': [14], 'tdm': np.array([[1], [1]])}
junction_info[8] = {'in': [11, 12], 'out': [15], 'tdm': np.array([[1], [1]])}
junction_info[9] = {'in': [13], 'out': [16], 'tdm': np.array([[1]])}
