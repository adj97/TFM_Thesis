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
# network[1] = {'length': 5, 'vmax': 90, 'dmax': 30, 'source': 1, 'sink': 1}
# def demand(rho): return (90*rho)*(rho <= 30) + 2700*(rho > 30)
# def supply(rho): return 2700*(rho <= 30) + (15*(30-rho)+2700)*(rho > 30)
# network[1]['demand'] = demand
# network[1]['supply'] = supply
# def demand_upstream(): return 0
# network[1]['source'] = demand_upstream
# def supply_downstream(): return 10000000000000
# network[1]['sink'] = supply_downstream

# Demand and Supply - density relations

def demand(rho): return (90*rho)*(rho <= 30) + 2700*(rho > 30)
def supply(rho): return 2700*(rho <= 30) + (15*(30-rho)+2700)*(rho > 30)

# Sink condition

def supply_downstream(t): return 1e12

# MOTORWAY ROAD SOURCE/SINKS

network[1] = {'length': 0.5, 'vmax': 112, 'dmax': 208, 'source': 0, 'sink': 1}
network[1]['demand'] = demand
network[1]['supply'] = supply
network[1]['sink'] = supply_downstream

network[2] = {'length': 0.5, 'vmax': 112, 'dmax': 208, 'source': 1, 'sink': 0}
network[2]['demand'] = demand
network[2]['supply'] = supply
def demand_upstream(t): return 0.2
network[2]['source'] = demand_upstream

network[3] = {'length': 0.5, 'vmax': 112, 'dmax': 208, 'source': 1, 'sink': 0}
network[3]['demand'] = demand
network[3]['supply'] = supply
def demand_upstream(t): return 0.1
network[3]['source'] = demand_upstream

network[4] = {'length': 0.5, 'vmax': 112, 'dmax': 208, 'source': 0, 'sink': 1}
network[4]['demand'] = demand
network[4]['supply'] = supply
network[4]['sink'] = supply_downstream

# MOTORWAY INNER SECTION

network[5] = {'length': 0.97, 'vmax': 112, 'dmax': 208, 'source': 0, 'sink': 0}
network[5]['demand'] = demand
network[5]['supply'] = supply

network[6] = {'length': 0.97, 'vmax': 112, 'dmax': 208, 'source': 0, 'sink': 0}
network[6]['demand'] = demand
network[6]['supply'] = supply

# SLIP ROADS

network[7] = {'length': 0.48, 'vmax': 112, 'dmax': 208, 'source': 0, 'sink': 0}
network[7]['demand'] = demand
network[7]['supply'] = supply

network[8] = {'length': 0.32, 'vmax': 112, 'dmax': 208, 'source': 0, 'sink': 0}
network[8]['demand'] = demand
network[8]['supply'] = supply

network[9] = {'length': 0.32, 'vmax': 112, 'dmax': 208, 'source': 0, 'sink': 0}
network[9]['demand'] = demand
network[9]['supply'] = supply

network[10] = {'length': 0.48, 'vmax': 112, 'dmax': 208, 'source': 0, 'sink': 0}
network[10]['demand'] = demand
network[10]['supply'] = supply

# A-ROAD SOURCE/SINKS

network[11] = {'length': 0.2, 'vmax': 80, 'dmax': 208, 'source': 1, 'sink': 0}
network[11]['demand'] = demand
network[11]['supply'] = supply
def demand_upstream(t): return 0.05
network[11]['source'] = demand_upstream

network[12] = {'length': 0.2, 'vmax': 80, 'dmax': 208, 'source': 0, 'sink': 1}
network[12]['demand'] = demand
network[12]['supply'] = supply
network[12]['sink'] = supply_downstream

network[13] = {'length': 0.2, 'vmax': 80, 'dmax': 208, 'source': 0, 'sink': 1}
network[13]['demand'] = demand
network[13]['supply'] = supply
network[13]['sink'] = supply_downstream

network[14] = {'length': 0.2, 'vmax': 80, 'dmax': 208, 'source': 1, 'sink': 0}
network[14]['demand'] = demand
network[14]['supply'] = supply
def demand_upstream(t): return 0.1
network[14]['source'] = demand_upstream

# ROUNDABOUT SECTIONS

network[15] = {'length': 0.13, 'vmax': 64, 'dmax': 208, 'source': 0, 'sink': 0}
network[15]['demand'] = demand
network[15]['supply'] = supply

network[16] = {'length': 0.03, 'vmax': 64, 'dmax': 208, 'source': 0, 'sink': 0}
network[16]['demand'] = demand
network[16]['supply'] = supply

network[17] = {'length': 0.08, 'vmax': 64, 'dmax': 208, 'source': 0, 'sink': 0}
network[17]['demand'] = demand
network[17]['supply'] = supply

network[18] = {'length': 0.03, 'vmax': 64, 'dmax': 208, 'source': 0, 'sink': 0}
network[18]['demand'] = demand
network[18]['supply'] = supply

network[19] = {'length': 0.13, 'vmax': 64, 'dmax': 208, 'source': 0, 'sink': 0}
network[19]['demand'] = demand
network[19]['supply'] = supply

network[20] = {'length': 0.03, 'vmax': 64, 'dmax': 208, 'source': 0, 'sink': 0}
network[20]['demand'] = demand
network[20]['supply'] = supply

network[21] = {'length': 0.08, 'vmax': 64, 'dmax': 208, 'source': 0, 'sink': 0}
network[21]['demand'] = demand
network[21]['supply'] = supply

network[22] = {'length': 0.03, 'vmax': 64, 'dmax': 208, 'source': 0, 'sink': 0}
network[22]['demand'] = demand
network[22]['supply'] = supply

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
# junction_info[id] = {'in': [ list of roads in by integer],
#                     'out': [ list of roads out by integer],
#                     'tdm': np.array([[first row of TDM elements], [second row of TDM elements], ...])}

# Junctions
junction_info[1] = {'in': [5,7], 'out': [1], 'tdm': np.array([[1], [1]])}
junction_info[2] = {'in': [2], 'out': [6,8], 'tdm': np.array([[0.9, 0.1]])}
junction_info[3] = {'in': [3], 'out': [5,9], 'tdm': np.array([[0.8, 0.2]])}
junction_info[4] = {'in': [6,10], 'out': [4], 'tdm': np.array([[1], [1]])}
junction_info[5] = {'in': [22], 'out': [7,15], 'tdm': np.array([[0.8, 0.2]])}
junction_info[6] = {'in': [8,15], 'out': [16], 'tdm': np.array([[1], [1]])}
junction_info[7] = {'in': [16], 'out': [13,17], 'tdm': np.array([[0.8, 0.2]])}
junction_info[8] = {'in': [14,17], 'out': [18], 'tdm': np.array([[1], [1]])}
junction_info[9] = {'in': [18], 'out': [10,19], 'tdm': np.array([[0.7, 0.3]])}
junction_info[10] = {'in': [9,19], 'out': [20], 'tdm': np.array([[1], [1]])}
junction_info[11] = {'in': [20], 'out': [12,21], 'tdm': np.array([[0.7, 0.3]])}
junction_info[12] = {'in': [11,21], 'out': [22], 'tdm': np.array([[1], [1]])}
