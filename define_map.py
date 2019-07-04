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
# def demand_upstream(t): return 0
# network[1]['source'] = demand_upstream
# def supply_downstream(t): return 10000000000000
# network[1]['sink'] = supply_downstream

# SOURCE ROADS

# road 1
network[1] = {'length': 5, 'vmax': 90, 'dmax': 150, 'source': 1, 'sink': 0}
def demand(rho): return (90*rho)*(rho <= 30) + 2700*(rho > 30)
def supply(rho): return 2700*(rho <= 30) + (15*(30-rho)+2700)*(rho > 30)
network[1]['demand'] = demand
network[1]['supply'] = supply
def demand_upstream(t): return 0
network[1]['source'] = demand_upstream

network[2] = {'length': 5, 'vmax': 90, 'dmax': 150, 'source': 0, 'sink': 0}
def demand(rho): return (90*rho)*(rho <= 30) + 2700*(rho > 30)
def supply(rho): return 2700*(rho <= 30) + (15*(30-rho)+2700)*(rho > 30)
network[2]['demand'] = demand
network[2]['supply'] = supply

network[3] = {'length': 5, 'vmax': 90, 'dmax': 150, 'source': 0, 'sink': 0}
def demand(rho): return (90*rho)*(rho <= 30) + 2700*(rho > 30)
def supply(rho): return 2700*(rho <= 30) + (15*(30-rho)+2700)*(rho > 30)
network[3]['demand'] = demand
network[3]['supply'] = supply

network[4] = {'length': 5, 'vmax': 90, 'dmax': 150, 'source': 0, 'sink': 0}
def demand(rho): return (90*rho)*(rho <= 30) + 2700*(rho > 30)
def supply(rho): return 2700*(rho <= 30) + (15*(30-rho)+2700)*(rho > 30)
network[4]['demand'] = demand
network[4]['supply'] = supply

network[5] = {'length': 5, 'vmax': 90, 'dmax': 150, 'source': 0, 'sink': 0}
def demand(rho): return (90*rho)*(rho <= 30) + 2700*(rho > 30)
def supply(rho): return 2700*(rho <= 30) + (15*(30-rho)+2700)*(rho > 30)
network[5]['demand'] = demand
network[5]['supply'] = supply

network[6] = {'length': 5, 'vmax': 90, 'dmax': 150, 'source': 0, 'sink': 1}
def demand(rho): return (90*rho)*(rho <= 30) + 2700*(rho > 30)
def supply(rho): return 2700*(rho <= 30) + (15*(30-rho)+2700)*(rho > 30)
network[6]['demand'] = demand
network[6]['supply'] = supply
def supply_downstream(t): return 100000000
network[6]['sink'] = supply_downstream

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
junction_info[1] = {'in': [1], 'out': [2, 3], 'tdm': np.array([[0.5, 0.5]])}
junction_info[2] = {'in': [2], 'out': [4], 'tdm': np.array([[1]])}
junction_info[3] = {'in': [3], 'out': [5], 'tdm': np.array([[1]])}
junction_info[4] = {'in': [4, 5], 'out': [6], 'tdm': np.array([[1], [1]])}
