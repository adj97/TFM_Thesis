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

from __future__ import print_function    # print no new line
import os                                # standard
import datetime
import numpy as np                       # numerical programming
from tqdm import tqdm                    # progressbar in time loop
import time                              # recording execution time
import json                              # read json format parameter file
import define_map                        # separate python file for network and junction info
import WENOReconstruction                # |----------"-----------| WENO reconstructions
import MUSCLReconstruction                # |----------"-----------| WENO reconstructions

# Start time
time0 = time.time()

# DEVELOPMENT
# PLAYGROUND
playground = False
if playground:
    pass
    exit()
##############

# From define_map.py file
#   read network and junction info
network = define_map.network
junction_info = define_map.junction_info

# Check for input errors in define_map.py
error_info = {'mis_attr': 0,
              'attr_typ': 0,
              'tdm_row_sum': 0,
              'tdm_shape': 0,
              'rd_id_er': 0}
error_description = {'mis_attr': 'missing attribute',
                     'attr_typ': 'wrong attribute type',
                     'tdm_row_sum': 'TDM row sum error',
                     'tdm_shape': 'TDM shape error',
                     'rd_id_er': 'road-junction ID mismatch'}
error_message = ''
jn_inout = {'in': [], 'out': []}

for junction in junction_info:

    # Check all junction attributes are present
    junction_keys = list(junction_info[junction].keys())
    junction_keys.sort()
    if junction_keys != ['in', 'out', 'tdm']:
        error_message += 'Badly defined junction {} \n'.format(junction)
        error_info['mis_attr'] += junction_keys != ['in', 'out', 'tdm']

    # Check junction info data types
    for in_road in junction_info[junction]['in']+junction_info[junction]['out']:  # in/out roads
        if (not isinstance(in_road, int)) or in_road > len(network) or in_road <= 0:
            error_info['attr_typ'] += 1
            error_message += 'Badly defined in/out road ID in junction {} \n'.format(junction)

    allowed_types = (int, float)
    for row in junction_info[junction]['tdm'].tolist():  # tdm elements
        for a in row:
                if not isinstance(a, allowed_types):
                    error_info['attr_typ'] += 1
                    error_message += 'Badly defined TDM element, {}, in junction {} \n'.format(a, junction)
        # Check row sum
        if abs(1 - sum(row)) >= 1e-3:
            error_info['tdm_row_sum'] += 1
            error_message += 'Incorrect TDM row sum, in junction {} \n'.format(junction)

    # Check TDM matches in/out roads
    shape = junction_info[junction]['tdm'].shape
    junction_roads = [len(junction_info[junction]['in']), len(junction_info[junction]['out'])]
    if (shape[0]!=junction_roads[0]) or (shape[1]!=junction_roads[1]):
        error_info['tdm_shape'] += 1
        error_message += 'Mismatch TDM shape and in/out roads, in junction {} \n'.format(junction)

    # In/Out list
    for in_id in junction_info[junction]['in']:
        jn_inout['in'].append(in_id)
    for out_id in junction_info[junction]['out']:
        jn_inout['out'].append(out_id)

for road in network:

    #Check all road attributes are present
    network_keys = list(network[road].keys())
    network_keys.sort()
    if network_keys != ['demand', 'dmax', 'length', 'sink', 'source', 'supply', 'vmax']:
        error_message += 'Badly defined road {} \n'.format(road)
        error_info['mis_attr'] += network_keys != ['demand', 'length', 'sink', 'source', 'supply', 'vmax']

    # Check correct road attribute types
    for attr in network_keys:
        if attr in ['demand', 'supply']:
            # Supply and demand should be functions of density
            if not callable(network[road][attr]):
                error_info['attr_typ'] += 1
                error_message += 'Badly defined {} in road {} \n'.format(attr, road)

        elif attr in ['length', 'vmax', 'dmax']:
            # Length and vmax should be real numbers
            allowed_types = (int, float)
            if (not isinstance(network[road][attr], allowed_types)) or network[road][attr] <= 0:
                error_info['attr_typ'] += 1
                error_message += 'Badly defined {} in road {} \n'.format(attr, road)

        elif attr in ['source', 'sink']:
            # Source and sink values can be either 0 or a function
            if not (network[road][attr]==0 or callable(network[road][attr])):
                error_info['attr_typ'] += 1
                error_message += 'Badly defined {} in road {} \n'.format(attr, road)

    # Check road-junction index matching
    if len(network)>1:
        if network[road]['source'] != 0:
            # Source roads only go IN to junctions
            if road in jn_inout['out']:
                error_info['rd_id_er'] += 1
                error_message += 'Bad ID for source road {}'.format(road)

            # Source road goes in to ONE junction only
            if jn_inout['in'].count(road) != 1:
                error_info['rd_id_er'] += 1
                error_message += 'Bad ID for source road {}'.format(road)

        elif network[road]['sink'] != 0:
            # Sink roads only come OUT of junctions
            if road in jn_inout['in']:
                error_info['rd_id_er'] += 1
                error_message += 'Bad ID for sink road {}'.format(road)

            # Sink road comes out of ONE junction only
            if jn_inout['out'].count(road) != 1:
                error_info['rd_id_er'] += 1
                error_message += 'Bad ID for sink road {}'.format(road)

        else:
            # Inner roads should be on both sides
            if jn_inout['in'].count(road)+jn_inout['out'].count(road) < 2:
                error_info['rd_id_er'] += 1
                error_message += 'Bad ID for internal road {}'.format(road)

# Get total number of all errors
error_total = 0
for error_type in error_info:
    error_total += error_info[error_type]

# Print input error and exit program
if error_total >= 1:
    print('ERROR : Found', error_total, 'error(s)')
    print('\nError breakdown:')
    for err in error_info:
        if error_info[err] != 0:
            print('{} {}'.format(error_info[err], error_description[err]))  # Breakdown
    print('\nError message(s):')
    print(error_message)  # Messages
    exit(1)
else:
    # No errors
    pass


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
total_length = 0                 # Total network length
sources = []                     # list of source-road indexes
sinks = []                       # list of sink-road indexes
# Road loop
for road in network:

    global_flows[road] = {}  # Each road is a dictionary with demand and supply

    V_max = max(V_max, network[road]['vmax'])  # compare all speeds

    total_length += network[road]['length']  # add up all road lengths

    # Create source and sink lists
    if network[road]['source'] != 0: sources.append(road)
    if network[road]['sink'] != 0: sinks.append(road)

# Geometry and junctions (map) defined, errors checked
time1 = time.time()

# Number of road sections
nb_link = len(network)

# Read parameter file
with open('params.txt') as file:
    params = json.load(file)

# Assign parameters
dx = params['dx']                             # Spatial resolution [km]
T = params['T']                               # Final time [hr]
CFL = params['CFL']                           # Courant-Friedrichs-Lewy (CFL) safety factor constraint
reconstruction = params['reconstruction']     # Reconstruction scheme

# Time resolution from CFL constraint
dt = dx/(CFL*V_max)  # [hr]

# Spatial grid
x = np.arange(dx/2, total_length+dx/2, dx)

# Initial density profile
Rho_0 = np.zeros(len(x))
#### INPUT ####
# Define initial global density profile
for i in range(0, len(x)):
    if x[i] < 1:
        Rho_0[i] = 10

# Set initial density as the first row of Rho array
n_t = len(np.arange(0, T, dt))
n_x = int(total_length/dx)
Rho = np.zeros((n_t, n_x))
Rho[0, ] = Rho_0

# WENO reconstruction definitions
weno3 = WENOReconstruction.weno3
weno5 = WENOReconstruction.weno5
weno7 = WENOReconstruction.weno7

# MUSCL reconstruction definitions
muscl2 = MUSCLReconstruction.muscl2
muscl3 = MUSCLReconstruction.muscl3

# Riemann solvers


# Specific timers
junction_time = 0
reconstruction_time = 0
RK_update_time = 0

# Get road start and end index function
def get_start_end(road_id):
    fn_start = 0
    for fn_road_index in range(1, road_id):
        fn_start += int(network[fn_road_index]['length'] / dx)
    fn_end = fn_start + int(network[road_id]['length'] / dx) - 1

    return [fn_start, fn_end]

# Minmod function
minmod = WENOReconstruction.minmod

# Junction flow function
def junction(network_section, A, rho_0, junction_number):

    # Get number of in and out roads
    n_in = len(junction_info[junction_number]['in'])
    n_out = len(junction_info[junction_number]['out'])

    # Initialise demand(in)/supply(out) array
    sup_dem = np.zeros(len(rho_0))

    # Get demand and supply functions
    # Evaluate at the appropriate boundary density
    # Only need the demand of last cells on in roads
    #           and supply of first cells on out roads
    # Road IN - demand f'n
    for df_i in range(0, n_in):
        road_identifier = junction_info[junction_number]['in'][df_i]
        sup_dem[df_i] = network_section[road_identifier]['demand'](rho_0[df_i])

    # Road OUT - supply f'n
    for sf_i in range(0, n_out):
        road_identifier = junction_info[junction_number]['out'][sf_i]
        sup_dem[n_in + sf_i] = network_section[road_identifier]['supply'](rho_0[n_in + sf_i])

    # Initialise empty flows output array
    flows = np.zeros(len(rho_0))

    # Fill with values
    for flow_i in range(0, len(flows)):
        # Fill each flows element

        if flow_i <= n_in-1:
            # These are the supply of in-roads

            # Smallest demands of all outgoing roads
            sup_min_proportion = sup_dem[n_in]/A[flow_i, 0]
            for out_i in range(1, n_out):
                sup_min_proportion = min(sup_min_proportion, sup_dem[out_i+n_in]/A[flow_i, out_i])

            # Calculate in-road new boundary flow
            flows[flow_i] = min(sup_dem[flow_i], sup_min_proportion)
        else:
            # These are the demands of out-roads

            # Sum of partial flow contribution from each in road
            flows[flow_i] = 0
            for in_i in range(0, n_in): # for all i in roads
                flows[flow_i] += A[in_i, flow_i-n_in]*flows[in_i]

    return flows

# Velocity model function u=u(rho)
def velocity_model(road_section, density):

    # Greenshield's
    u_f = road_section['vmax']
    d_m = road_section['dmax']
    velocity_value = u_f*(1-density/d_m)

    return velocity_value

# Compute the numerical flux
def compute_flux(CdL, CdR, FlL, FlR):

    # CdL/CdR - conserved left/right density
    # FlL/FlR - physical left/right flux

    riemann_solver = 'LaxFriedrichs'

    if riemann_solver == 'LaxFriedrichs':

        # Lax-Friedrichs
        flux_value = 0.5*(FlL+FlR)-0.5*(dx/dt)*(CdR-CdL)

    elif riemann_solver == 'Rusanov':

        # Rusanov
        print('Rusanov Solver')
        exit()

    elif riemann_solver == 'Murman-Roe':

        # Murman-Roe
        print('Murman-Roe Solver')
        exit()

    elif riemann_solver == 'HLL':

        # HLL
        print('HLL Solver')
        exit()

    else:  # Wrong flux specifier

        print('ERROR: Wrong flux specifier')
        exit(1)


    return flux_value

# Spatial reconstruction
def spatial_reco(reconstruction_type):

    # Create array
    output_array = np.zeros((len(CellFluxes)+1, 2))

    # Loop over all cells
    for cell_id in range(0, len(output_array)):

        if reconstruction_type == 'FirstOrder':
            # First order reconstruction
            left = rho_ghost[cell_id]
            right = rho_ghost[cell_id]

        elif reconstruction_type == 'SecondOrder':
            # 2nd order minmod reconstruction

            # Right
            xmr = rho_ghost[cell_id] - rho_ghost[cell_id - 1]
            ymr = rho_ghost[cell_id + 1] - rho_ghost[cell_id]
            right = rho_ghost[cell_id] + 0.5 * minmod(xmr, ymr)

            # Left
            xml = rho_ghost[cell_id] - rho_ghost[cell_id-1]
            yml = rho_ghost[cell_id + 1] - rho_ghost[cell_id]
            left = rho_ghost[cell_id] - 0.5 * minmod(xml, yml)

        elif reconstruction_type == 'WENO3':
            # 5th-Order Weighted Essentially Non-Oscillatory reconstruction

            [left, right] = weno3(cell_id, rho_ghost)

        elif reconstruction_type == 'WENO5':
            # 5th-Order Weighted Essentially Non-Oscillatory reconstruction

            [left, right] = weno5(cell_id, rho_ghost)

        elif reconstruction_type =='WENO7':
            # 7th-Order Weighted Essentially Non-Oscillatory reconstruction

            [left, right] = weno7(cell_id, rho_ghost)

        elif reconstruction_type == 'MUSCL2':
            # 2nd-Order Monotonic Upwind reconstruction Scheme for Conservation Laws

            [left, right] = muscl2(cell_id, rho_ghost)

        elif reconstruction_type == 'MUSCL3':
            # 3rd-Order Monotonic Upwind reconstruction Scheme for Conservation Laws

            [left, right] = muscl3(cell_id, rho_ghost)

        else:  # Wrong reconstruction specifier

            print('ERROR: Wrong reconstruction specifier')
            exit(1)

        output_array[cell_id] = [left, right]
    return output_array

# Network global demand
def net_glob_demand(road, t):

    if road in sources:
        n_g_d = global_flows[road]['demand'](t)
    else:
        n_g_d = global_flows[road]['demand']

    return n_g_d

# Network global supply
def net_glob_supply(road, t):

    if road in sinks:
        n_g_s = global_flows[road]['supply'](t)
    else:
        n_g_s = global_flows[road]['supply']

    return n_g_s

# Functions and arrays created
time2 = time.time()

# Time loop
for t in tqdm(np.arange(dt, T, dt)):  # loop with progress bar
# for t in np.arange(dt, T, dt):          # loop without progress bar

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

        # Fill these boundary densities
        #   with the end of in-roads and the start of out-roads
        for rho_0_index in range(0, len(rho_0)):

            # Get road value
            road_identifier = (roads_in+roads_out)[rho_0_index]

            # Get road start and end indexes
            [start, end] = get_start_end(road_identifier)

            # extract appropriate boundary elements
            if rho_0_index <= len(roads_in)-1:
                # in-road end
                rho_0[rho_0_index] = Rho[i, end]
            else:
                # out-road start
                rho_0[rho_0_index] = Rho[i, start]

        # Log junction time
        junction_time -= time.time()

        # Call junction to get input/output flows
        local_flows = junction(network, A, rho_0, junction_index)

        # Log junction time
        junction_time += time.time()

        # Store each road's new supply/demand in global_flows
        for local_flow_i in range(0, len(local_flows)):
            local_flow_val = local_flows[local_flow_i]

            if local_flow_i < len(roads_in):
                # in-road
                road_identifier = junction_info[junction_index]['in'][local_flow_i]
                global_flows[road_identifier]['supply'] = local_flow_val
            else:
                # out-road
                road_identifier = junction_info[junction_index]['out'][local_flow_i-len(roads_in)]
                global_flows[road_identifier]['demand'] = local_flow_val

    # Fill the network source/sink supply/demand values in global_flows
    for source in sources:
        global_flows[source]['demand'] = network[source]['source']
    for sink in sinks:
        global_flows[sink]['supply'] = network[sink]['sink']

    # Each road separately
    for road in network:

        # Get start and end indexes
        [start, end] = get_start_end(road)

        # Update the new 'initial' density as the previous-time-step solution
        rho_0 = Rho[i, start:end+1]

        # number of x points
        n_x = int(network[road]['length'] / dx)

        # initialise density array
        rho1 = Rho[i + 1, start:end + 1]
        rho2 = Rho[i + 1, start:end + 1]
        rho3 = Rho[i + 1, start:end + 1]

        # Get current road supply and demand
        supply = network[road]['supply']
        demand = network[road]['demand']

        # types list for callable function or non-callable number
        allowed_types = (np.float64, int, np.int64)

        # incoming flow
        f_demand_upstream = net_glob_demand(road, t)
        f_supply_downstream = net_glob_supply(road, t)

        # Runge-Kutta steps
        for RK in range(0, 3):

            # Chose the most updated array for flux calculation
            if RK == 0:
                rho_flux = rho_0
            elif RK == 1:
                rho_flux = rho1
            elif RK == 2:
                rho_flux = rho2

            # Add in ghost BCs
            rho_ghost = np.append(rho_flux, [rho_flux[-2:-5:-1],   # end ghosts
                                             rho_flux[3:0:-1]])    # start ghosts

            # Initialise cell flux array
            CellFluxes = np.zeros(len(rho_0))

            # Log reconstruction time
            reconstruction_time -= time.time()

            # Save all reconstructed L and R states in an array size len(CellFluxes) X 2
            reconstructed = spatial_reco(reconstruction)

            # Log reconstruction time
            reconstruction_time += time.time()

            # Compute each cell RHS interface flux
            for cell_id in range(0,len(CellFluxes)):

                # Reconstruction : Get L and R conserved(C) densities
                CL = reconstructed[cell_id, 1]
                CR = reconstructed[cell_id+1, 0]

                # Get L and R physical flux(F)
                FL = CL*velocity_model(network[road], CL)
                FR = CR*velocity_model(network[road], CR)

                # Calculate flux
                CellFluxes[cell_id] = compute_flux(CL, CR, FL, FR)

            # Update

            # Log Runge-Kutta update time
            RK_update_time -= time.time()

            # first cell
            j = 0
            temp_demand_upstream = f_demand_upstream if isinstance(f_demand_upstream, allowed_types) else f_demand_upstream((i+1)*dt)
            use = temp_demand_upstream - CellFluxes[j]

            if RK == 0:
                rho1[j] = rho_0[j] + (dt / dx) * use
            elif RK == 1:
                rho2[j] = 0.75 * rho_0[j] + 0.25 * rho1[j] + 0.25 * (dt / dx) * use
            elif RK == 2:
                rho3[j] = (1 / 3) * rho_0[j] + (2 / 3) * rho2[j] + (2 / 3) * (dt / dx) * use
            if road in sources:
                rho1[j] = global_flows[road]['demand'](0)
                rho2[j] = global_flows[road]['demand'](0)
                rho3[j] = global_flows[road]['demand'](0)

            # internal cells
            for j in range(1, n_x-1):
                use = CellFluxes[j - 1] - CellFluxes[j]
                if RK == 0:
                    rho1[j] = rho_0[j] + (dt / dx) * use
                elif RK == 1:
                    rho2[j] = 0.75*rho_0[j] + 0.25*rho1[j] + 0.25*(dt / dx) * use
                elif RK == 2:
                    rho3[j] = (1/3) * rho_0[j] + (2/3) * rho2[j] + (2/3) * (dt / dx) * use

            # last cell
            j = n_x-1
            temp_supply_downstream = f_supply_downstream if isinstance(f_supply_downstream, allowed_types) else f_supply_downstream((i + 1) * dt)
            outflow = min(CellFluxes[j], f_supply_downstream)
            use = CellFluxes[j - 1] - outflow
            if RK == 0:
                rho1[j] = rho_0[j] + (dt / dx) * use
            elif RK == 1:
                rho2[j] = 0.75 * rho_0[j] + 0.25 * rho1[j] + 0.25 * (dt / dx) * use
            elif RK == 2:
                rho3[j] = (1 / 3) * rho_0[j] + (2 / 3) * rho2[j] + (2 / 3) * (dt / dx) * use

            # Log Runge-Kutta update time
            RK_update_time += time.time()


        Rho[i+1, start:end+1] = rho3

# Time loop completed
time3 = time.time()

# Chose to save or not
do_print = True

# Save density profile
if do_print:
    # Create output folder
    if not os.path.exists('Simulation_Results'): os.mkdir('Simulation_Results')
    path = 'Simulation_Results/'+str(datetime.datetime.now().strftime("%d-%m-%Y_%H-%M-%S"))
    os.mkdir(path)

    # Update the saved density profile in pwd
    density_output = path + '/density.txt'
    np.savetxt(density_output, Rho)
else:
    # Developing runs
    os.remove('density.txt')
    np.savetxt('density.txt', Rho)

# Results saved - program complete
time4 = time.time()

# Print program times
if do_print:
    # Overwrite file (comment out)
    filename = path + '/simulation_info.txt'

    # Write and open file
    info_txt = open(filename, 'w+')
    info_txt.write('TFM Simulation Information \n\n')
    now = datetime.datetime.now()
    info_txt.write(now.strftime("%d/%m/%Y %H:%M:%S \n"))

    # Read parameter file
    with open('params.txt') as file:
        params = json.load(file)

    # Find required slope limiter
    chosen_limiter = params['limiter']
    print_limiter = False  # default

    # Write reconstruction phrases
    if reconstruction == 'FirstOrder':
        reconstruction = '1st Order single cell average'
    elif reconstruction == 'SecondOrder':
        reconstruction = '2nd Order minmod interpolation'
    elif reconstruction == 'WENO3':
        reconstruction = '3rd Order WENO'
    elif reconstruction == 'WENO5':
        reconstruction = '5th Order WENO'
    elif reconstruction == 'WENO7':
        reconstruction = '7th Order WENO monotonicity preserving bounds'
    elif reconstruction == 'MUSCL2':
        reconstruction = '2nd Order MUSCL'
        print_limiter = True
    elif reconstruction == 'MUSCL3':
        reconstruction = '3rd Order MUSCL'
        print_limiter = True

    # Method information
    info_txt.write('\n\nMethodology:\n')
    info_txt.write('{:30s} {:.6e}\n'.format('Spatial step', dx))
    info_txt.write('{:30s} {:.6e}\n'.format('Final time', T))
    info_txt.write('{:30s} {:.6e}\n'.format('CFL constraint', CFL))
    info_txt.write('{:30s} {:10s}\n'.format('Reconstruction method', reconstruction))
    if print_limiter: info_txt.write('{:30s} {:10s}\n'.format('Slope limiter', chosen_limiter))

    # Time breakdown
    info_txt.write('\n\nTime breakdown:\n')
    info_txt.write('{:30s} {:.10s}\n'.format('Code Section', 'Time [s]'))
    info_txt.write('----------------------------------------------\n')
    info_txt.write('{:30s} {:.6e}\n'.format('Defining map and error check', time1-time0))
    info_txt.write('{:30s} {:.6e}\n'.format('Initialising', time2-time1))
    info_txt.write('{:30s} {:.6e}\n'.format('Time loop', time3-time2))
    info_txt.write('- - - - - - - - - - - - - - - - - - - - - - - \n')
    info_txt.write('{:30s} {:.6e}\n'.format('Junction solver', junction_time))
    info_txt.write('{:30s} {:.6e}\n'.format('Spatial reconstruction', reconstruction_time))
    info_txt.write('{:30s} {:.6e}\n'.format('Runge-Kutta updates', RK_update_time))
    info_txt.write('- - - - - - - - - - - - - - - - - - - - - - - \n')
    info_txt.write('{:30s} {:.6e}\n'.format('Save results', time4-time3))
    info_txt.write('----------------------------------------------\n')
    info_txt.write('----------------------------------------------\n')
    info_txt.write('{:30s} {:.6e}\n'.format('Total program time', time4-time0))
    info_txt.write('----------------------------------------------\n')

    # File code line count
    num_lines = {}
    files = ['main.py', 'define_map.py', 'params.txt', 'MUSCLReconstruction.py', 'WENOReconstruction.py']
    for file in files:
        num_lines[file] = 0
        with open(file, 'r') as f:
            for line in f:
                num_lines[file] += 1
    num_lines['Total'] = 0

    info_txt.write('\n\nLine Count:\n')
    info_txt.write('{: <30} {: <20}\n'.format('File', 'Number of Lines'))
    info_txt.write('----------------------------------------------\n')
    for file in num_lines:
        if file=='Total':
            info_txt.write('----------------------------------------------\n')
            info_txt.write('{: <30} {: <20}\n'.format(file, num_lines[file]))
            info_txt.write('----------------------------------------------\n')
        else:
            info_txt.write('{: <30} {: <20}\n'.format(file, num_lines[file]))
        num_lines['Total'] += num_lines[file]

    # Write memory information
    info_txt.write('\n\nOutput File Memory:\n')
    info_txt.write('{:30} {:}\n'.format('File', 'Size [MB]'))
    info_txt.write('----------------------------------------------\n')
    info_txt.write('{0:30} {1:5.3f}\n'.format('density.txt',os.stat('density.txt').st_size / 1000000))
    info_txt.write('----------------------------------------------\n')

    # Close file
    info_txt.close()

time.sleep(1)
print('Done')
exit()
