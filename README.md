# Traffic Flow Modelling Thesis
Source code for python based traffic flow modelling simulation for Cranfield University individual research project

## Model Capabilities
### - Network
- Any number of junctions 
- Junctions are n in, m out
- Any number of sources and sinks
- Error messages for incorrect usage
- Write simulation run info text file
- Split `main.py` to input map file and parameter file `params.txt`
### - Numerical
- Spatial reconstruction
  - 1st and 2nd order
  - MUSCL 2nd and 3rd order, with 15 slope limiter options
  - WENO 3rd, 5th and 7th (monotonicity preserving bounds) order
- Riemann solvers / numerical flux calculations
  - Lax-Friedrichs
  - Rusanov
  - HLL
  - Murman-Roe
- Classic 4th order Runge-Kutta update scheme

## Ideas / Future Capabilities
- Traffic distribution matrix could be a function of time for different preferences throughout a working day
- Runge-Kutta error adaptive global time step size
- Density gradient adaptive local spatial step size

## Test cases used
1. Example test case description, *`Author, A.N., Paper Name, 2019`*

## File Breakdown
- `.idea` \(folder) : Python project files
- `venv` \(folder) : Python virtual environment libraries
- `MUSCLReconstruction.py` : MUSCL reconstruction functions
- `README.md` : Project and repository information markdown file
- `TFM_Network_Pseudocode.txt` : A general structure plan for the whole program
- `WENOReconstruction.py` : WENO reconstruction functions
- `define_map.py` : User input define road network and junctions (Python)
- `density_plot.m` : Visualise and verify solution/numerical calculations (MATLAB)
- `main.py` : Flow modelling program (Python)
- `params.txt` : Simulation JSON parameter file read in by `main.py`
