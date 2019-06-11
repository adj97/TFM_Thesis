# Traffic Flow Modelling Thesis
Source code for python based traffic flow modelling simulation for Cranfield University individual research project

## Model Capabilities
- [x] Two single roads joined by a *ghost* junction (03/06/2019)
- [x] One - Two junction (05/06/2019)
- [ ] Any number of junctions 
- [ ] Junctions are n in, m out
- [ ] Any number of sources and sinks
- [ ] Create a parameter file `params.py` to be read in `main.py` (T, dx, methodology choices, etc.)
- [ ] Add error messages for incorrect usage
- [ ] Add extra numerical schemes/methods/limiters/reconstruction

### Ideas \(Future Capabilities)
- Traffic distribution matrix could be a function of time for different preferences throughout a working day

## Test cases used
1. Example test case description, *`Author, A.N., Paper Name, 2019`*

## File Breakdown
- `TFM_Network_Pseudocode.txt` : A general structure plan for the whole program
- `density.txt` : Saved density profile over the whole network (most recently simulated profile)
- `density_plot.m` : Visualise and verify solution/numerical calculations (MATLAB)
- `main.py` : Whole flow modelling program code (Python)
