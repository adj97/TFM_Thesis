# Traffic Flow Modelling Thesis
Source code for python based traffic flow modelling simulation for Cranfield University individual research project

## Model Capabilities
- [x] Two single roads joined by a *ghost* junction \(03/06/2019)
- [x] One - Two junction \(05/06/2019)
- [x] Create a parameter file `params.txt` to be read in `main.py` \(17/06/2019)
- [x] Add error messages for incorrect usage \(19/06/2019)
- [x] Tested on a diverge and merge network - < > - \(19/06/2019)
- [x] Write simulation run info text file \(19/06/2019)
- [x] Split `main.py` to input map file \(19/06/2019)
- [ ] Any number of junctions 
- [ ] Junctions are n in, m out
- [ ] Any number of sources and sinks
- [ ] Add extra numerical schemes/methods/limiters/reconstruction

## Ideas \ Future Capabilities
- Traffic distribution matrix could be a function of time for different preferences throughout a working day
- Split output object to hold each road density seperate

## Test cases used
1. Example test case description, *`Author, A.N., Paper Name, 2019`*

## File Breakdown
- `.idea` \(folder) : Python project files
- `Presentation 1 13:06` \(folder) : Presentation slides and material 
- `venv` \(folder) : Python virtual environment libraries
- `README.md` : Project and repository information markdown file
- `TFM_Network_Pseudocode.txt` : A general structure plan for the whole program
- `define_map.py` : User input define road network and junctions (Python)
- `density.txt` : Saved density profile over the whole network (most recently simulated profile)
- `density_plot.m` : Visualise and verify solution/numerical calculations (MATLAB)
- `main.py` : Flow modelling program (Python)
- `params.txt` : Simulation JSON parameter file read in by `main.py`
- `simulaion_info.txt` : Automatically written simulation run information output file
