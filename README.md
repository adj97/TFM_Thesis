# Traffic Flow Modelling Thesis
Source code for python based traffic flow modelling simulation for Cranfield University individual research project

## Model Capabilities
- [x] Two single roads joined by a *ghost* junction (03/06/2019)
- [x] One - Two junction (05/06/2019)
- [ ] Any number of junctions 
- [ ] Junctions are n in, m out
- [ ] Any number of sources and sinks
- [ ] \(Optional) Split the source code in `main.py` into multiple files for cleaner easier program

## File Breakdown
- `TFM_Network_Pseudocode.txt` : A general structure plan for the whole program
- `density.txt` : Saved density profile over the whole network (most recently simulated profile)
- `density_plot.m` : Visualise and verify solution/numerical calculations (MATLAB)
- `main.py` : Whole flow modelling program code (Python)
