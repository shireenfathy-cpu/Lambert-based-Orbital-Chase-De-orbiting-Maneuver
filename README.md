**Lambert-based Orbital Chase & De-orbiting Maneuver** 

## Overview
This project simulates an *orbital chase maneuver* using the Lambert problem, followed by a *de-orbiting maneuver* to remove space debris.  
It includes MATLAB scripts for orbit propagation, Lambert transfer, and visualization of the Earth with orbits.  

## Features
- Solve Lambert problem to compute transfer orbits.  
- Perform orbital chase between two spacecraft.  
- De-orbiting sequence using low-thrust model.  
- 3D visualization of Earth and orbital trajectories.  

## Requirements
- MATLAB R2021a or later.  
- Functions included in this repo:
  - `stateVecFromOE` – Convert orbital elements to state vectors.  
  - `OEFromStateVec` – Convert state vectors to orbital elements.  
  - `solveLambert` – Solve Lambert’s problem.  
  - `solveKepler` – Solve Kepler’s equation.  
  - `getOrbitInputs` – User orbital input interface.  
  - `plotEarthSphere` (optional) – Earth 3D visualization.  
  - `getJulianDate` (optional) – Convert calendar date to Julian date.  

## How to Run
1. Open MATLAB and run the main script:  
   ```matlab
   main.m
