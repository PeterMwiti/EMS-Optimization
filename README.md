# Microgrid Sizing and MPC-based Energy Management System

This project provides a comprehensive, two-stage workflow for designing and operating a solar PV and Battery Energy Storage System (BESS) for Commercial & Industrial (C&I) applications. 

The system architecture is broken down into system sizing via Particle Swarm Optimization (PSO) and hourly operational dispatch using Model Predictive Control (MPC).

## Project Workflow

1. **System Sizing (Python)**: The project begins with a PSO algorithm that determines the optimal physical hardware parameters for the microgrid. It simulates operations at 1-minute intervals to find the most cost-effective PV capacity, battery capacity, inverter rating, and grid peak-shaving limit over a 25-year project lifespan.
2. **Receding-Horizon Dispatch (MATLAB)**: The operational phase uses a Mixed-Integer Linear Programming (MILP) approach to manage the microgrid dynamically. It utilizes the hardware constraints calculated by the PSO phase to minimize daily energy costs, demand charges, and grid reliance.

## Core File Overview

### Sizing Optimization
* **`sizing_pso.py`**: A Python script that utilizes Particle Swarm Optimization (PSO) to output the ideal microgrid specifications. It evaluates initial capital costs, operation and maintenance (O&M), battery/inverter replacement costs, and heavy penalty weights for grid peak violations or curtailment.

### Energy Management System (EMS)
* **`EMS_input_data.m`**: A MATLAB function that generates 24-hour cyclic data arrays for PV power, C&I load, and tariff pricing. It handles real-world KPLC Time-of-Use (TOU) commercial tariff schedules and establishes a soft-cap penalty threshold for grid draw.
* **`EMS_MILP_Design.m`**: The core optimization engine that formulates the power balance, battery dynamics, and mutual-exclusivity constraints (preventing simultaneous charging and discharging). It solves a 24-hour lookahead step using MATLAB's `intlinprog` solver to minimize costs.
* **`MILP_MPCwrapper.m`**: The main execution script simulating a 72-hour operating window with an 8% forecast error noise to mimic real-world unpredictability. It implements the receding-horizon control loop: viewing 24 hours ahead, locking in only the immediate next hour's charge/discharge command, updating the true battery state, and repeating the process.

### Data
* ** `pvLoadPriceData_New.mat`**: Binary MATLAB data files containing the time-series arrays for clear, cloudy, and partly cloudy irradiance day-types, alongside cost and variable load profiles.

## Dependencies

**Python Requirements:**
* `numpy`
* `scipy`
* `matplotlib`

**MATLAB Requirements:**
* MATLAB R2024b 
* Optimization Toolbox (required for `intlinprog`)

## How to Use

1. **Run the Sizing Optimizer:** Execute `sizing_pso.py` to calculate the optimal PV (kWp), Battery (kWh), and Inverter (kW) capacities. The script will generate a visual plot (`optimization_results.png`) showing power flow and state-of-charge dynamics.
2. **Update the EMS:** Take the capacity outputs from the Python script and update the configuration parameters (`p.Pch_max`, `p.Ecap`, `data_params`) in `MILP_MPCwrapper.m`.
3. **Run the Simulation:** Execute `MILP_MPCwrapper.m` in MATLAB. The script will output a comparison against a baseline (grid-only) setup, detailing total cost reductions and generating operational plots for the 72-hour window. Finally, it exports a `.mat` file formatted for Simscape physical validation.
