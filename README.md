# Microgrid Sizing and MPC-based Energy Management System

This project provides a comprehensive, two-stage workflow for designing and operating a solar PV and Battery Energy Storage System (BESS) for Commercial & Industrial (C&I) applications. 
The system architecture is broken down into system sizing via Particle Swarm Optimization (PSO) and hourly operational dispatch using Mixed Integer Linear Programming (MILP) with a Model Predictive Control (MPC) wrapper.

## Project Workflow

* **System Sizing (Python):** The project begins with a PSO algorithm that determines the optimal physical hardware parameters for the microgrid. It simulates operations at 1-minute intervals to find the most cost-effective PV capacity, battery capacity, inverter rating, and grid peak-shaving limit over a 25-year project lifespan.
* **Receding-Horizon Dispatch (MATLAB):** The operational phase uses a Mixed-Integer Linear Programming (MILP) approach to manage the microgrid dynamically. It utilizes the hardware constraints calculated by the PSO phase to minimize daily energy costs, demand charges, and grid reliance.The Receding-Horizon (RH) approach is actualized using the MPC wrapper. The approach,at every hour, solves the 24-hour-ahead MILP using the forecast data, but only keeps the battery command for the very next hour. It then applies that command to the "actual" system, recomputes the real grid draw, updates the real SOC, and moves the clock forward one hour before re-solving.

## Core File Overview

### Sizing Optimization

* **`sizing_pso.py`**: A Python script that utilizes Particle Swarm Optimization (PSO) to output the ideal microgrid specifications. It evaluates initial capital costs, operation and maintenance (O&M), battery/inverter replacement costs, and heavy penalty weights for grid peak violations or curtailment.

### Energy Management System (EMS)

* **`EMS_input_data2.m`**: A MATLAB function that generates 24-hour cyclic data arrays for PV power, C&I load, and tariff pricing. It handles real-world KPLC Time-of-Use (TOU) commercial tariff schedules and establishes a soft-cap penalty threshold for grid draw. Currently, it utilizes two different CI load profile datasets: 'THIWASCO' and 'MAT2'. it also utilizes an equator-based solar irradiance profile of which Kenya is a part of.
* **`EMS_MILP_Design.m`**: The core optimization engine that formulates the equations and tables for steps 1 to 7: Decision variables at every hour of the horizon, power balance with curtailment, battery energy dynamics plus SOC bounds, binary mutual-exclusivity (ensures simultaneous charging and discharging does not occur),the peak-demand epigraph in kVA, the objective function (energy cost + demand charge + penalty),linear penalty mechanism or soft cap on planned grid draw.
* **`MILP_MPCwrapper.m`**: This is step 8: The MPC or Receding Horizon wrapper. It is the main execution script using MATLAB's `intlinprog` solver to simulate both a 72-hour operating window for 'MAT2' dataset and a 168-hour operating window for 'THIWASCO' dataset separately with the use of a switch. It makes use of an 8% forecast error noise to mimic real-world unpredictability and implements the receding-horizon control loop: viewing 24 hours ahead, locking in only the immediate next hour's charge/discharge command, updating the true battery state, and repeating the process.

## Data

The project utilizes two pairs of distinct datasets tailored for the different temporal resolutions required by the sizing and optimization stages:

* **`pvLoadCostData.mat` or **: Contains **1-minute** time-series data. This high-resolution dataset is used by the sizing algorithm (`sizing_pso.py`) to accurately capture rapid fluctuations in load and irradiance for optimal hardware sizing.
* **`pvLoadPriceData.mat` or THIWASCO_real_load_Apr19to25_2026**: Contains **hourly** time-series data. This dataset includes arrays for clear, cloudy, and partly cloudy irradiance day-types, alongside cost and variable load profiles. It is used by the MATLAB MPC scripts to calculate hourly operational dispatch.

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
3. **Run the Simulation:** Execute `MILP_MPCwrapper.m` in MATLAB. The script will output a comparison against a baseline (grid-only) setup, detailing total cost reductions and generating operational plots for either the 'MAT2' load profile or the 'THIWASCO' load profile. Finally, it exports a `.mat` file formatted for Simscape physical validation.
