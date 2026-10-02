import numpy as np
import scipy.io
import matplotlib.pyplot as plt

def load_data():
    data = scipy.io.loadmat('pvLoadPriceData_New.mat')
    # load profile 2 is the variable load
    # Added to a 350kW constant load
    # Values in loadData seem to be in Watts
    load_profile_2 = data['loadData'][:, 1]
    total_load_kw = (350000 + load_profile_2) / 1000.0
    
    # Solar irradiance of clear day
    irradiance = data['clearDay'].flatten() # W/m^2
    
    # time is 1441 points, assuming 1 minute interval
    return total_load_kw, irradiance

def simulate_system(pv_cap, bat_cap, inv_cap, grid_limit, load_profile, irradiance_profile):
    """
    Simulates the PV-BESS system for 24 hours.
    pv_cap: kWp
    bat_cap: kWh
    inv_cap: kW (max charge/discharge)
    grid_limit: kW (peak shaving threshold)
    """
    dt = 1/60.0 # 1 minute in hours
    n = len(load_profile)
    
    # Efficiencies
    eta_conv = 0.95
    eta_ch = 0.95
    eta_dis = 0.95
    
    # SoC limits
    soc_min = 0.15 * bat_cap
    soc_max = 0.90 * bat_cap
    soc = soc_min # Start at min SoC
    
    p_grid = np.zeros(n)
    p_pv_gen = (pv_cap * irradiance_profile / 1000.0) * eta_conv # kW AC equivalent
    
    soc_history = np.zeros(n)
    curtailed_pv = 0
    
    for t in range(n):
        p_load = load_profile[t]
        p_pv = p_pv_gen[t]
        
        # Net power before battery
        p_net = p_load - p_pv
        
        p_bat = 0 # Power from battery (+ is discharging, - is charging)
        
        if p_net > grid_limit:
            # Need to shave peak
            required_from_bat = p_net - grid_limit
            # Limit by inverter and battery capacity
            available_dis = min(inv_cap, (soc - soc_min) * eta_dis / dt)
            p_bat = min(required_from_bat, available_dis)
            soc -= (p_bat / eta_dis) * dt
        elif p_net > 0:
            # Load exceeds PV but below grid limit
            # Priority: Maximize self-consumption
            required_from_bat = p_net
            available_dis = min(inv_cap, (soc - soc_min) * eta_dis / dt)
            p_bat = min(required_from_bat, available_dis)
            soc -= (p_bat / eta_dis) * dt
        elif p_net < 0:
            # PV exceeds load, charge battery
            excess_pv = -p_net
            available_ch = min(inv_cap, (soc_max - soc) / (eta_ch * dt))
            p_ch_actual = min(excess_pv, available_ch)
            soc += p_ch_actual * eta_ch * dt
            curtailed_pv += (excess_pv - p_ch_actual) * dt
            p_bat = -p_ch_actual
            
        p_grid[t] = p_load - p_pv - p_bat
        soc_history[t] = soc
        
    return p_grid, curtailed_pv, soc_history

def objective_function(x, load_profile, irradiance_profile):
    pv_cap, bat_cap, inv_cap, grid_limit = x
    
    if pv_cap < 0 or bat_cap < 0 or inv_cap < 0 or grid_limit < 0:
        return 1e18
    
    p_grid, curtailed, soc_hist = simulate_system(pv_cap, bat_cap, inv_cap, grid_limit, load_profile, irradiance_profile)
    
    # 25 year period Costs in KSH
    
    # Initial Capital Costs
    c_cap_pv = 65000 * pv_cap
    c_cap_bat = 55000 * bat_cap
    c_cap_inv = 19000 * inv_cap
    tic = c_cap_pv + c_cap_bat + c_cap_inv
    
    # O&M Costs (25 years)
    c_om_pv = 1300 * pv_cap * 25
    c_om_bat = 1000 * bat_cap * 25
    c_om_inv = 400 * inv_cap * 25
    tom = c_om_pv + c_om_bat + c_om_inv
    
    # Replacement Costs (one replacement cycle)
    c_rep_bat = 30000 * bat_cap
    c_rep_inv = 13000 * inv_cap
    trc = c_rep_bat + c_rep_inv
    
    # Grid Costs (25 years)
    daily_grid_energy = np.sum(np.maximum(0, p_grid)) / 60.0
    tgc = daily_grid_energy * 28.24 * 365 * 25
    
    # Penalties
    # 1. Peak shaving penalty: if p_grid exceeds grid_limit
    peak_violation = np.maximum(0, p_grid - grid_limit)
    penalty_peak = np.sum(peak_violation) * 1e7 # Heavy penalty
    
    # 2. Grid dependence penalty
    penalty_grid_reliance = daily_grid_energy * 100 * 365 * 25 
    
    # 3. Self-consumption reward (via curtailment penalty)
    penalty_curtailment = curtailed * 1000 * 365 * 25 

    # 4. Minimize the peak limit itself (to "shave as much as possible")
    # We want the smallest constant peak limit.
    cost_peak_limit = grid_limit * 1e6 # Weight to favor lower limits
    
    total_cost = tic + tom + trc + tgc + penalty_peak + penalty_grid_reliance + penalty_curtailment + cost_peak_limit
    
    return total_cost

def pso(load_profile, irradiance_profile):
    n_particles = 50
    max_iter = 1000
    dim = 4 # pv_cap, bat_cap, inv_cap, grid_limit
    
    # Bounds: [pv_cap, bat_cap, inv_cap, grid_limit]
    lb = np.array([0, 0, 0, 0])
    ub = np.array([2000, 5000, 1000, 800])
    
    pos = np.random.uniform(lb, ub, (n_particles, dim))
    vel = np.zeros((n_particles, dim))
    
    pbest_pos = pos.copy()
    pbest_val = np.array([objective_function(p, load_profile, irradiance_profile) for p in pos])
    
    gbest_idx = np.argmin(pbest_val)
    gbest_pos = pbest_pos[gbest_idx].copy()
    gbest_val = pbest_val[gbest_idx]
    
    w = 0.7
    c1 = 1.5
    c2 = 1.5
    
    print(f"Starting PSO... Initial GBest: {gbest_val:.2e}")
    
    for i in range(max_iter):
        r1, r2 = np.random.rand(n_particles, dim), np.random.rand(n_particles, dim)
        vel = w * vel + c1 * r1 * (pbest_pos - pos) + c2 * r2 * (gbest_pos - pos)
        pos += vel
        pos = np.clip(pos, lb, ub)
        
        for p_idx in range(n_particles):
            val = objective_function(pos[p_idx], load_profile, irradiance_profile)
            if val < pbest_val[p_idx]:
                pbest_val[p_idx] = val
                pbest_pos[p_idx] = pos[p_idx].copy()
                if val < gbest_val:
                    gbest_val = val
                    gbest_pos = pos[p_idx].copy()
        
        if (i+1) % 100 == 0:
            print(f"Iteration {i+1}/{max_iter}, GBest: {gbest_val:.2e}")
            
    return gbest_pos

if __name__ == "__main__":
    load_profile, irradiance_profile = load_data()
    best_sizing = pso(load_profile, irradiance_profile)
    
    pv_cap, bat_cap, inv_cap, grid_limit = best_sizing
    
    pv_efficiency = 0.21
    pv_area = pv_cap / (1.0 * pv_efficiency)
    
    print("\n" + "="*35)
    print("OPTIMIZATION RESULTS (FINAL)")
    print("="*35)
    print(f"PV Capacity:      {pv_cap:10.2f} kWp")
    print(f"PV Panel Area:     {pv_area:10.2f} m^2")
    print(f"Battery Capacity:  {bat_cap:10.2f} kWh")
    print(f"Inverter Capacity: {inv_cap:10.2f} kW")
    print(f"Grid Peak Limit:   {grid_limit:10.2f} kW")
    print("-" * 35)
    
    p_grid, curtailed, soc_hist = simulate_system(pv_cap, bat_cap, inv_cap, grid_limit, load_profile, irradiance_profile)
    daily_load_energy = np.sum(load_profile) / 60.0
    daily_grid_energy = np.sum(np.maximum(0, p_grid)) / 60.0
    sc_ratio = ((daily_load_energy - daily_grid_energy) / daily_load_energy) * 100
    
    print(f"Daily Load Energy: {daily_load_energy:10.2f} kWh")
    print(f"Daily Grid Energy: {daily_grid_energy:10.2f} kWh")
    print(f"Self-Consumption:  {sc_ratio:10.2f} %")
    print(f"Max Grid Peak:     {np.max(p_grid):10.2f} kW")
    print("="*35)

    time_axis = np.linspace(0, 24, len(load_profile))
    plt.figure(figsize=(12, 10))
    
    plt.subplot(3, 1, 1)
    plt.plot(time_axis, load_profile, label='Total Load', color='black')
    plt.plot(time_axis, (pv_cap * irradiance_profile / 1000.0) * 0.95, label='PV Output (AC)', color='orange')
    plt.plot(time_axis, p_grid, label='Grid Power', color='blue', linestyle='--')
    plt.axhline(y=grid_limit, color='red', linestyle=':', label='Grid Limit Target')
    plt.ylabel('Power (kW)')
    plt.legend()
    plt.title('Power Profiles')
    
    plt.subplot(3, 1, 2)
    plt.fill_between(time_axis, 0, soc_hist, color='green', alpha=0.3)
    plt.plot(time_axis, soc_hist, color='green')
    plt.ylabel('Battery SoC (kWh)')
    plt.title('Battery State of Charge')
    
    plt.subplot(3, 1, 3)
    p_bat = load_profile - (pv_cap * irradiance_profile / 1000.0) * 0.95 - p_grid
    plt.plot(time_axis, p_bat, color='purple')
    plt.ylabel('Battery Power (kW)')
    plt.axhline(0, color='black', lw=0.5)
    plt.xlabel('Time (Hours)')
    plt.title('Battery Power Flow (+: Discharging, -: Charging)')
    
    plt.tight_layout()
    plt.savefig('optimization_results.png')
    print("Results plot saved as 'optimization_results.png'")
