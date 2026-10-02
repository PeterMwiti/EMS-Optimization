% Step 8:MPC/receding-horizon wrapper
%
% What this script does, in plain terms:
%   1. Builds an "actual" weather/load/tariff timeline and a slightly
%      noisier "forecast" version of the same timeline (real forecasts
%      are never perfect -- this is what makes the receding-horizon
%      re-solving worth doing).
%   2. At every hour, solves the 24-hour-ahead MILP (EMS_MILP_Design.m)
%      using the FORECAST data, but only keeps the battery command for
%      the very next hour.
%   3. Applies that command to the "actual" system, recomputes the real
%      grid draw from the real power balance, updates the real SOC, and
%      moves the clock forward one hour before re-solving.
%   4. Compares the resulting cost against the grid-only baseline
%      (Section 3.6) and prints/plots the results.

% clear; clc; close all;
clearvars; clc; close all;

% %% ---------------- configuration ----------------
% p = struct();
% p.dt        = 1;          % timestep (h)
% p.eta_ch    = 0.95;        % charge efficiency (ηch)
% p.eta_dis   = 0.95;        % discharge efficiency (ηdis)
% p.Pch_max   = 250;         % kW, max charge power (Pch(t))
% p.Pdis_max  = 250;         % kW, max discharge power (Pdis(t))
% p.Ecap      = 1000;        % kWh, nominal battery capacity
% p.Emin      = 0.20*p.Ecap; % SOC floor 20% (200 kWh)
% p.Emax      = 0.90*p.Ecap; % SOC ceiling 90% (900 kWh)
% 
% %CI1; (Consumption threshold: >15,000 KWh/Month,
%        % Peak TOU: 13.44 KES/kWh, Off-peak TOU: 6.72 KES/kWh
%        % Demand charge: 1,100 KES/kVA)
% 
% p.Pgrid_max = 2000;        % kW, contracted/physical grid import limit (Pgrid(t)) 
% p.CD        = 1100;        % KES/kVA, demand charge rate (CD)
% p.rho       = 15;          % KES/kWh, linear penalty weight on overshoot s(t)
%                             % (should sit above lambda's peak rate (13.44 KES/kWh) so the
%                             % EMS prefers staying under Pcap when it can,
%                             % but never so high the problem behaves like
%                             % a hard, potentially-infeasible constraint)
% 
% data_params = struct('Ppv_rated', 500, 'Pload_base', 150, 'Pload_peak', 600);
% 
% T_horizon      = 24;   % hours the MILP looks ahead at each solve
% Nsim           = 72;   % hours to actually simulate (3 days)
% forecast_noise = 0.08; % 8% forecast error, similar order to the MAPE
%                         % used for LMP forecasting in Wang et al. (2024)

%% ---------------- configuration ----------------
p = struct();
p.dt        = 1;          % timestep (h)
p.eta_ch    = 0.95;       % charge efficiency (ηch)
p.eta_dis   = 0.95;       % discharge efficiency (ηdis)

% Inverter Capacity (from PSO: 612.23 kW)
p.Pch_max   = 612.23;     % kW, max charge power (Pch(t))
p.Pdis_max  = 612.23;     % kW, max discharge power (Pdis(t))

% Battery Capacity (from PSO: 5000.00 kWh)
p.Ecap      = 5000;       % kWh, nominal battery capacity
p.Emin      = 0.20*p.Ecap;% SOC floor 20% (1000 kWh)
p.Emax      = 0.90*p.Ecap;% SOC ceiling 90% (4500 kWh)

% CI1 Tariff Configuration & Physical Grid Limit
p.Pgrid_max = 2000;       % kW, contracted/physical grid import limit (Pgrid(t)) 
p.CD        = 1100;       % KES/kVA, demand charge rate (CD)
p.rho       = 15;         % KES/kWh, linear penalty weight on overshoot s(t)

% Microgrid Sizing Parameters (from PSO: 1360.61 kWp PV / 695 kW Peak Load)
data_params = struct('Ppv_rated', 1360.61, 'Pload_base',350, 'Pload_peak', 895);

T_horizon      = 24;   % hours the MILP looks ahead at each solve
Nsim           = 72;   % hours to actually simulate (3 days)
forecast_noise = 0.08; % 8% forecast error

%% ---------------- build actual and forecast data ----------------
N_total = Nsim + T_horizon - 1; %(covers the full lookahead window needed by the last hour in the simulation))
[Ppv_actual, Pload_actual, lambda_actual, Pcap_actual] = ...
    EMS_input_data(N_total, data_params);

rng(42);  % reproducible demo noise (ensures that we always get the same random noise patterns in every run)

%Roughen the true PV and load values by about 8% hour by hour and ensure
%they don't go to a negative value...this simulates the unpredictable load
%and pv patterns;
Ppv_forecast   = max(Ppv_actual   .* (1 + forecast_noise*randn(N_total,1)), 0);
Pload_forecast = max(Pload_actual .* (1 + forecast_noise*randn(N_total,1)), 0);

% Tariff and the soft-cap rule are treated as known in advance (they are
% published schedules, not forecasted quantities).
lambda_forecast = lambda_actual;
Pcap_forecast   = Pcap_actual;

%% ---------------- receding-horizon loop ----------------
E_actual   = zeros(Nsim+1,1); %(blank notebook on the battery's actual energy level per hour)
E_actual(1) = 0.5*(p.Emin+p.Emax);   % start at mid-SOC

Pgrid_actual = zeros(Nsim,1); %(blank notebook on the grid power pulled per hour)
Pch_cmd_log  = zeros(Nsim,1); %(blank notebook on the battery charging command issued per hour)
Pdis_cmd_log = zeros(Nsim,1); %(blank notebook on the battery discharging command issued per hour)
curtailed    = zeros(Nsim,1); %(blank notebook on how much surplus PV power had to be thrown away per hour)

% Ppeak_running = 0; %(Initial peak demand before running)
Ppeak_running = 576.00; % Initial target grid peak from PSO sizing (kW)

fprintf('Running receding-horizon MILP EMS for %d hours ...\n', Nsim);
for k = 1:Nsim %(from 1 to 72)
    window = k : k+T_horizon-1; %(k:k+24-1; Makes the MILP look 24hrs ahead at every time step)

    sol = EMS_MILP_Design( ...
        Ppv_forecast(window), Pload_forecast(window), ...
        lambda_forecast(window), Pcap_forecast(window), ...
        E_actual(k), Ppeak_running, p);

 % Flags NaN values in the run as errors while stating at what timestep
 % they occured
    if any(isnan(sol.Pgrid))
        error('EMS solve failed at hour %d -- see warning above.', k);
    end

    Pch_cmd  = sol.Pch(1); %pulls out the first entry(charging power) of that 24 hour long list(sol.Pch)
    Pdis_cmd = sol.Pdis(1); %pulls out the first entry(discharging power) of that 24 hour long list(sol.Pdis)


    % --- real power balance using the ACTUAL (not forecast) PV/load ---
    Pgrid_needed = Pload_actual(k) - Ppv_actual(k) - Pdis_cmd + Pch_cmd;
    if Pgrid_needed < 0
        % More generation than needed after the battery command -- with
        % no export/net-metering modeled yet (see Section 2.2.4), surplus
        % is curtailed. Swap this block out if/when export credit is added.
        curtailed(k) = -Pgrid_needed;
        Pgrid_needed = 0;
    end
    Pgrid_actual(k) = Pgrid_needed;

    % --- update the real battery state using the committed command ---
    E_next = E_actual(k) + Pch_cmd*p.eta_ch*p.dt - (Pdis_cmd/p.eta_dis)*p.dt;
    E_actual(k+1) = min(max(E_next, p.Emin), p.Emax);

    Ppeak_running = max(Ppeak_running, Pgrid_actual(k));

    Pch_cmd_log(k)  = Pch_cmd;
    Pdis_cmd_log(k) = Pdis_cmd;
end
fprintf('Done.\n\n');

%% ---------------- cost accounting (proposed system) ----------------
lambda_win = lambda_actual(1:Nsim);
energy_cost_ems = sum(lambda_win .* Pgrid_actual .* p.dt);
demand_cost_ems = p.CD * Ppeak_running;
total_cost_ems  = energy_cost_ems + demand_cost_ems;

%% ---------------- baseline (grid-only) ----------------
base = baseline_EMS_dispatch(Pload_actual(1:Nsim), lambda_win, p.CD, p.dt);

%% ---------------- summary metrics (Section 3.6) ----------------
pct_reduction = 100*(base.total_cost - total_cost_ems)/base.total_cost;
par_ems  = max(Pgrid_actual)/mean(Pgrid_actual);
par_base = max(base.Pgrid)/mean(base.Pgrid);
served_pv_bess = 1 - sum(Pgrid_actual)/sum(Pload_actual(1:Nsim));

fprintf('=================== RESULTS (%d-hour window) ===================\n', Nsim);
fprintf('%-32s %14s %14s\n','','Baseline','Proposed EMS');
fprintf('%-32s %14.1f %14.1f\n','Energy cost (KES)', base.energy_cost, energy_cost_ems);
fprintf('%-32s %14.1f %14.1f\n','Demand charge (KES)', base.demand_cost, demand_cost_ems);
fprintf('%-32s %14.1f %14.1f\n','Total cost (KES)', base.total_cost, total_cost_ems);
fprintf('%-32s %14.1f %14.1f\n','Peak grid demand (kW)', base.peak, Ppeak_running);
fprintf('%-32s %14.2f %14.2f\n','Peak-to-average ratio', par_base, par_ems);
fprintf('------------------------------------------------------------------\n');
fprintf('Cost reduction vs. baseline: %.1f%%\n', pct_reduction);
fprintf('Share of load served by PV+BESS: %.1f%%\n', 100*served_pv_bess);
fprintf('===================================================================\n');

%% ---------------- plots ----------------
t = (1:Nsim)';
figure('Name','EMS Dispatch');
subplot(3,1,1);
plot(t, base.Pgrid, 'g--', t, Pgrid_actual, 'b-', 'LineWidth', 1.3);
yline(Ppeak_running,'b:','EMS peak'); yline(base.peak,'g:','Baseline peak');
legend('Baseline grid draw','EMS grid draw','Location','best');
ylabel('kW'); title('Grid Power Drawn'); grid on;

subplot(3,1,2);
plot(t, Pch_cmd_log, 'g-', t, -Pdis_cmd_log, 'r-', 'LineWidth', 1.3);
legend('Charge','Discharge (negative)','Location','best');
ylabel('kW'); title('Battery Charge/Discharge Command'); grid on;

subplot(3,1,3);
plot(t, 100*E_actual(1:Nsim)/p.Ecap, 'm-', 'LineWidth', 1.3);
yline(100*p.Emin/p.Ecap,'k:'); yline(100*p.Emax/p.Ecap,'k:');
ylabel('SOC (%)'); xlabel('Hour'); title('Battery State of Charge'); grid on;

%% ---------------- export for Simscape (Section 3.5) ----------------
% Package the realized dispatch as timeseries objects that a
% "From Workspace" block can feed directly into the Simscape PV/BESS
% model for the physical feasibility/safety validation stage.
ts_Pgrid = timeseries(Pgrid_actual, t*3600, 'Name','Pgrid_EMS');
ts_Pch   = timeseries(Pch_cmd_log,  t*3600, 'Name','Pch_cmd');
ts_Pdis  = timeseries(Pdis_cmd_log, t*3600, 'Name','Pdis_cmd');
save('EMS_dispatch_for_simscape.mat', 'ts_Pgrid','ts_Pch','ts_Pdis');
fprintf('Saved EMS_dispatch_for_simscape.mat for use in the Simscape model.\n');