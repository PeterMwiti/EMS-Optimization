% MILP_MPCwrapper.m
% Receding-horizon MILP EMS -- works with EITHER the THIWASCO real-load
% dataset OR the pvLoadPriceData MAT-file dataset (MILP_2), switched
% with ONE line below (data_params.data_source). Everything else in
% this script, and every function it calls, is identical either way.
%
% WHAT THIS MERGES FROM THE TWO ORIGINAL PIPELINES:
%   - From MILP_THIWASCO: curtailment handling, the off-peak battery
%     lockout, kVA peak tracking with Ppeak_floor carry-forward, the
%     cost-comparison bar charts, and the solar-vs-grid power plot.
%   - From MILP_2 (test2MPCwrapper.m): real GTI-based PV generation
%     (now used for BOTH datasets), the solar
%     curtailment tracking/plot, and PSO-informed default sizing.
%   - data_params.data_source is the single switch point. PF is
%     now one scalar (p.PF) used consistently everywhere, instead of
%     THIWASCO's per-hour vector or MAT2's split kW/kVA handling.
%
% What this script does, in plain terms:
%   1. Loads load/PV/tariff data for the selected dataset, plus a
%      slightly noisier "forecast" version (real forecasts are never
%      perfect -- this is what makes the receding-horizon re-solving
%      worth doing).
%   2. At every hour, solves the 24-hour-ahead MILP using the FORECAST
%      data, but only keeps the battery command for the very next hour.
%   3. Applies that command to the "actual" system, recomputes the real
%      grid draw, updates the real SOC, and moves the clock forward one
%      hour before re-solving.
%   4. Compares the resulting cost against the grid-only baseline, and
%      confirms the mutual-exclusivity (SCD) guarantee held throughout.

clear; clc; close all;

%% ================= THE SWITCH: pick a dataset =================
data_source = 'MAT2';   % 'THIWASCO' or 'MAT2'
% ======================================================================

%% ---------------- configuration ----------------
p = struct();
p.PF        = 0.92;        % scalar power factor, used everywhere now (Step 5)
p.dt        = 1;           % timestep (h)
p.eta_ch    = 0.95;        % charge efficiency (ηch)
p.eta_dis   = 0.95;        % discharge efficiency (ηdis)
p.Ecap      = 1000;        % kWh, nominal battery capacity
p.Emin      = 0.20*p.Ecap; % SOC floor 20%
p.Emax      = 0.90*p.Ecap; % SOC ceiling 90%
p.Pgrid_max = 2000;        % kW, contracted/physical grid import limit
p.CD        = 1100;        % KES/kVA, demand charge rate (CI1)
p.rho       = 15;          % KES/kWh, linear penalty weight on overshoot s(t)
                            % (should sit above the tariff's peak rate so the
                            % EMS prefers staying under Pcap when it can, but
                            % never so high the problem risks infeasibility)

% Battery/inverter sizing and dataset-specific defaults. Nsim/T_horizon
% differ by source: THIWASCO is a FIXED real 168-hour week (no
% extrapolation beyond real data), while MAT2's profiles are cyclic
% 24-hour templates and can run for any length.
switch upper(data_source)
  case 'THIWASCO'
    p.Pch_max  = 5000;   p.Pdis_max = 5000;
    data_params = struct('Ppv_rated', 1457.6, 'Pload_base', 350, 'Pload_peak', 692.98, ...
                          'data_source', 'THIWASCO');
    T_horizon = 24;
    Nsim      = 145;    % Nsim + T_horizon - 1 = 168, the full real week
    Ppeak_initial = 0;  % start the running peak from scratch

  case 'MAT2'
    p.Pch_max  = 612.23;  p.Pdis_max = 612.23;   % PSO-sized inverter
    p.Ecap     = 5000;    p.Emin = 0.20*p.Ecap;  p.Emax = 0.90*p.Ecap;  % PSO-sized battery
    data_params = struct('Ppv_rated', 1360.61, 'Pload_base', 350, 'Pload_peak', 695, ...
                          'data_source', 'MAT2', 'pv_day_type', 'clear');
    T_horizon = 24;
    Nsim      = 72;     % 3 simulated days; the cyclic templates support any length
    Ppeak_initial = 576; % PSO-derived target peak used as the starting floor (kW-equivalent)

  otherwise
    error('MILP_MPCwrapper: bad_source', 'data_source must be ''THIWASCO'' or ''MAT2''.');
end

forecast_noise = 0.08; % 8% forecast error, similar order to the MAPE
                        % used for LMP forecasting in Wang et al. (2024)

%% ---------------- build actual and forecast data ----------------
N_total = Nsim + T_horizon - 1;
[Ppv_actual, Pload_actual, lambda_actual, Pcap_actual, Papparent_actual] = ...
    EMS_input_data2(N_total, data_params);

fprintf('Dataset: %s | N_total=%d hours | Ppv_rated=%.1f kWp | Pch_max=%.1f kW | Ecap=%.0f kWh\n', ...
    data_source, N_total, data_params.Ppv_rated, p.Pch_max, p.Ecap);

rng(42);  % reproducible demo noise
Ppv_forecast   = max(Ppv_actual   .* (1 + forecast_noise*randn(N_total,1)), 0);
Pload_forecast = max(Pload_actual .* (1 + forecast_noise*randn(N_total,1)), 0);
% Tariff and the soft-cap rule are treated as known in advance (both are
% published/derived schedules, not forecasted quantities), regardless
% of data source.
lambda_forecast = lambda_actual;
Pcap_forecast   = Pcap_actual;

%% ---------------- receding-horizon loop ----------------
E_actual   = zeros(Nsim+1,1);
E_actual(1) = 0.5*(p.Emin+p.Emax);

Pgrid_actual   = zeros(Nsim,1);   % kW, active power from grid
Pgrid_kVA_log  = zeros(Nsim,1);   % kVA, apparent power from grid
Pch_cmd_log    = zeros(Nsim,1);
Pdis_cmd_log   = zeros(Nsim,1);
curtailed      = zeros(Nsim,1);   % kW, REAL curtailed solar each hour
ych_log        = zeros(Nsim,1);
ydis_log       = zeros(Nsim,1);

Ppeak_running = Ppeak_initial;

fprintf('Running receding-horizon MILP EMS for %d hours ...\n', Nsim);
for k = 1:Nsim
    window = k : k+T_horizon-1;

    sol = EMS_MILP_Design( ...
        Ppv_forecast(window), Pload_forecast(window), ...
        lambda_forecast(window), Pcap_forecast(window), ...
        E_actual(k), Ppeak_running, p);

    if any(isnan(sol.Pgrid))
        error('EMS solve failed at hour %d -- see warning above.', k);
    end

    Pch_cmd  = sol.Pch(1);
    Pdis_cmd = sol.Pdis(1);
    ych_log(k)  = sol.ych(1);
    ydis_log(k) = sol.ydis(1);

    % --- real power balance using the ACTUAL (not forecast) PV/load ---
    Pgrid_needed = Pload_actual(k) - Ppv_actual(k) - Pdis_cmd + Pch_cmd;
    if Pgrid_needed < 0
        % Real surplus beyond what charging absorbed -- no export/net-
        % metering modeled, so it's curtailed (thrown away).
        curtailed(k) = -Pgrid_needed;
        Pgrid_needed = 0;
    end
    Pgrid_actual(k) = Pgrid_needed;
    Pgrid_kVA_log(k) = Pgrid_needed / p.PF;

    % --- update the real battery state using the committed command ---
    E_next = E_actual(k) + Pch_cmd*p.eta_ch*p.dt - (Pdis_cmd/p.eta_dis)*p.dt;
    E_actual(k+1) = min(max(E_next, p.Emin), p.Emax);

    Ppeak_running = max(Ppeak_running, Pgrid_kVA_log(k));

    Pch_cmd_log(k)  = Pch_cmd;
    Pdis_cmd_log(k) = Pdis_cmd;
end
fprintf('Done.\n\n');

%% ---------------- SCD / mutual-exclusivity confirmation (Step 4) ----------------
n_idle    = sum(ych_log == 0 & ydis_log == 0);
n_charge  = sum(ych_log == 1 & ydis_log == 0);
n_dis     = sum(ych_log == 0 & ydis_log == 1);
n_both    = sum(ych_log == 1 & ydis_log == 1);   % must always be 0

fprintf('=================== SCD / MUTUAL-EXCLUSIVITY CHECK ===================\n');
fprintf('Idle hours (both 0, battery doing nothing):      %d\n', n_idle);
fprintf('Charging hours (ych=1, ydis=0):                  %d\n', n_charge);
fprintf('Discharging hours (ych=0, ydis=1):                %d\n', n_dis);
fprintf('Simultaneous charge+discharge hours (must be 0): %d\n', n_both);
assert(n_both == 0, 'MILP_MPCwrapper_unified:SCD_violation', ...
    'Simultaneous charging and discharging occurred -- this should be impossible.');
fprintf('CONFIRMED: charge and discharge were never both 1 in the same hour.\n');
fprintf('========================================================================\n\n');

%% ---------------- cost accounting (proposed system) ----------------
lambda_win = lambda_actual(1:Nsim);
energy_cost_ems = sum(lambda_win .* Pgrid_actual .* p.dt);
demand_cost_ems = p.CD * Ppeak_running;
total_cost_ems  = energy_cost_ems + demand_cost_ems;

%% ---------------- baseline (grid-only) ----------------
base = baseline_EMS_dispatch(Pload_actual(1:Nsim), lambda_win, p.CD, p.dt, ...
    p.PF, Papparent_actual);
if base.peak_is_measured
    fprintf('Baseline peak uses REAL measured apparent power (%s).\n', data_source);
else
    fprintf('Baseline peak is a PF-based ESTIMATE (%s has no measured kVA column).\n', data_source);
end

%% ---------------- summary metrics ----------------
pct_reduction = 100*(base.total_cost - total_cost_ems)/base.total_cost;
par_ems  = max(Pgrid_kVA_log)/mean(Pgrid_kVA_log);
par_base = base.peak / (base.energy_cost/sum(lambda_win*p.dt));  % rough PAR fallback if needed
par_base = max(base.Pgrid) / mean(base.Pgrid);                   % kW-based PAR, source-agnostic
served_pv_bess = 1 - sum(Pgrid_actual)/sum(Pload_actual(1:Nsim));

total_pv_kwh        = sum(Ppv_actual(1:Nsim)) * p.dt;
total_curtailed_kwh = sum(curtailed) * p.dt;
if total_pv_kwh > 0
    solar_utilization_rate = 100 * (1 - total_curtailed_kwh/total_pv_kwh);
else
    solar_utilization_rate = 100;
end

fprintf('=================== RESULTS (%d-hour window, %s) ===================\n', Nsim, data_source);
fprintf('%-32s %14s %14s\n','','Baseline','Proposed EMS');
fprintf('%-32s %14.1f %14.1f\n','Energy cost (KES)', base.energy_cost, energy_cost_ems);
fprintf('%-32s %14.1f %14.1f\n','Demand charge (KES)', base.demand_cost, demand_cost_ems);
fprintf('%-32s %14.1f %14.1f\n','Total cost (KES)', base.total_cost, total_cost_ems);
fprintf('%-32s %14.1f %14.1f\n','Peak grid demand (kVA)', base.peak, Ppeak_running);
fprintf('%-32s %14.2f %14.2f\n','Peak-to-average ratio (kW)', par_base, par_ems);
fprintf('----------------------------------------------------------------------------------\n');
fprintf('Cost reduction vs. baseline: %.1f%%\n', pct_reduction);
fprintf('Share of load served by PV+BESS: %.1f%%\n', 100*served_pv_bess);
fprintf('Total available solar generation: %.1f kWh\n', total_pv_kwh);
fprintf('Total curtailed (wasted) solar:   %.1f kWh\n', total_curtailed_kwh);
fprintf('Solar energy utilization rate:    %.1f%%\n', solar_utilization_rate);
fprintf('====================================================================================\n');

%% ---------------- plots ----------------
t = (1:Nsim)';

figure('Name','EMS Dispatch and Curtailment');
subplot(4,1,1);
plot(t, Pgrid_kVA_log, 'b-', 'LineWidth', 1.3);
yline(Ppeak_running,'b:','EMS peak (kVA)'); yline(base.peak,'g:','Baseline peak (kVA)');
ylabel('kVA'); title('Grid Apparent Power Drawn (demand-charge basis)'); grid on;

subplot(4,1,2);
plot(t, Pch_cmd_log, 'g-', t, -Pdis_cmd_log, 'r-', 'LineWidth', 1.3);
legend('Charge','Discharge (negative)','Location','best');
ylabel('kW'); title('Battery Charge/Discharge Command'); grid on;

subplot(4,1,3);
plot(t, 100*E_actual(1:Nsim)/p.Ecap, 'm-', 'LineWidth', 1.3);
yline(100*p.Emin/p.Ecap,'k:'); yline(100*p.Emax/p.Ecap,'k:');
ylabel('SOC (%)'); title('Battery State of Charge'); grid on;

subplot(4,1,4);
area(t, curtailed, 'FaceColor', [0.85 0.33 0.10], 'FaceAlpha', 0.6, 'EdgeColor', 'r');
ylabel('kW'); xlabel('Hour');
title(sprintf('Solar Curtailment (Total Wasted: %.1f kWh | Utilization: %.1f%%)', ...
      total_curtailed_kwh, solar_utilization_rate));
grid on;

% ---- cost comparison bar charts (Baseline vs Proposed EMS) ----
figure('Name','Cost Comparison');
cats = categorical({'Baseline','Proposed EMS'}, {'Baseline','Proposed EMS'});

subplot(1,2,1);
bar(cats, [base.energy_cost, energy_cost_ems]);
ylabel('KES'); title('Energy Cost Comparison'); grid on;

subplot(1,2,2);
bar(cats, [base.total_cost, total_cost_ems]);
ylabel('KES'); title('Total Cost Comparison'); grid on;

% ---- solar vs grid active power over the simulated window ----
figure('Name','Solar vs Grid Power');
plot(t, Ppv_actual(1:Nsim), 'Color', [0.93 0.69 0.13], 'LineWidth', 1.3); hold on;
plot(t, Pgrid_actual, 'b-', 'LineWidth', 1.3);
legend('Solar (PV) power','Grid (active) power','Location','best');
xlabel('Hour'); ylabel('kW');
title(sprintf('Solar vs. Grid Active Power (%s)', data_source)); grid on;

%% ---------------- export for Simscape (Section 3.5) ----------------
ts_Pgrid = timeseries(Pgrid_actual, t*3600, 'Name','Pgrid_EMS');
ts_Pch   = timeseries(Pch_cmd_log,  t*3600, 'Name','Pch_cmd');
ts_Pdis  = timeseries(Pdis_cmd_log, t*3600, 'Name','Pdis_cmd');
save('EMS_dispatch_for_simscape.mat', 'ts_Pgrid','ts_Pch','ts_Pdis');
fprintf('Saved EMS_dispatch_for_simscape.mat for use in the Simscape model.\n');
