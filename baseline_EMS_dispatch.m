
function base = baseline_EMS_dispatch(Pload, lambda, CD, dt, pf)
% BASELINE_DISPATCH  Grid-only, no-PV, no-BESS scenario.
%
% This is the simplest possible comparator: every kW of load is served
% directly from the grid, billed under the same TOU + demand-charge
% structure as the proposed system, so the ONLY difference between this
% and the EMS result is the presence/absence of the PV-BESS system and
% its optimized control.
%
% INPUTS
%   Pload  - actual load series (kW), column vector
%   lambda - actual TOU tariff series (KES/kWh), same length as Pload
%   CD     - Maximum demand charge rate (KES/kVA)
%   dt     - timestep (h), normally 1
%   pf     - site power factor (e.g., 0.920) [OPTIONAL, default = 0.920]
%
% OUTPUT base (struct)
%   Pgrid       - grid draw each hour (kW) = Pload (no PV/BESS to offset it)
%   energy_cost - total peak TOU energy cost (KES)
%   demand_cost - maximum demand charge on the single peak (KES)
%   total_cost  - energy_cost + demand_cost (KES)
%   peak        - the peak demand recorded (kW)
%   peak_kva    - the peak demand recorded in kVA

    % Fallback to 0.920 power factor if not passed explicitly
    if nargin < 5 || isempty(pf)
        pf = 0.920;
    end

    base.Pgrid       = Pload(:);
    base.energy_cost = sum(lambda(:) .* base.Pgrid .* dt);
    
    % Peak active power (kW) and apparent power (kVA)
    base.peak        = max(base.Pgrid);          % kW
    base.peak_kva    = base.peak / pf;           % kVA = kW / PF
    
    % Demand charge calculation using kVA
    base.demand_cost = CD * base.peak_kva;       % (KES/kVA) * kVA
    base.total_cost  = base.energy_cost + base.demand_cost;
end