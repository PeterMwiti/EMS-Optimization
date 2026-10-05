function base = baseline_EMS_dispatch(Pload, lambda, CD, dt, PF, Papparent_meas)
% baseline_EMS_dispatch Grid-only, no-PV, no-BESS scenario.
% Works for either data source: uses the REAL measured apparent power
% as the ground-truth peak when it's available (THIWASCO), and falls
% back to a PF-based estimate when it isn't (MAT2, which has no
% measured kVA column at all).
%
% INPUTS
%   Pload          - actual load series (kW), column vector
%   lambda         - actual tariff series (KES/kWh), same length
%   CD             - demand charge rate (KES/kVA)
%   dt             - timestep (h), normally 1
%   PF             - scalar power factor, used ONLY as the fallback
%                    estimate when Papparent_meas is unavailable
%   Papparent_meas - real measured apparent power (kVA), same length
%                    as Pload, OR [] / omitted if not available for
%                    this dataset
%
% OUTPUT base (struct)
%   Pgrid      - grid draw each hour (kW) = Pload (no PV/BESS to offset it)
%   energy_cost- total TOU energy cost (KES)
%   demand_cost- demand charge on the single peak (KES)
%   total_cost - energy_cost + demand_cost (KES)
%   peak       - the peak APPARENT demand recorded or estimated (kVA)
%   peak_is_measured - true if peak came from real data, false if
%                      it's a PF-based estimate (useful for labeling
%                      plots/tables so the two aren't mixed up)

    if nargin < 6
        Papparent_meas = [];
    end

    base.Pgrid = Pload(:);
    base.energy_cost = sum(lambda(:) .* base.Pgrid .* dt);

    if ~isempty(Papparent_meas)
        base.peak = max(Papparent_meas(:));   % real measured kVA peak
        base.peak_is_measured = true;
    else
        base.peak = max(base.Pgrid) / PF;      % PF-based estimate
        base.peak_is_measured = false;
    end

    base.demand_cost = CD * base.peak;
    base.total_cost  = base.energy_cost + base.demand_cost;
end
