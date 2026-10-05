function [Pload, PF, Papparent_meas] = real_CI_load_profile(N, csv_path)
% REAL_CI_LOAD_PROFILE  Load the real, meter-measured C&I load profile.
%
% Reads the pre-cleaned CSV built from the Hexing meter export
% ("Load Prof" report, 19-25 April 2026). See the header comments in
% the extraction script this CSV was produced from for how the
% ambiguous 12-hour meter timestamps were reconstructed into true
% 24-hour time, and how the small number of missing readings (5 of
% 168 hours) were linearly interpolated -- both are flagged in the
% CSV's own 'interpolated' column for transparency.
%
% CHANGE LOG (this version):
%   - PF is now a FLAT constant, 0.92, for every hour -- no longer
%     read from the CSV's own PowerFactor column (which varied
%     0.847-0.941 across the real week). This is a deliberate
%     simplification: using one fixed, known-in-advance PF removes
%     the real-vs-forecast PF mismatch as a source of variability,
%     at the cost of the small amount of realism the per-hour PF
%     captured. Papparent_meas (the baseline's ground-truth kVA) is
%     UNCHANGED -- it's read directly from its own measured column,
%     never derived via PF, so it isn't affected by this change.
%
% INPUTS
%   N        - number of hourly points requested. Must equal the
%              number of rows in the CSV (168, i.e. exactly the
%              7-day/168-hour real window) -- this function does NOT
%              extrapolate beyond the real measured week, since
%              fabricating load data beyond what was actually metered
%              would undermine the whole point of using real data.
%   csv_path - path to KPLC_real_load_Apr19to25_2026.csv
%
% OUTPUTS (column vectors, length N)
%   Pload          - real active power demand (kW) -- this is what
%                    feeds the MILP's power balance (Step 2)
%   PF             - FLAT 0.92 for every hour -- used to convert the
%                    MILP's optimized Pgrid(t) [kW] into an equivalent
%                    apparent power [kVA] for demand-charge tracking
%                    (Step 5). Kept as a length-N vector (rather than
%                    a bare scalar) purely so EMS_MILP_Design.m and
%                    MILP_MPCwrapper.m need no changes to their own
%                    interfaces -- every element is just 0.92.
%   Papparent_meas - real MEASURED apparent power (kVA) at the meter --
%                    the historical, no-BESS ground truth, used by
%                    baseline_EMS_dispatch.m. Independent of PF.

    T = readtable(csv_path);
    if height(T) ~= N
        error(['real_CI_load_profile: CSV has %d rows but %d were ' ...
               'requested. This dataset covers exactly one real ' ...
               '168-hour week -- set Nsim so that Nsim+T_horizon-1 ' ...
               '== 168 (e.g. Nsim=145, T_horizon=24).'], height(T), N);
    end
    Pload          = T.Pload_kW;
    Papparent_meas = T.Papparent_kVA;      % real measured, unaffected by the PF change below
    PF             = 0.92 * ones(N,1);     % flat, constant power factor (was: T.PowerFactor)

    n_interp = sum(T.interpolated);
    if n_interp > 0
        fprintf(['real_CI_load_profile: %d of %d hours were linearly ' ...
                 'interpolated (missing in the source meter export).\n'], ...
                 n_interp, N);
    end
end