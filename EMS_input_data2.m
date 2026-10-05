function [Ppv, Pload, lambda, Pcap, Papparent_meas] = EMS_input_data2(N, params)
% EMS_Input_Data2  Build hourly PV, load, tariff and soft-cap
% arrays for EITHER the THIWASCO real-load dataset OR the
% pvLoadPriceData MAT-file dataset (MILP_2), selected via
% params.data_source, so the SAME downstream solver/wrapper/baseline
% code works unchanged with either data source.
%
% WHAT'S UPDATED:
%   - PV is the REAL solar irradiance profile from
%     pvLoadPriceData.mat (clearDay/cloudyDay/partlyCloudyDay, a
%     repeating 24-hour GTI template, W/m^2), converted to kW via
%     Ppv = Ppv_rated * (GTI/1000). This REPLACES the old synthetic
%     raised-cosine PV curve previously used for THIWASCO -- that was
%     always a documented placeholder, and real irradiance data is
%     now available, so there is no reason to keep both approaches.
%   - Load is dataset-specific, since the two sources store it
%     differently:
%       * 'THIWASCO': THIWASCO_real_load_Apr19to25_2026.csv's
%         Pload_kW column is already the site's TOTAL load (constant
%         + variable combined) -- used as-is, exactly 168 real hours,
%         no extrapolation beyond the real measured week.
%       * 'MAT2'    : pvLoadPriceData.mat's loadData_kW is a
%         VARIABLE-only profile; a constant baseload
%         (params.Pload_base, default 350 kW) is added on top to get
%         the total, and the 24-hour template repeats cyclically for
%         any N.
%   - Tariff (lambda) is also dataset-specific:
%       * 'THIWASCO': built from the CI1 TOU schedule (13.44 peak /
%         6.72 off-peak + a placeholder FCC/FERFA adder).
%       * 'MAT2'    : read directly from pvLoadPriceData.mat's
%         costData, which ALREADY sums the peak CI1 rate (13.44) and
%         the additional charges (FCC, FERMA, etc.) -- used as one
%         combined value per hour, with no further peak/off-peak
%         split applied on top of it.
%   - Power factor (PF) is NO LONGER returned from here. Both
%     datasets use one constant, known power factor in practice, so
%     it now lives as a single scalar, p.PF, in the parameter struct
%     used by the solver and wrapper (see
%     EMS_MILP_Design.m / MILP_MPCwrapper.m). This
%     removes the earlier vector-vs-scalar inconsistency between the
%     two original pipelines.
%   - Papparent_meas (real measured apparent power, kVA) is only
%     available for 'THIWASCO' (the CSV has a genuine measured kVA
%     column); for 'MAT2' there is no such measurement, so this
%     output is returned as [] (empty).
%     baseline_EMS_dispatch.m accepts that and falls back to
%     a PF-based estimate when it's empty, so the SAME baseline
%     function works for both.
%
% INPUTS
%   N      - number of hourly points to generate.
%            'THIWASCO': MUST equal 168 (the full real week) -- this
%            dataset is not cyclic/extendable; see real_CI_load_profile.m.
%            'MAT2': any N -- the 24-hour templates repeat cyclically.
%   params - struct, fields used here:
%     .data_source  - 'THIWASCO' (default) or 'MAT2'
%     .Ppv_rated    - kWp, PV array rating (applies to either source;
%                     default 1360.61 if omitted)
%     .pv_day_type  - 'clear' (default) | 'cloudy'/'overcast' |
%                     'partlycloudy'/'mixed' -- which GTI template to use
%     .Pload_base   - kW, constant baseload. For 'THIWASCO' this is
%                     used ONLY for the Pcap floor fallback (the CSV
%                     load is already total, so no offset is added to
%                     it -- see EMS_input_data2.m's change log for why
%                     a naive add would double-count). For 'MAT2' it
%                     is genuinely added on top of the variable profile.
%     .Pload_peak   - kW, used for the Pcap floor; defaults to the
%                     loaded data's own max if omitted.
%     .load_col     - 'MAT2' only: which column of loadData_kW to use
%                     (default 1).
%
% OUTPUTS (all column vectors, length N, hourly resolution, except
% Papparent_meas which may be [])
%   Ppv            - REAL irradiance-derived PV power (kW)
%   Pload          - total site load (kW)
%   lambda         - grid energy tariff (KES/kWh)
%   Pcap           - soft "planned ceiling" on grid draw (kW) -- Step 7
%   Papparent_meas - real measured apparent power (kVA), THIWASCO only;
%                    [] for MAT2

    if nargin < 2 || isempty(params)
        params = struct();
    end
    data_source = 'THIWASCO';
    if isfield(params, 'data_source') && ~isempty(params.data_source)
        data_source = params.data_source;
    end

    hours   = (0:N-1)';
    hod     = mod(hours, 24);
    hod_idx = hod + 1;   % GTI/load/cost templates are 24-long

    % ---- 1. PV generation: ALWAYS real irradiance, from pvLoadPriceData.mat ----
    mat_path = fullfile(fileparts(mfilename('fullpath')), 'pvLoadPriceData.mat');
    if ~exist(mat_path, 'file')
        error('EMS_input_data2:missing_mat', ...
            ['pvLoadPriceData.mat not found next to this file -- required ' ...
             'for PV generation under both data sources now.']);
    end
    pvdata = load(mat_path);

    day_type = 'clear';
    if isfield(params, 'pv_day_type') && ~isempty(params.pv_day_type)
        day_type = params.pv_day_type;
    end
    switch lower(day_type)
        case 'clear'
            gti_profile = pvdata.clearDay;
        case {'cloudy', 'overcast'}
            gti_profile = pvdata.cloudyDay;
        case {'partlycloudy', 'partly_cloudy', 'mixed'}
            gti_profile = pvdata.partlyCloudyDay;
        otherwise
            warning('EMS_input_data_unified:bad_day_type', ...
                'Unrecognized pv_day_type "%s". Defaulting to "clear".', day_type);
            gti_profile = pvdata.clearDay;
    end

    Ppv_rated = 1360.61;   % kWp, default; override via params.Ppv_rated
    if isfield(params, 'Ppv_rated') && ~isempty(params.Ppv_rated)
        Ppv_rated = params.Ppv_rated;
    end
    Ppv = Ppv_rated .* (gti_profile(hod_idx) ./ 1000);   % GTI [W/m^2] -> kW
    Ppv = Ppv(:);

    % ---- 2. Load + tariff: dataset-specific ----
    Pload_base = 350;
    if isfield(params, 'Pload_base') && ~isempty(params.Pload_base)
        Pload_base = params.Pload_base;
    end

    switch upper(data_source)
      case 'THIWASCO'
        csv_path = fullfile(fileparts(mfilename('fullpath')), ...
            'THIWASCO_real_load_Apr19to25_2026.csv');
        [Pload, ~, Papparent_meas] = real_CI_load_profile(N, csv_path);
        % Pload_kW is already the TOTAL load -- used as-is. The second
        % output (PF) is discarded here: PF now lives in p.PF instead.

        peak_hours = (hod >= 6) & (hod < 22);
        lambda = zeros(N,1);
        lambda(peak_hours)  = 13.44;   % KES/kWh, CI1 peak
        lambda(~peak_hours) = 6.72;    % KES/kWh, CI1 off-peak
        fcc_ferfa_adder = 0;           % KES/kWh, placeholder -- update when available
        lambda = lambda + fcc_ferfa_adder;

      case 'MAT2'
        load_col = 1;
        if isfield(params, 'load_col') && ~isempty(params.load_col)
            load_col = params.load_col;
        end
        if isfield(pvdata, 'loadData_kW')
            load_base = pvdata.loadData_kW(:, load_col);
        else
            load_base = pvdata.loadData(:, load_col);
        end
        Pload = Pload_base + load_base(hod_idx);   % variable profile + constant baseload
        Pload = Pload(:);
        Papparent_meas = [];   % not available for this source -- no measured kVA column

        % lambda ALREADY sums the peak CI1 rate (13.44) and the
        % additional charges (FCC, FERMA, etc.) -- used directly, no
        % further peak/off-peak split applied on top of it.
        lambda = pvdata.costData(hod_idx);
        lambda = lambda(:);

      otherwise
        error('EMS_input_data2:bad_source', ...
            'Unknown params.data_source "%s" -- use ''THIWASCO'' or ''MAT2''.', data_source);
    end

    % ---- 3. Soft planned-draw ceiling P_cap(t) (same rule, either source) ----
    net_need = max(Pload - Ppv, 0);
    Pcap = 1.10 * net_need;

    peak_load = max(Pload);   % robust default: the data's own max
    if isfield(params, 'Pload_peak') && ~isempty(params.Pload_peak)
        peak_load = params.Pload_peak;
    end
    Pcap = max(Pcap, 0.05*peak_load);
    Pcap = Pcap(:);
end