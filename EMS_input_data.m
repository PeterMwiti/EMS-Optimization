function [Ppv, Pload, lambda, Pcap] = EMS_input_data(N, params)

% Assign defaults if inputs are omitted
    if nargin < 1 || isempty(N)
        N = 24; % Default to 24-hour horizon
    end
    if nargin < 2 || isempty(params)
        params = struct(); % Initialize empty struct to use internal defaults
    end

% ---- 1. Load Datasets from MAT-file ----
    if ~exist('pvLoadPriceData.mat', 'file')
        error('pvLoadPriceData.mat not found. Ensure the dataset file exists in the MATLAB path.');
    end
    data = load('pvLoadPriceData.mat');

    % Generate 1-based 24-hour cyclic indices for N simulation hours
    hod_idx = mod(0:N-1, 24) + 1;

    % ---- 2. Extract Solar Day Type Irradiance (W/m^2) ----
    day_type = 'clear'; % Default selection
    if isfield(params, 'pv_day_type') && ~isempty(params.pv_day_type)
        day_type = params.pv_day_type;
    end

    switch lower(day_type)
        case 'clear'
            gti_profile = data.clearDay;
        case {'cloudy', 'overcast'}
            gti_profile = data.cloudyDay;
        case {'partlycloudy', 'partly_cloudy', 'mixed'}
            gti_profile = data.partlyCloudyDay;
        otherwise
            warning('Unrecognized pv_day_type "%s". Defaulting to "clear".', day_type);
            gti_profile = data.clearDay;
    end

    % ---- 3. Convert Irradiance (W/m^2) to Active PV Power (kW) ----
    % Default to 1360.61 kWp if Ppv_rated is not provided in params struct
    Ppv_rated = 1360.61; 
    if isfield(params, 'Ppv_rated') && ~isempty(params.Ppv_rated)
        Ppv_rated = params.Ppv_rated;
    end

    % Formula: P_pv (kW) = Ppv_rated (kWp) * (GTI Irradiance / 1000 W/m^2)
    pv_kw = gti_profile .* (Ppv_rated / 1000);

    % Expand to N hours
    Ppv = pv_kw(hod_idx);

    % ---- 4. Real C&I Load Profile (kW) ----
    load_col = 2; % Default to load profile column 1
    if isfield(params, 'load_col') && ~isempty(params.load_col)
        load_col = params.load_col;
    end
    
    % Get constant baseload offset (default to 350 kW)
    Pload_base = 350;
    if isfield(params, 'Pload_base') && ~isempty(params.Pload_base)
        Pload_base = params.Pload_base;
    end
    
    if isfield(data, 'loadData_kW')
        load_base = data.loadData_kW(:, load_col);
    else
        load_base = data.loadData(:, load_col);
    end

    % Reconstruct total gross load: constant baseload + variable profile
    Pload = Pload_base + load_base(hod_idx); % Expands and shifts profile to (195.81 kW -> 695.05 kW)

    % ---- 5. Real Tariff Schedule (KES/kWh) ----
    lambda = data.costData(hod_idx); % Expand to N hours % costData = peakTOU + Levies (FCC, WERMA)

    % ---- 6. Soft Planned-Draw Ceiling Pcap(t) ----
    net_need = max(Pload - Ppv, 0);
    Pcap = 1.10 * net_need;
    
    peak_load = max(Pload);
    if isfield(params, 'Pload_peak') && ~isempty(params.Pload_peak)
        peak_load = params.Pload_peak;
    end
    Pcap = max(Pcap, 0.05 * peak_load);

    % Ensure outputs are strict N x 1 column vectors
    Ppv    = Ppv(:);
    Pload  = Pload(:);
    lambda = lambda(:);
    Pcap   = Pcap(:);

end