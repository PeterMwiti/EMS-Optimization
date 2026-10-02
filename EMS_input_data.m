function [Ppv, Pload, lambda, Pcap] = EMS_input_data(N, params)
% GENERATE_INPUT_DATA  Build hourly PV, load, tariff and soft-cap arrays.
%
% THIS IS A PLACEHOLDER. Section 3.3 of the proposal calls for real
% Kenyan GHI data, a real/published C&I load curve, and the actual KPLC
% CI1-CI6 tariff schedule. Replace the body of this function with code
% that loads and resamples your real .csv/.xlsx data to hourly
% resolution, keeping the same four output arrays (Ppv, Pload, lambda,
% Pcap), and everything downstream (solver, RH driver (MPC), plotting) will
% keep working unchanged.
%
% INPUTS
%   N       - number of hourly points to generate (must cover
%             Nsim + T_horizon - 1 hours, see MILP_MPCwrapper.m)
%   params  - struct with at least: Ppv_rated (kWp), Pload_base (kW),
%             Pload_peak (kW)
%
% OUTPUTS (all column vectors, length N, hourly resolution)
%   Ppv    - forecast/actual PV power (kW)
%   Pload  - forecast/actual C&I load (kW)
%   lambda - grid energy tariff (KES/kWh), TOU + FCC/FERFA surcharge
%   Pcap   - soft "planned ceiling" on grid draw (kW), used only in the
%            penalty term (Step 7 of the design) -- NOT a hard limitz

% hours = (0:N-1)'; %(0:N-1) ^ T(Transpose) - How far into the run are we?
% hod   = mod(hours, 24);          % hour-of-day, 0-23 - What time of day is it?
% day   = floor(hours/24); %which day are we on? 
% 
% % ---- 1. PV generation: raised-cosine daylight curve, 6:00-18:00 ----
% 
% daylight = (hod >= 6) & (hod <= 18); % daylight is between 6 am and 6 pm
% theta = pi * (hod - 6) / 12;                 % 0 at 6:00, 0.5 * pi at noon, pi at 18:00 - converts clock time into an angle
% pv_shape = zeros(N,1); % Initializes the pv_shape matrix with zeros
% pv_shape(daylight) = sin(theta(daylight)).^1.2; %overwrites the daylight slots with the actual values leaving the night time slots at zero
% cloud_factor = 0.75 + 0.25*sin(2*pi*day/7 + 1); % mild day-to-day cloud variability - adds day to day cloudiness
% cloud_factor = max(cloud_factor, 0.4); % creates a cloudiness floor that ensures the cloud variability doesn't fall below it
% Ppv = params.Ppv_rated .* pv_shape .* cloud_factor; %converts the pv_shape and cloud factor into actual KW using the solar systems rated capacity.
% 
% % ---- 2. C&I load: baseline + two-peak industrial shape + noise ----
% 
% morning = exp(-((hod-9).^2)/(2*2.2^2)); %creates a smooth curve with the peak occuring at 9 am; 2.2 determines how wide the hump is
% afternoon = exp(-((hod-15).^2)/(2*2.5^2));  %creates a smooth curve with the peak occuring at 3 pm; 2.5 determines how wide the hump is
% shape = 0.35 + 0.65*(morning + afternoon)/max(morning+afternoon); % provides the daily silhouette with a floor of 0.35 and a ceiling of 1.0 
% noise = 1 + 0.04*sin(3*hours/7); % creates a slight hourly variation in load
% Pload = params.Pload_base + (params.Pload_peak - params.Pload_base) .* shape .* noise; %converts to actual KW
% 
% % ---- 3. KPLC-style TOU tariff (KES/kWh), from Section 3.3 ranges ----
% % Peak window: 06:00-22:00 -> 13.44 KES/kWh for CI1
% % Off-peak:    22:00-06:00 -> 6.72 KES/kWh for CI1
% % FCC/FERFA surcharge added as a flat adder (replace with real value)
% peak_hours = (hod >= 6) & (hod < 22);
% lambda = zeros(N,1);
% lambda(peak_hours)  = 13.44;   % KES/kWh, CI1
% %lambda(~peak_hours) = 6.72;    % KES/kWh, CI1
% %fcc_ferfa_adder = 8.65;        % KES/kWh,
% lambda = 22.09; %+ fcc_ferfa_adder
% 
% % ---- 4. Soft planned-draw ceiling P_cap(t) ----
% % Placeholder rule: 110% of the forecast net grid need (load - PV),
% % floored at a small minimum so Pcap is never zero/negative.
% 
% net_need = max(Pload - Ppv, 0); %before the battery gets involved, how much power would need to come from the grid at this hour
% Pcap = 1.10 * net_need; %planned ceiling with a 10% extra breathing room
% Pcap = max(Pcap, 0.05*params.Pload_peak); % makes sure the ceiling never goes to near zero - can occur when net_need is an almost zero value 
% end



% Assign defaults if inputs are omitted
    if nargin < 1 || isempty(N)
        N = 24; % Default to 24-hour horizon
    end
    if nargin < 2 || isempty(params)
        params = struct(); % Initialize empty struct to use internal defaults
    end

% ---- 1. Load Datasets from MAT-file ----
    if ~exist('pvLoadPriceData_New.mat', 'file')
        error('pvLoadPriceData_New.mat not found. Ensure the dataset file exists in the MATLAB path.');
    end
    data = load('pvLoadPriceData_New.mat');

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
    load_col = 1; % Default to load profile column 1
    if isfield(params, 'load_col') && ~isempty(params.load_col)
        load_col = params.load_col;
    end

    if isfield(data, 'loadData_kW')
        load_base = data.loadData_kW(:, load_col);
    else
        load_base = data.loadData(:, load_col);
    end

    Pload = load_base(hod_idx); % Expand to N hours

    % ---- 5. Real Tariff Schedule (KES/kWh) ----
    lambda = data.costData(hod_idx); % Expand to N hours

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