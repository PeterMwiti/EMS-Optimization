
% % This implements, exactly, Steps 1-7 of the EMS design:
% %   Step 1 - What the EMS has to actually decide every hour {24 hours}
% %   Step 2 - power balance (equality)
% %   Step 3 - battery energy dynamics + SOC bounds (equality + bounds)
% %   Step 4 - binary mutual-exclusivity of charge/discharge (MILP core)
% %   Step 5 - peak-demand epigraph trick (P_peak >= P_grid(t), all t)
% %   Step 6 - objective function (energy cost + demand charge + penalty)
% %   Step 7 - linear penalty / soft cap on planned grid draw
% %
% % INPUTS (all forecasts over the horizon, length T = length(Ppv))
% %   Ppv, Pload, lambda, Pcap : column vectors, length
% %   E0          : battery energy (kWh) at the START of the horizon
% %                 (this must be the ACTUAL measured SOC, not a forecast)
% %   Ppeak_floor : the highest grid demand (kW) already recorded so far
% %                 in the current billing period (0 at the start of a
% %                 new billing month) - ensures that the peak demand is
% %                 cumulative
% %   p           : parameter struct, see MILP_MPCwrapper.m for the
% %                 full field list (dt, eta_ch, eta_dis, Emin, Emax,
% %                 Pch_max, Pdis_max, Pgrid_max, CD, rho)
% %
% % OUTPUT sol (struct)
% %   Pgrid, Pch, Pdis, E, s : column vectors, length T
% %   Ppeak    : scalar, the horizon's tracked peak demand (kW)
% %   exitflag, fval : intlinprog diagnostics
% 
%     T   = length(Ppv);
%     dt  = p.dt;
% 
% 
%
% 
%  

function sol = EMS_MILP_Design(Ppv, Pload, lambda, Pcap, E0, Ppeak_floor, p)
% SOLVE_EMS_MILP  Build and solve one receding-horizon MILP dispatch.

    T   = length(Ppv);
    dt  = p.dt;

    % ---------------- variable layout ----------------
    % x = [Pgrid(1:T); Pch(1:T); Pdis(1:T); E(1:T); ych(1:T); ydis(1:T); s(1:T); Pcurt(1:T); Ppeak]
    nvar     = 8*T + 1;  % (length of x column vector)- 193
    idxPg    = 1:T;        % Grid power import [kW] (1 to T)
    idxPch   = T+1:2*T;    % Battery charging power [kW] (T+1 to 2T)
    idxPd    = 2*T+1:3*T;  % Battery discharging power [kW] (2T+1 to 3T)
    idxE     = 3*T+1:4*T;  % Battery energy level [kWh] (3T+1 to 4T)
    idxYch   = 4*T+1:5*T;  % Binary indicator: 1 if charging, 0 otherwise
    idxYd    = 5*T+1:6*T;  % Binary indicator: 1 if discharging, 0 otherwise
    idxS     = 6*T+1:7*T;  % Overshoot slack variable above Pcap [kW]
    idxPcurt = 7*T+1:8*T;  % Solar power curtailment [kW]
    idxPk    = 8*T+1;      % Maximum peak grid demand over the horizon [kW]

    % ---------------- objective (Step 6) ----------------
    f = zeros(nvar,1);
    f(idxPg) = lambda(:) * dt;   % energy cost term
    f(idxS)  = p.rho;            % linear penalty term
    f(idxPk) = p.CD;             % maximum demand-charge term
    CD_effective = p.CD / p.pf; % Effective demand charge per kW of peak demand (KES/kW)

    % Assign CD_effective to the peak variable index in objective vector f
    f(idxPk) = CD_effective;  % Ensure this uses CD_effective instead of p.CD

    % ---------------- equality constraints ----------------
    nEq  = 2*T; 
    Aeq  = zeros(nEq, nvar); 
    beq  = zeros(nEq, 1); 

    % Row block 1 (t=1..T): power balance with solar curtailment slack
    % Pgrid(t) + Pdis(t) - Pch(t) - Pcurt(t) = Pload(t) - Ppv(t)
    for t = 1:T
        r = t;
        Aeq(r, idxPg(t))    =  1;
        Aeq(r, idxPd(t))    =  1;
        Aeq(r, idxPch(t))   = -1;
        Aeq(r, idxPcurt(t)) = -1;
        beq(r) = Pload(t) - Ppv(t);
    end

    % Row block 2 (t=1..T): battery energy dynamics (Step 3)
    for t = 1:T
        r = T + t;
        Aeq(r, idxE(t))   =  1;
        Aeq(r, idxPch(t)) = -p.eta_ch * dt;
        Aeq(r, idxPd(t))  =  dt / p.eta_dis;
        if t == 1
            beq(r) = E0;
        else
            Aeq(r, idxE(t-1)) = -1;
            beq(r) = 0;
        end
    end

    % ---------------- inequality constraints ----------------
    nIneq = 4*T + T + 1;   
    A = zeros(nIneq, nvar); 
    b = zeros(nIneq, 1); 
    r = 0;

    for t = 1:T   % (3) Pch(t) - Pch_max*ych(t) <= 0
        r = r+1;
        A(r, idxPch(t)) = 1;
        A(r, idxYch(t)) = -p.Pch_max;
        b(r) = 0;
    end

    for t = 1:T   % (4) Pdis(t) - Pdis_max*ydis(t) <= 0
        r = r+1;
        A(r, idxPd(t)) = 1;
        A(r, idxYd(t))  = -p.Pdis_max;
        b(r) = 0;
    end

    for t = 1:T   % (5) ych(t) + ydis(t) <= 1 % Battery cannot charge and discharge at the same time
        r = r+1;
        A(r, idxYch(t)) = 1;
        A(r, idxYd(t))  = 1;
        b(r) = 1;
    end

    for t = 1:T   % (6) Pgrid(t) - Ppeak <= 0 % Peak tracking
        r = r+1;
        A(r, idxPg(t)) = 1;
        A(r, idxPk)    = -1;
        b(r) = 0;
    end

    for t = 1:T   % (7) Pgrid(t) - s(t) <= Pcap(t) % Soft cap penalty tracking
        r = r+1;
        A(r, idxPg(t)) = 1;
        A(r, idxS(t))  = -1;
        b(r) = Pcap(t);
    end

    r = r+1;      % (8) peak floor: Ppeak >= Ppeak_floor % Historical peak floor constraint
    A(r, idxPk) = -1;
    b(r) = -Ppeak_floor;

    % ---------------- bounds ----------------
    lb = zeros(nvar,1);    % Default lower bounds to 0
    ub = inf(nvar,1);      % Default upper bounds to infinity
    
    ub(idxPg)    = p.Pgrid_max;   % Transformer/grid connection limit
    ub(idxPch)   = p.Pch_max;     % Inverter max charge rating
    ub(idxPd)    = p.Pdis_max;    % Inverter max discharge rating
    lb(idxE)     = p.Emin;        % Battery minimum SOC (kWh)
    ub(idxE)     = p.Emax;        % Battery maximum capacity (kWh)
    ub(idxYch)   = 1;             % Binary flag upper bound
    ub(idxYd)    = 1;             % Binary flag upper bound
    lb(idxPk)    = Ppeak_floor;   % Peak demand cannot fall below historical monthly peak
    lb(idxPcurt) = 0;             % Curtailment cannot be negative
    ub(idxPcurt) = inf;           % Curtailment capped by available solar

    intcon = [idxYch, idxYd];     % Forces ych and ydis to be integers (0 or 1)

    % ---------------- solve ----------------
    opts = optimoptions('intlinprog', 'Display', 'off');
    [x, fval, exitflag] = intlinprog(f, intcon, A, b, Aeq, beq, lb, ub, opts); 

    sol = struct(); 
    if isempty(x) 
        warning('solve_EMS_MILP:infeasible', ...
            'MILP returned no solution (exitflag=%d). Check Pcap/rho/limits.', exitflag);
        sol.exitflag = exitflag; sol.fval = NaN;
        sol.Pgrid = nan(T,1); sol.Pch = nan(T,1); sol.Pdis = nan(T,1);
        sol.E = nan(T,1); sol.s = nan(T,1); sol.Pcurt = nan(T,1); sol.Ppeak = NaN;
        sol.ych = nan(T,1); sol.ydis = nan(T,1);
      
    % ---Performance Analysis Metrics to sol ---
    sol.wasted_solar_kWh = sum(sol.Pcurt) * dt;
    total_pv_available   = sum(Ppv);
    if total_pv_available > 0
        sol.solar_utilization_pct = 100 * (1 - (sum(sol.Pcurt) / total_pv_available));
    else
        sol.solar_utilization_pct = 100;
    end
 
        return
    end

    % outputs the real answers obtained
    sol.Pgrid = x(idxPg);
    sol.Pch   = x(idxPch);
    sol.Pdis  = x(idxPd);
    sol.E     = x(idxE);
    sol.s     = x(idxS);
    sol.ych   = round(x(idxYch));   % round(): intlinprog returns exact
    sol.ydis  = round(x(idxYd));    % 0/1 up to solver tolerance, not literal integers
    sol.Pcurt = x(idxPcurt);
    sol.Ppeak = x(idxPk);
    sol.exitflag = exitflag;
    sol.fval = fval;
end