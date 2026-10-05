function sol = EMS_MILP_Design(Ppv, Pload, lambda, Pcap, E0, Ppeak_floor, p)
% EMS_MILP_Design  Build and solve one receding-horizon MILP
% dispatch. Works identically regardless of which dataset produced
% Ppv/Pload/lambda/Pcap -- see EMS_input_data_unified.m.
%
% Steps 1-7 of the EMS design:
%   Step 1 - decision variables, every hour of the horizon
%   Step 2 - power balance (equality), WITH CURTAILMENT
%   Step 3 - battery energy dynamics + SOC bounds
%   Step 4 - binary mutual-exclusivity of charge/discharge
%   Step 5 - peak-demand epigraph trick, IN kVA
%   Step 6 - objective function (energy cost + demand charge + penalty)
%   Step 7 - linear penalty / soft cap on planned grid draw
%
% UNIFICATION CHANGE: PF is no longer a per-hour input vector (as it
% was for THIWASCO) or something re-derived only at reporting time (as
% it was for MAT2/test2). It is now a single scalar, p.PF, used
% consistently both inside this solver's peak-tracking constraint and
% by the wrapper's own kVA reporting -- one definition, one place,
% works the same regardless of data source.
%
% (Off-peak battery lockout, solar curtailment, and the Ppeak_floor
% carry-forward mechanism are unchanged from the THIWASCO version of
% this file -- see EMS_MILP_Design.m's own change log for why each
% exists.)
%
% INPUTS (all forecasts over the horizon, length T = length(Ppv))
%   Ppv, Pload, lambda, Pcap : column vectors, length T
%   E0          : battery energy (kWh) at the START of the horizon
%                 (this must be the ACTUAL measured SOC, not a forecast)
%   Ppeak_floor : highest grid demand already recorded this billing
%                 period, in kVA (0 at the start of a new period)
%   p           : parameter struct -- dt, eta_ch, eta_dis, Emin, Emax,
%                 Pch_max, Pdis_max, Pgrid_max, CD, rho, PF
%
% OUTPUT sol (struct)
%   Pgrid, Pch, Pdis, E, s, Pcurt : column vectors, length T
%   ych, ydis : binary charge/discharge flags, column vectors, length T
%   Ppeak    : scalar, the horizon's tracked peak APPARENT demand (kVA)
%   exitflag, fval : intlinprog diagnostics

    T   = length(Ppv);
    dt  = p.dt;

    % ---------------- variable layout ----------------
    % x = [Pgrid(1:T); Pch(1:T); Pdis(1:T); E(1:T); ych(1:T); ydis(1:T); s(1:T); Pcurt(1:T); Ppeak]
    nvar    = 8*T + 1; %193
    idxPg   = 1:T; % 1-24
    idxPch  = T+1:2*T; % 25-48
    idxPd   = 2*T+1:3*T; % 49-72
    idxE    = 3*T+1:4*T; % 73-96
    idxYch  = 4*T+1:5*T; % 97-120
    idxYd   = 5*T+1:6*T; % 121-144
    idxS    = 6*T+1:7*T; % 145-168
    idxCurt = 7*T+1:8*T; % 169-192
    idxPk   = 8*T+1; %193

    % ---------------- objective (Step 6) ----------------
    f = zeros(nvar,1);
    f(idxPg) = lambda(:) * dt;
    f(idxS)  = p.rho;
    f(idxPk) = p.CD;
    % f(idxCurt) intentionally left at 0

    % ---------------- equality constraints ----------------
    nEq  = 2*T;
    Aeq  = zeros(nEq, nvar);
    beq  = zeros(nEq, 1);
    for t = 1:T
        r = t;
        Aeq(r, idxPg(t))   =  1;
        Aeq(r, idxPd(t))   =  1;
        Aeq(r, idxPch(t))  = -1;
        Aeq(r, idxCurt(t)) = -1;
        beq(r) = Pload(t) - Ppv(t);
    end
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
    % (3) Pch(t)  - Pch_max*ych(t)  <= 0
    % (4) Pdis(t) - Pdis_max*ydis(t)<= 0
    % (5) ych(t) + ydis(t)          <= 1
    % (6) Pgrid(t)/PF - Ppeak       <= 0   (PF is now a SCALAR -- same
    %                                       coefficient every hour)
    % (7) Pgrid(t) - s(t)           <= Pcap(t)
    % (8) -Ppeak                    <= -Ppeak_floor
    nIneq = 4*T + T + 1;
    A = zeros(nIneq, nvar);
    b = zeros(nIneq, 1);
    r = 0;

    for t = 1:T   % (3)
        r = r+1;
        A(r, idxPch(t)) = 1;
        A(r, idxYch(t)) = -p.Pch_max;
        b(r) = 0;
    end
    for t = 1:T   % (4)
        r = r+1;
        A(r, idxPd(t)) = 1;
        A(r, idxYd(t)) = -p.Pdis_max;
        b(r) = 0;
    end
    for t = 1:T   % (5)
        r = r+1;
        A(r, idxYch(t)) = 1;
        A(r, idxYd(t))  = 1;
        b(r) = 1;
    end
    for t = 1:T   % (6)
        r = r+1;
        A(r, idxPg(t)) = 1 / p.PF;
        A(r, idxPk)    = -1;
        b(r) = 0;
    end
    for t = 1:T   % (7)
        r = r+1;
        A(r, idxPg(t)) = 1;
        A(r, idxS(t))  = -1;
        b(r) = Pcap(t);
    end
    r = r+1;      % (8)
    A(r, idxPk) = -1;
    b(r) = -Ppeak_floor;

    % ---------------- bounds ----------------
    lb = zeros(nvar,1);
    ub = inf(nvar,1);
    ub(idxPg)   = p.Pgrid_max;
    ub(idxPch)  = p.Pch_max;
    ub(idxPd)   = p.Pdis_max;
    lb(idxE)    = p.Emin;
    ub(idxE)    = p.Emax;
    ub(idxYch)  = 1;  ub(idxYd) = 1;
    ub(idxCurt) = Ppv(:);        % can't curtail more solar than was ever generated
    lb(idxPk)   = Ppeak_floor;

    % Force the battery idle during off-peak hours. Portable heuristic:
    % if a dataset's lambda never dips below 10 (e.g. MAT2's combined,
    % already-surcharged tariff), this simply never triggers -- that's
    % correct behaviour, not a bug, since such a dataset has no
    % off-peak discount to exploit in the first place.
    is_offpeak = lambda(:) < 10;
    ub(idxYch(is_offpeak)) = 0;
    ub(idxYd(is_offpeak))  = 0;

    intcon = [idxYch, idxYd];

    % ---------------- solve ----------------
    opts = optimoptions('intlinprog', 'Display', 'off');
    [x, fval, exitflag] = intlinprog(f, intcon, A, b, Aeq, beq, lb, ub, opts);

    sol = struct();
    if isempty(x)
        warning('EMS_MILP_Design_unified:infeasible', ...
            'MILP returned no solution (exitflag=%d). Check Pcap/rho/limits.', exitflag);
        sol.exitflag = exitflag; sol.fval = NaN;
        sol.Pgrid = nan(T,1); sol.Pch = nan(T,1); sol.Pdis = nan(T,1);
        sol.E = nan(T,1); sol.s = nan(T,1); sol.Ppeak = NaN;
        sol.ych = nan(T,1); sol.ydis = nan(T,1); sol.Pcurt = nan(T,1);
        return
    end
    sol.Pgrid = x(idxPg);
    sol.Pch   = x(idxPch);
    sol.Pdis  = x(idxPd);
    sol.E     = x(idxE);
    sol.s     = x(idxS);
    sol.Pcurt = x(idxCurt);
    sol.ych   = round(x(idxYch));
    sol.ydis  = round(x(idxYd));
    sol.Ppeak = x(idxPk);
    sol.exitflag = exitflag;
    sol.fval = fval;

    assert(all(sol.ych + sol.ydis <= 1), ...
        'EMS_MILP_Design_unified:SCD_violation', ...
        'Simultaneous charge+discharge detected -- Step 4 constraint was violated.');
end
