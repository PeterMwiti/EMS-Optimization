% function sol = EMS_MILP_Design(Ppv, Pload, lambda, Pcap, E0, Ppeak_floor, p)
% % SOLVE_EMS_MILP  Build and solve one receding-horizon MILP dispatch.
% %
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
%     % ---------------- variable layout ----------------
%     % x = [Pgrid(1:T); Pch(1:T); Pdis(1:T); E(1:T); ych(1:T); ydis(1:T); s(1:T); Ppeak]
%     nvar   = 7*T + 1; % (total length of that column vector)- 169
%     idxPg  = 1:T; %(grid power value per hour)- 1 to 24
%     idxPch = T+1:2*T; %(charging power value per hour) - 25 to 48
%     idxPd  = 2*T+1:3*T; %(discharging power value per hour)- 49 to 72
%     idxE   = 3*T+1:4*T; %(battery energy level per hour)- 73 to 96
%     idxYch = 4*T+1:5*T; %("is it charging?" binary flag per hour)- 97 to 120
%     idxYd  = 5*T+1:6*T; %("is it discharging?" binary flag per hour)- 121 to 144
%     idxS   = 6*T+1:7*T; %(penalty/overshoot slack value per hour)- 145 to 168
%     idxPk  = 7*T+1; %(peak demand value)- 169
% 
%     % ---------------- objective (Step 6) ----------------
%     f = zeros(nvar,1);
%     f(idxPg) = lambda(:) * dt;   % energy cost term
%     f(idxS)  = p.rho;            % linear penalty term
%     f(idxPk) = p.CD;             % maximum demand-charge term
% 
%     % ---------------- equality constraints ----------------
%     % Row block 1 (t=1..T): power balance  (Step 2)
%     %   Pgrid(t) + Pdis(t) - Pch(t) = Pload(t) - Ppv(t)
%     % Row block 2 (t=1..T): battery energy dynamics (Step 3)
%     %   E(t) - eta_ch*dt*Pch(t) + (dt/eta_dis)*Pdis(t) - E(t-1) = 0
%     %   (E(0) replaced by the known constant E0 for t=1)
%     nEq  = 2*T; %(total number of equality equations; 2*24=48 - The '2' is as a result of the power balance and battery energy dynamics equations)
%     Aeq  = zeros(nEq, nvar); % (creates a spreadsheet of 48 by 169 filled with zeros)
%     beq  = zeros(nEq, 1); % (creates a 48 by 1 spreadsheet filled with zeros to be later replaced with the solved values for the power balance and battery energy dynamics equations)
% 
%     for t = 1:T
%         r = t;
%         Aeq(r, idxPg(t)) =  1;
%         Aeq(r, idxPd(t)) =  1;
%         Aeq(r, idxPch(t)) = -1;
%         beq(r) = Pload(t) - Ppv(t);
%     end
%     for t = 1:T
%         r = T + t;
%         Aeq(r, idxE(t))   =  1;
%         Aeq(r, idxPch(t)) = -p.eta_ch * dt;
%         Aeq(r, idxPd(t))  =  dt / p.eta_dis;
%         if t == 1
%             beq(r) = E0;
%         else
%             Aeq(r, idxE(t-1)) = -1;
%             beq(r) = 0;
%         end
%     end
% 
%     % ---------------- inequality constraints ----------------
%     % (3) Pch(t)  - Pch_max*ych(t)  <= 0        mutual-exclusivity cap (Step 4)
%     % (4) Pdis(t) - Pdis_max*ydis(t)<= 0        mutual-exclusivity cap (Step 4)
%     % (5) ych(t) + ydis(t)          <= 1        mutual exclusivity (Step 4)
%     % (6) Pgrid(t) - Ppeak          <= 0        peak-tracking (Step 5)
%     % (7) Pgrid(t) - s(t)           <= Pcap(t)  soft cap / penalty (Step 7)
%     % (8) -Ppeak                    <= -Ppeak_floor   carry-forward peak floor
%     nIneq = 4*T + T + 1;   % (3)+(4)+(5)+(6) each T rows, (7) T rows, (8) 1 row; total = 121 rows
%     A = zeros(nIneq, nvar); %( 121 by 169 spreadsheet filled with zeros)
%     b = zeros(nIneq, 1); % ( 121 by 1 spreadsheet that provides the solved values for the inequality constraints equations)
%     r = 0;
% 
%     for t = 1:T   % (3)
%         r = r+1;
%         A(r, idxPch(t)) = 1;
%         A(r, idxYch(t)) = -p.Pch_max;
%         b(r) = 0;
%     end
%     for t = 1:T   % (4)
%         r = r+1;
%         A(r, idxPd(t)) = 1;
%         A(r, idxYd(t)) = -p.Pdis_max;
%         b(r) = 0;
%     end
%     for t = 1:T   % (5)
%         r = r+1;
%         A(r, idxYch(t)) = 1;
%         A(r, idxYd(t))  = 1;
%         b(r) = 1;
%     end
%     for t = 1:T   % (6)
%         r = r+1;
%         A(r, idxPg(t)) = 1;
%         A(r, idxPk)    = -1;
%         b(r) = 0;
%     end
%     for t = 1:T   % (7)
%         r = r+1;
%         A(r, idxPg(t)) = 1;
%         A(r, idxS(t))  = -1;
%         b(r) = Pcap(t);
%     end
%     r = r+1;      % (8) peak floor: Ppeak >= Ppeak_floor
%     A(r, idxPk) = -1;
%     b(r) = -Ppeak_floor;
% 
%     % ---------------- bounds ----------------
%     lb = zeros(nvar,1); %(lower bound - every variable floor is zero initially (169 by 1))
%     ub = inf(nvar,1); %(upper bound - every variable ceiling is infinity initially (169 by 1))
%     ub(idxPg)  = p.Pgrid_max;
%     ub(idxPch) = p.Pch_max;
%     ub(idxPd)  = p.Pdis_max;
%     lb(idxE)   = p.Emin;
%     ub(idxE)   = p.Emax;
%     ub(idxYch) = 1;  ub(idxYd) = 1;   % binaries in {0,1}
%     lb(idxPk)  = Ppeak_floor;         % Ppeak can't record less than the past
% 
%     intcon = [idxYch, idxYd]; %(Forces the variables to be whole numbers)
% 
%     % ---------------- solve ----------------
%     opts = optimoptions('intlinprog', 'Display', 'off');
%     [x, fval, exitflag] = intlinprog(f, intcon, A, b, Aeq, beq, lb, ub, opts); %solver ( x- actual answer,
%                                                                                         %fval-resulting total cost at that answer
%                                                                                         %exitflag-status code saying how it went.i.e., found the best answer, ran out of time, found nothing that satisfies the constraint)
% 
%     sol = struct(); %creates a container for packaging the results instead of the 169 column list
% 
%     if isempty(x) % handles the situation where no answer exists (NaN)
%         warning('solve_EMS_MILP:infeasible', ...
%             'MILP returned no solution (exitflag=%d). Check Pcap/rho/limits.', exitflag);
%         sol.exitflag = exitflag; sol.fval = NaN;
%         sol.Pgrid = nan(T,1); sol.Pch = nan(T,1); sol.Pdis = nan(T,1);
%         sol.E = nan(T,1); sol.s = nan(T,1); sol.Ppeak = NaN;
%         return
%     end
%     % outputs the real answers obtained
%     sol.Pgrid = x(idxPg);
%     sol.Pch   = x(idxPch);
%     sol.Pdis  = x(idxPd);
%     sol.E     = x(idxE);
%     sol.s     = x(idxS);
%     sol.Ppeak = x(idxPk);
%     sol.exitflag = exitflag;
%     sol.fval = fval;
% end



function sol = EMS_MILP_Design(Ppv, Pload, lambda, Pcap, E0, Ppeak_floor, p)
% SOLVE_EMS_MILP  Build and solve one receding-horizon MILP dispatch.

    T   = length(Ppv);
    dt  = p.dt;

    % ---------------- variable layout ----------------
    % x = [Pgrid(1:T); Pch(1:T); Pdis(1:T); E(1:T); ych(1:T); ydis(1:T); s(1:T); Pcurt(1:T); Ppeak]
    nvar     = 8*T + 1; 
    idxPg    = 1:T;
    idxPch   = T+1:2*T;
    idxPd    = 2*T+1:3*T;
    idxE     = 3*T+1:4*T;
    idxYch   = 4*T+1:5*T;
    idxYd    = 5*T+1:6*T;
    idxS     = 6*T+1:7*T;
    idxPcurt = 7*T+1:8*T;
    idxPk    = 8*T+1;

    % ---------------- objective (Step 6) ----------------
    f = zeros(nvar,1);
    f(idxPg) = lambda(:) * dt;   % energy cost term
    f(idxS)  = p.rho;            % linear penalty term
    f(idxPk) = p.CD;             % maximum demand-charge term

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

    for t = 1:T   % (5) ych(t) + ydis(t) <= 1
        r = r+1;
        A(r, idxYch(t)) = 1;
        A(r, idxYd(t))  = 1;
        b(r) = 1;
    end

    for t = 1:T   % (6) Pgrid(t) - Ppeak <= 0
        r = r+1;
        A(r, idxPg(t)) = 1;
        A(r, idxPk)    = -1;
        b(r) = 0;
    end

    for t = 1:T   % (7) Pgrid(t) - s(t) <= Pcap(t)
        r = r+1;
        A(r, idxPg(t)) = 1;
        A(r, idxS(t))  = -1;
        b(r) = Pcap(t);
    end

    r = r+1;      % (8) peak floor: Ppeak >= Ppeak_floor
    A(r, idxPk) = -1;
    b(r) = -Ppeak_floor;

    % ---------------- bounds ----------------
    lb = zeros(nvar,1); 
    ub = inf(nvar,1); 

    ub(idxPg)    = p.Pgrid_max;
    ub(idxPch)   = p.Pch_max;
    ub(idxPd)    = p.Pdis_max;
    lb(idxE)     = p.Emin;
    ub(idxE)     = p.Emax;
    ub(idxYch)   = 1;  
    ub(idxYd)    = 1;   
    lb(idxPk)    = Ppeak_floor; 
    lb(idxPcurt) = 0;
    ub(idxPcurt) = inf;

    intcon = [idxYch, idxYd]; 

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
        return
    end

    % outputs the real answers obtained
    sol.Pgrid = x(idxPg);
    sol.Pch   = x(idxPch);
    sol.Pdis  = x(idxPd);
    sol.E     = x(idxE);
    sol.s     = x(idxS);
    sol.Pcurt = x(idxPcurt);
    sol.Ppeak = x(idxPk);
    sol.exitflag = exitflag;
    sol.fval = fval;
end