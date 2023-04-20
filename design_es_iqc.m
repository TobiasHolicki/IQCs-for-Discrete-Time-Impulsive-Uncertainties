%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% File     : design_es_iqc.m                                              %
%                                                                         %
% Author   : Tobias Holicki                                               %
% Date     : 20.04.2023                                                   %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% This function implements the synthesis criteria for finding a 
% non-impulsive controller (an LTI estimator) such that the robust energy 
% gain of the closed-loop interconnection with the following impulsive 
% open-loop plant is minimized: 
%      x(t+1) =  A x(t) +  Bw w(t) +  Bd d(t)                  
%        z(t) = Cz x(t) + Dzw w(t) + Dzd d(t)
%        e(t) = Ce x(t) + Dew w(t) + Ded d(t) -u(t) 
%        y(t) = Cy x(t) + Dyw w(t) + Dyd d(t)
% where 
%     Del(z)(t) = z(t) if t = t(k) for some k and Del(z)(t) = 0  otherwise
% with (t(k))_k being a sequence that is merely known to satisfy a minimum 
% or range dwell-time condition.
%
% The synthesis criteria result from a generalization of the IQC theorem 
% from [1] with an IQC for the operator Del that relies on the discrete-
% time lifting procedure.
%
% [1] C. W. Scherer, J. Veenman, Stability analysis by dynamic dissipation
%     inequalities: On merging frequency-domain techniques with time domain
%     conditions, 2018
% ----- Input ---------------------------------------------------------- 
% sys               - System describing the linear part of the 
%                     open-loop plant
% inp               - Partition of the input signals [w, d, u]
% out               - Partition of the output signals [z, e, y]
% T                 - Assumed dwell-time conditions, i.e., 
%                           t(k+1) - t(k) - 1 \in [T(1), T(2)] 
%                     where T(2) can be Inf.    
% opt               - Struct
%    .psi           - The dynamic filter involved in the IQC for the 
%                     impulse operator.
%    .opt           - LMIlab solver options
%    .resolve_bound - Resolving the synthesis LMIs with a lower bound on
%                     the performance level is a simple way which can help
%                     to get better conditioned controller matrices.
%    .reconstruct   - Can be a time-saver if no controller reconstruction
%                     is desired.
% ----- Output ---------------------------------------------------------
% ga                - Numerically determined optimal upper bound on the 
%                     robust energy gain achieved by the designed
%                     estimator.
% E                 - The designed estimator
%
function [ga, E] = design_es_iqc(sys, inp, out, T, opt)
    % Some sanity checks
    arguments
        sys        {mustBeA(sys, "ss")}
        inp (1, 3) {mustBeInteger, mustBeNonnegative}
        out (1, 3) {mustBeInteger, mustBeNonnegative}
        T   (1, 2) {mustBeNonnegative}
        opt.psi {mustBeA(opt.psi, "ss")} = ...
                  basis_filter(1, 2*out(1), sampling_time=sys.Ts, ...
                               pole=0, type=2)
        opt.opt (1, 5) double = [1e-3, 400, 1e7, 50, 1]
        opt.resolve_bound (1, 1) {mustBeNonnegative} = 1.1
        opt.reconstruct {mustBeA(opt.reconstruct, "logical")} = true
    end

    %% Abbreviations

    % State dimensions
    lx  = size(sys.a, 1);     % of the original system 
    lxI = size(opt.psi.a, 1); % of the filter in the FB multiplier
    lxf = lxI + lx;           % of the filtered system

    % Descriptive names for several channel dimensions
    dis = inp(2); % generalized disturbance
    err = out(2); % error signal
    act = inp(3); % control signal
    mea = out(2); % measured output
    
    % Inner multiplier/scaling matrix dimension
    lm = size(opt.psi, 1);

    % This is to handle minimum and range dwell-time at the same time
    L  = max(T(T<Inf)); 

    %% Build outer factors
    
    % *Factors for the system LMIs*
    % Annihilator
    [~, ~, C3, D312] = ssdata(sys(end-mea+1:end, 1:end-act));
    V  = blkdiag(eye(lxI), null([C3, D312]));

    % Extended system without control channel
    sye = [eye(inp(1)+inp(2)); sys(1:end-mea, 1:end-act)]; 
    pnl = outerfactor_permutation(inp(1:2), out(1:2), 0);
    sye = blkdiag(opt.psi, eye(err+dis)) * sye(pnl, :); % Filtered system
    OY  = factor(sye);
    OX  = OY * V;

    % *Factors for the LMIs corresponding to the impulsive operator* 
    e = eye(out(1));    % An abbreviation
    z = zeros(out(1));  % Another abbreviation
    for i = T(1)+1:L+1 
       % Build lifted filter outer factor
       OZ{i} = factor(lifted_system(opt.psi, 2*out(1), lm, i));
       OZ{i} = OZ{i} * blkdiag(eye(lxI), kron(eye(i-1), [e; z]), [e; e]);   
   end
   if T(2) == Inf % Minimum dwell-time special case
       OM = factor(opt.psi);   
   end

    %% Define variables

    setlmis([]);
    
    % Performance level gamma
    [ga, ~,  sga] = lmivar(1, [  1, 1]); 
    
    % Lyapunov certificates
    [X, ~, sX] = lmivar(1, [lxf, 1]); % of primal system LMI
    [Y, ~, sY] = lmivar(1, [lxf, 1]); % of dual system LMI
    [Z, ~, sZ] = lmivar(1, [lxI, 1]); % of uncertainty LMI

    % Multiplier/scaling matrices
    [M, ~, sM] = lmivar(1, [lm, 1]);

    % *Inner terms*
    IX  = lmivar(3, blkdiag(sX, -sX, sM, zeros(err), -sga * eye(dis)));   
    IXr = blkdiag(zeros(2*lxf+lm), eye(err), zeros(dis));
    IY  = lmivar(3, blkdiag(sY, -sY, sM, zeros(err), -sga * eye(dis)));
    
    %% Constraints

    % *Positivity/coupling*
    O = eye(lxI, lxf);
    k = newlmi;
    lmiterm([-k, 1, 1, X],  1, 1);
    lmiterm([ k, 1, 1, Z], O', O);
    lmiterm([-k, 1, 2, Y],  1, 1);
    lmiterm([ k, 1, 2, Z], O', O);
    lmiterm([-k, 2, 2, Y],  1, 1);
    lmiterm([ k, 2, 2, Z], O', O);

    % *Primal system LMI*
    k = newlmi;
    lmiterm([k, 1, 1, IX], OX', OX);
    lmiterm([k, 1, 1,  0], OX' * IXr * OX);
    
    % *Dual system LMI (this is written in primal form)*
    k = newlmi;
    lmiterm([k, 1, 1, IY], OY', OY);
    
    % *Uncertainty LMI*
    for i = T(1)+1:L+1 
       IZI = lmivar(3, blkdiag(sZ, -sZ, kron(eye(i), sM)));
       k   = newlmi;
       lmiterm([-k, 1, 1, IZI], OZ{i}', OZ{i});   
   end
   if T(2) == Inf % Minimum dwell-time special case
       IZI = lmivar(sZ, -sZ, sM); 
       k   = newlmi;
       lmiterm([-k, 1, 1, IZI], OM', OM);   
   end

    %% Solve problem

    lmis = getlmis;
    ndec = decnbr(lmis) - 1;
 
    [gao, xfeas] = mincx(lmis, [1, zeros(1, ndec)], opt.opt);

    % Stop if no estimator reconstruction is desired
    if ~opt.reconstruct
        ga = sqrt(gao); % Get optimal performance gamma
        E  = [];
        return
    end

    %% Reconstruct estimator

    % Resolve LMIs with lower bound on determined performance level. 
    % This often yields numerically nicer matrices.
    if opt.resolve_bound ~= 1
        setlmis(lmis)

        k = newlmi;
        lmiterm([-k, 1, 1, ga], 1, 1);
        lmiterm([ k, 1, 1,  0], gao*opt.resolve_bound);

        lmis = getlmis;
        
        [gao, xfeas] = mincx(lmis, [1, zeros(1, ndec)], opt.opt);
    end
    ga = sqrt(gao); % Get optimal performance gamma


    % Get relevant variables
    X = dec2mat(lmis, xfeas, X);
    Y = dec2mat(lmis, xfeas, Y);
    M = dec2mat(lmis, xfeas, M);

    % Build extended Lyapunov certificate
    U  = Y - X;
    Xc = [X, U; U, -U];

    % Middle matrix
    P = blkdiag(Xc, -Xc, M, eye(err), -gao*eye(dis));

    % Build U, V and W (this is somewhat awkward/annoying)
    W = [sye.a, zeros(lxf), sye.b; zeros(lxf, 2*lxf+inp(1)+dis); ...
         eye(2*lxf), zeros(2*lxf, inp(1)+dis); ...
         sye.c, zeros(lm+err+dis, lxf), sye.d];
    U = [blkodiag(zeros(lxf, act), eye(lxf)); zeros(2*lxf+lm, lxf+act); ...
         blkodiag(-eye(act), zeros(dis, lxf))]';
    V = [blkdiag(zeros(lxf, lxI), C3), blkdiag(eye(lxf), D312)];

    % Apply elimination lemma from Helmersson to reconstruct the desired 
    % estimator
    E = mat2ss(elimi(P, U, V, W), lxf, sys.Ts);    
end


