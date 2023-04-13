%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% File     : ana_perf_clock.m                                             %
%                                                                         %
% Author   : Tobias Holicki                                               %
% Date     : 13.04.2023                                                   %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% This function implements the clock-based H_infty analysis criteria for a 
% linear discrete-time impulsive system
%   x(t+1) = A x(t) + B d(t)            x(t(k)+1) = AJ x(t(k)) + BJ d(t(k))
%     e(t) = C x(t) + D d(t)              e(t(k)) = CJ x(t(k)) + DJ d(t(k))
% where the sequence t(k) satisfies either a minimum or range dwell-time
% condition. These criteria are a variation of the ones in [1].
%
% [1] C. Briat, Convex conditions for robust stability analysis and 
%     stabilization of linear aperiodic impulsive and sampled-data systems 
%     under dwell-time constraints, 2013
%
% ----- Input ---------------------------------------------------------- 
%   sysF     - State-space model for the system's flow component
%   sysJ     - State-space model for the system's jump component
%   T        - Assumed dwell-time conditions, i.e., 
%                       t(k+1) - t(k) - 1 \in [T(1), T(2)] 
%              where T(2) can be Inf.
%   opt      - Struct
%      .opt  - LMIlab solver options
% ----- Output ---------------------------------------------------------
%   gao      - Numerically determined upper bound on the system's energy 
%              gain.
function gao = ana_perf_clock(sysF, sysJ, T, opt)
    % Some sanity checks
    arguments
        sysF {mustBeA(sysF, "ss")}
        sysJ {mustBeA(sysJ, "ss")}
        T    (1, 2) {mustBeNonnegative}
        opt.opt (1, 5) double = [1e-3, 400, 1e9, 50, 1]
    end
    %% Abbreviations
    
    lx = size(sysF.a, 1); % State dimension
    L  = max(T(T<Inf));   % This is to handle minimum and range dwell-time
                          % conditions at the same time

    dis = size(sysF, 2);  % Dimension of generalized disturbance
    err = size(sysF, 1);  % Dimension of error signal

    %% Define variables

    setlmis([]);
    
    % Upper bound on energy gain
    [~, ~, sga] = lmivar(1, [1, 1]);

    % Lyapunov certificates
    for i = 1 : L+1
        [X{i}, ~,  sX{i}] = lmivar(1, [lx, 1]);    
    end
    
    %% Constraints
  
    % *Positivity*
    for i = 1 : length(X)
        k = newlmi;
        lmiterm([-k, 1, 1, X{i}], 1, 1);
    end

    % *Remaining LMIs*
    % Outer factors
    [OXF, ~] = outerfactor(sysF, dis, err, 'ana');
    [OXJ, ~] = outerfactor(sysJ, dis, err, 'ana');

    % Flow LMIs
    IXr = blkdiag(zeros(2*lx), eye(err), zeros(dis));
    for i = 1 : L
        IX = lmivar(3, blkdiag(sX{i+1}, -sX{i}, zeros(err), ...
                               -sga*eye(dis))); % Inner term
        k  = newlmi;
        lmiterm([k, 1, 1, IX], OXF', OXF);
        lmiterm([k, 1, 1,  0], OXF' * IXr * OXF);
    end
    if T(2) == Inf % Minimum dwell-time special case
        IX = lmivar(3, blkdiag(sX{end}, -sX{end}, zeros(err), ...
                               -sga*eye(dis))); % Inner term
        k  = newlmi;
        lmiterm([k, 1, 1, IX], OXF', OXF);
        lmiterm([k, 1, 1,  0], OXF' * IXr * OXF);
    end
    % Jump LMIs
    for i = T(1) : L
        IX = lmivar(3, blkdiag(sX{1}, -sX{i+1}, zeros(err), ...
                               -sga*eye(dis))); % Inner term
        k  = newlmi;
        lmiterm([k, 1, 1, IX], OXJ', OXJ);
        lmiterm([k, 1, 1,  0], OXJ' * IXr * OXJ);
    end
    
    %% Try to solve LMI problem

    lmis = getlmis;
    ndec = decnbr(lmis) - 1;

    [gao, ~] = mincx(lmis, [1, zeros(1, ndec)], opt.opt);
    gao      = sqrt(gao);        
end
    
    
    
    
