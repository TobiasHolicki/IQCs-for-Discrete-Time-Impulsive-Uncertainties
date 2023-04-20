%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% File     : design_es_svar.m                                             %
%                                                                         %
% Author   : Tobias Holicki                                               %
% Date     : 18.04.2023                                                   %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% This function implements the synthesis criteria for finding a 
% non-impulsive controller (an LTI estimator) such that the robust energy 
% gain of the closed-loop interconnection with the following impulsive 
% open-loop plant is minimized: 
% x(t+1) =  A x(t) +  Bd d(t)        x(tk+1) =  AJ x(tk) +  BJd d(tk)     
%   e(t) = Ce x(t) + Ded d(t) -u(t)    e(tk) = CJe x(tk) + DJed d(tk)-u(tk)
%    y(t) = Cy x(t) + Dyd d(t)         y(tk) = CJy x(tk) + DJyd d(tk)
% with (t(k))_k being a sequence that is merely known to satisfy a minimum 
% or range dwell-time condition.
%
% The synthesis criteria result from a generalization of the clock-based 
% analysis criteria for impulsive systems from [1] and from incorporating 
% slack-variables [2] similarly as has been done in [3].
%
% [1] C. Briat, Convex conditions for robust stability analysis and 
%     stabilization of linear aperiodic impulsive and sampled-data systems 
%     under dwell-time constraints, 2013
% [2] M.C. de Oliveira, J. Bernussou, J.C. Geromel, A new discrete-time
%     robust stability condition, 1999
% [3] J. C. Geromel, M. C. de Oliveira, J. Bernussou, Robust filtering for
%     discrete-time linear systems with parameter dependent Lyapunof
%     functions, 2002
%
% ----- Input ---------------------------------------------------------- 
% sys               - System describing the above flow component 
% sysJ              - System describing the above jump component
% inp               - Partition of the input signals [d, u]
% out               - Partition of the output signals [e, y]
% T                 - Assumed dwell-time conditions, i.e., 
%                           t(k+1) - t(k) - 1 \in [T(1), T(2)].
% opt               - Struct
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
function [ga, E] = design_es_svar(sys, sysJ, inp, out, T, opt)
    % Some sanity checks
    arguments
        sys {mustBeA(sys, "ss")}
        sysJ {mustBeA(sysJ, "ss")}
        inp (1, 2) {mustBeInteger, mustBeNonnegative}
        out (1, 2) {mustBeInteger, mustBeNonnegative}
        T   (1, 2) {mustBeNonnegative}
        opt.opt (1, 5) double = [1e-3, 400, 1e7, 50, 1]
        opt.resolve_bound (1, 1) {mustBeNonnegative} = 1.1
        opt.reconstruct {mustBeA(opt.reconstruct, "logical")} = true
    end
    
    %% Abbreviations
    
    % State dimensions and descriptive names for several channel dimensions
    lx  = size(sys.a, 1); % State dimension
    dis = inp(1); % generalized disturbance
    err = out(1); % error signal
    act = inp(2); % control signal
    mea = out(2); % measured output
    
    %% Build outer factors
    % This is somewhat tiresome here

    Pl = blkdiag(eye(2*lx, 4*lx+dis)', eye(err));
    Pr = [zeros(2*lx+dis, 2*lx), eye(2*lx+dis), zeros(2*lx+dis, err)];
    
    [A, B, C, D] = sssdata(sys, inp, out);
    O  = [A, A, B{1}; A, A, B{1}; C{1}, C{1}, D{1, 1}];
    Ol = [zeros(lx, lx+act); blkdiag(eye(lx), -eye(act))];
    Or = [eye(lx, 2*lx+dis); C{2}, C{2}, D{2,1}];
    
    
    [AJ, BJ, CJ, DJ] = sssdata(sysJ, inp, out);
    OJ  = [AJ, AJ, BJ{1}; AJ, AJ, BJ{1}; CJ{1}, CJ{1}, DJ{1, 1}];
    OJl = [zeros(lx, lx+act); blkdiag(eye(lx), -eye(act))];
    OJr = [eye(lx, 2*lx+dis); CJ{2}, CJ{2}, DJ{2,1}];
    
    %% Define variables

    setlmis([]);
    
    [ga, ~,  sga] = lmivar(1, [1, 1]);  % Performance level gamma
    
    % Lyapunov certificate
    for i = 1 : T(2)+1
       [X{i}, ~, sX{i}] = lmivar(1, [2*lx, 1]); 
    end
    
    % Components of the slack variables for the flow component
    for i = 1 : T(2)
       [G{i}, ~, sG{i}] = lmivar(2, [lx, lx]);
       [H{i}, ~, sH{i}] = lmivar(2, [lx, lx]);
    end
    % S - G = SJ - GJ =: C is required to be constant
    [C, ~, sC] = lmivar(2, [lx, lx]);
    
    % Components of the slack variables for the jump component
    for i = T(1) : T(2)
        [GJ{i}, ~, sGJ{i}] = lmivar(2, [lx, lx]);
        [HJ{i}, ~, sHJ{i}] = lmivar(2, [lx, lx]);
    end
    
    % KLMN matrices
    [KLMN, ~, ~] = lmivar(2, [lx+act, lx+mea]);
    
    % *Inner terms*
    for i = 1 : T(2)
        IX{i} = lmivar(3, blkdiag(sX{i+1}, -sX{i}, -sga*eye(dis), ...
                                  zeros(err)));
        IG{i} = lmivar(3, blkdiag([sH{i}, sH{i}; sG{i}, sG{i}], ...
                                  zeros(2*lx+dis+err)));
        IA{i} = lmivar(3, blkdiag(sH{i}, sG{i}, zeros(err)));
    end
    IC = lmivar(3, blkdiag([zeros(lx, 2*lx); sC, zeros(lx)], ...
                           zeros(2*lx+dis+err)));
    IXr = blkdiag(zeros(4*lx+dis), eye(err));
    IAr = blkdiag(zeros(2*lx), eye(err));
    
    for i = T(1) : T(2)
        IXJ{i} = lmivar(3, blkdiag(sX{1}, -sX{i+1}, -sga*eye(dis), ...
                                  zeros(err)));
        IGJ{i} = lmivar(3, blkdiag([sHJ{i}, sHJ{i}; sGJ{i}, sGJ{i}], ...
                                  zeros(2*lx+dis+err)));
        IAJ{i} = lmivar(3, blkdiag(sHJ{i}, sGJ{i}, zeros(err)));
    end
    
    %% Constraints

    % *Positivity/coupling*
    for i = 1 : T(2)+1
        k = newlmi;
        lmiterm([-k, 1, 1, X{i}], 1, 1);
    end
    
    % *Flow LMI*
    for i = 1 : T(2)
        k = newlmi;
        % Diagonal term (Xb(k+1) - Gb - Gb', -Xb(k), -ga^2I, -I)
        lmiterm([k, 1, 1, IX{i}], 1, 1);
        lmiterm([k, 1, 1, IG{i}], -1, 1, 's');
        lmiterm([k, 1, 1,    IC], -1, 1, 's');
        lmiterm([k, 1, 1,     0], -IXr);
        % Remaining term involving slack variables
        lmiterm([k, 1, 1, IA{i}], Pl, O * Pr, 's');
        lmiterm([k, 1, 1,     0], (Pl * IAr * O * Pr) + ...
                                  (Pl * IAr * O * Pr)');
        % KLMN variables
        lmiterm([k, 1, 1, KLMN], Pl * Ol, Or * Pr, 's');
    end
    
    % *Jump LMI*
    for i = T(1) : T(2)
        k = newlmi;
        % Diagonal term (Xb(0) - Gb - Gb', -Xb(k), -ga^2I, -I)
        lmiterm([k, 1, 1, IXJ{i}], 1, 1);
        lmiterm([k, 1, 1, IGJ{i}], -1, 1, 's');
        lmiterm([k, 1, 1,     IC], -1, 1, 's');
        lmiterm([k, 1, 1,      0], -IXr);
        % Remaining term involving slack variables
        lmiterm([k, 1, 1, IAJ{i}], Pl, OJ * Pr, 's');
        lmiterm([k, 1, 1,      0], (Pl * IAr * OJ * Pr) + ...
                                   (Pl * IAr * OJ * Pr)');
        % KLMN variables
        lmiterm([k, 1, 1, KLMN], Pl * OJl, OJr * Pr, 's');
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
    C    = dec2mat(lmis, xfeas, C);
    KLMN = dec2mat(lmis, xfeas, KLMN);
    
    % Undo convexifying parameter transform
    E = mat2ss(blkdiag(C, eye(act)) \ KLMN, lx, sys.Ts);
    
end


