%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% File     : ana_iqc.m                                                    %
%                                                                         %
% Author   : Tobias Holicki                                               %
% Date     : 13.04.2023                                                   %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% This function implements stability analysis criteria for the feedback 
% interconnection
%  x(t+1) = A x(t) + B w(t),  z(t) = C x(t) + D w(t),     w(t) = Del(z)(t)
% with impulsive operator
%  Del(z)(t) = z(t)  if t = t(k) for some k  and  Del(z)(t) = 0  otherwise
% where (t(k))_k is a sequence that is merely known to satisfy a minimum or
% range dwell-time condition.
%
% The criteria result from a generalization of the IQC theorem from [1] and
% make use of some of the available ideas on discrete-time impulsive
% systems. This function implements three variants based on different IQCs
% motivated by:
% - The classical discrete-time lifting procedure.
% - The clock-based approach from [2].
% - The path-based approach from [3].
%
% [1] C. W. Scherer, J. Veenman, Stability analysis by dynamic dissipation 
%     inequalities: On merging frequency-domain techniques with time domain
%     conditions, 2018
% [2] C. Briat, Convex conditions for robust stability analysis and 
%     stabilization of linear aperiodic impulsive and sampled-data systems 
%     under dwell-time constraints, 2013
% [3] W. Xiang, H.-D. Tran, T. T. Johnson, Nonconservative lifted convex 
%     conditions for stability of discrete-time switched systems under 
%     minimum dwell-time constraint, 2019
% ----- Input ---------------------------------------------------------- 
%   sys      - System describing the linear part of the interconnection
%   T        - Assumed dwell-time conditions, i.e., 
%                       t(k+1) - t(k) - 1 \in [T(1), T(2)] 
%              where T(2) can be Inf.
%   opt      - Struct
%      .type - The type of IQCs conditions that are considered.
%              Implemented are the lifting, the clock-based and a variant
%              of the path-based approach which involves lifted matrices.
%      .psi  - The dynamic filter involved in the IQC
%      .opt  - LMIlab solver options
%      .L    - Path length in the path-based approach
% ----- Output ---------------------------------------------------------
%   t        - If this scalar is negative the solver has found a solution
%              to the involved LMI problem and the impulsive system is
%              stable.
function t = ana_iqc(sys, T, opt)
    % Some sanity checks
    arguments
        sys        {mustBeA(sys, "ss")}
        T   (1, 2) {mustBeNonnegative}
        opt.type  char {mustBeMember(opt.type, ...
                                   {'Lifting','Clock','Path'})} = 'Lifting'
        opt.opt (1, 5) double = [1e-3, 400, 1e9, 50, 1]
        opt.psi {mustBeA(opt.psi, "ss")} = ...
                  basis_filter(1, sum(size(sys)), sampling_time=sys.Ts, ...
                               pole=0, type=2)
        opt.L (1, 1) {mustBeInteger, mustBeNonnegative} = max(T(T<Inf)) + 1
    end
    %% Abbreviations
    
    % State dimensions
    lx  = size(sys.a, 1);      % of the original system 
    lxI = size(opt.psi.a, 1);  % of the filter in the FB multiplier
    lxf = lx + lxI;            % of the filtered system
    
    % Inner multiplier/scaling matrix dimension
    lmI = size(opt.psi, 1);
    
    % Filter
    ps = opt.psi;
    
    % This is to handle minimum and range dwell-time at the same time
    L  = max(T(T<Inf));                  

    % Generate admissible and postadmissible paths for the path approach.
    if strcmp(opt.type, 'Path')
        [aP, paP] = paths(T, opt.L);
    end
    
    %% Build outer factors
    
    % For the system LMI
    sye = ps * [sys; eye(size(sys, 2))]; % Filtered system
    OX  = factor(sye);

    %% Define variables

    setlmis([]);
    
    % Lyapunov certificates
    [X, ~,  sX] = lmivar(1, [lxf, 1]); % for the filtered system LMI
    
    % Certificates and multiplier matrices for the impulsive part
    [Z, sZ, M, sM] = define_variables();
    
    % *Inner terms*
    IX  = lmivar(3, blkdiag(sX, -sX, sM)); % System LMI

    %% Constraints
  
    % *Positivity/Coupling*
    OP = eye(lxI, lxI + lx);
    for i = 1 : length(Z)
        k = newlmi;
        lmiterm([-k, 1, 1,    X],   1,  1);
        lmiterm([ k, 1, 1, Z{i}], OP', OP);
    end

    % *System LMI*
    k = newlmi;
    lmiterm([k, 1, 1, IX], OX', OX);

    % *Constraints for the impulsive part*
    add_multiplier_constraints(Z, sZ, M, sM);
    
    %% Try to solve LMI problem

    lmis = getlmis;
    
    [t, ~] = feasp(lmis, opt.opt);

    %% Auxilliary functions

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Function : define_variables
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % This functions introduces the variables appearing in the three
    % approaches.
    function [Z, sZ, M, sM] = define_variables()
        % Lyapunov certificates
        if strcmp(opt.type, 'Lifting')
            [Z{1}, ~, sZ{1}] = lmivar(1, [lxI, 1]);
        elseif strcmp(opt.type, 'Clock')
            for ii = 1 : L+1
                [Z{ii}, ~, sZ{ii}] = lmivar(1, [lxI, 1]);
            end
        elseif strcmp(opt.type, 'Path')
            for ii = 1 : length(aP)
                [Z{ii}, ~, sZ{ii}] = lmivar(1, [lxI, 1]);
            end
        end

        % Multiplier matrix
        [M, ~, sM] = lmivar(1, [lmI, 1]);
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Function : add_multiplier_constraints
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % This functions produces the LMI constraints in the respective
    % approaches.
    function add_multiplier_constraints(Z, sZ, M, sM)
        e = eye(size(sys, 1));     % An abbreviation
        z = zeros(size(sys, 1));   % Another abbreviation
        
        if strcmp(opt.type, 'Lifting')
           [a, b, c, d] = ssdata(ps); % Filter data
           for ii = T(1)+1:L+1 
               % Build lifted filter outer factor
               [A, B, C, D] = lifted_matrices(a, b, c, d, ii);
               OZ  = [A, B; eye(size([A, B])); C, D];
               OZ  = OZ * ...
                        blkdiag(eye(lxI), kron(eye(ii-1), [e; z]), [e; e]);
               % Build inner term
               IZI = lmivar(3, blkdiag(sZ{1}, -sZ{1}, kron(eye(ii), sM)));
               k   = newlmi;
               lmiterm([-k, 1, 1, IZI], OZ', OZ);
           end
           if T(2) == Inf % Minimum dwell-time special case
                OZ  = factor(ps);                % Outer factor
                IZI = lmivar(sZ{1}, -sZ{1}, sM); % Inner term
                k   = newlmi;
                lmiterm([-k, 1, 1, IZI], OZ', OZ);
           end
        elseif strcmp(opt.type, 'Clock')
            % Build outer factors
            OZ  = factor(ps);
            OZF = OZ * blkdiag(eye(lxI), [e; z]);
            OZJ = OZ * blkdiag(eye(lxI), [e; e]);
            % Flow
            for ii = 1 : L
                IZI = lmivar(3, blkdiag(sZ{ii+1}, -sZ{ii}, sM));
                k   = newlmi;
                lmiterm([-k, 1, 1, IZI], OZF', OZF);
            end
            if T(2) == Inf % Minimum dwell-time special case
               IZI = lmivar(3, blkdiag(sZ{end}, -sZ{end}, sM));
               k   = newlmi;
               lmiterm([-k, 1, 1, IZI], OZF', OZF);
            end
            % Jump
            for ii = T(1) : L
               IZI = lmivar(3, blkdiag(sZ{1}, -sZ{ii+1}, sM));
               k   = newlmi;
               lmiterm([-k, 1, 1, IZI], OZJ', OZJ);
            end
        elseif strcmp(opt.type, 'Path')
            oL = opt.L; % Abbreviation
            [a, b, c, d] = ssdata(ps); % Filter data
            [A, B, C, D] = lifted_matrices(a, b, c, d, oL);
            
           for ii = 1:length(aP)
               % Outer factor
               OZ = mat2cell([ones(1, oL); aP{ii}'], 2, ones(1, oL));
               OZ = blkdiag(eye(lxI), ...
                            kron(blkdiag(OZ{:}), eye(size(sys, 1))));
               OZ = [A, B; eye(size([A, B])); C, D] * OZ;
               % Inner term
               for jj = 1 : length(paP{ii})
                   IZI = lmivar(3, blkdiag(sZ{paP{ii}{jj}}, -sZ{ii}, ...
                                           kron(eye(oL), sM)));
                   k   = newlmi;
                   lmiterm([-k, 1, 1, IZI], OZ', OZ);
               end
           end
        end       
    end % add_multiplier_constraints
    
end
    
    
    
    
