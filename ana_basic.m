%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% File     : ana_basic.m                                                  %
%                                                                         %
% Author   : Tobias Holicki                                               %
% Date     : 13.04.2023                                                   %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% This function implements stability analysis criteria for a linear
% discrete-time impulsive system
%   x(t+1) = A x(t)      x(t(k)+1) = AJ x(t(k))
% where the sequence t(k) satisfies either a minimum or range dwell-time
% condition.
% One can choose among three types of criteria:
% - The classical one based on the discrete-time lifting procedure.
% - A variant of the clock-based approach from [1].
% - A variant of the path-based approach from [2].
%
% [1] C. Briat, Convex conditions for robust stability analysis and 
%     stabilization of linear aperiodic impulsive and sampled-data systems 
%     under dwell-time constraints, 2013
% [2] W. Xiang, H.-D. Tran, T. T. Johnson, Nonconservative lifted convex 
%     conditions for stability of discrete-time switched systems under 
%     minimum dwell-time constraint, 2019
% ----- Input ---------------------------------------------------------- 
%   A        - Describing matrix of the system's flow component
%   AJ       - Describing matrix of the system's jump component
%   T        - Assumed dwell-time conditions, i.e., 
%                       t(k+1) - t(k) - 1 \in [T(1), T(2)] 
%              where T(2) can be Inf.
%   opt      - Struct
%      .type - The type of analysis conditions that are considered.
%              Implemented are the lifting, the clock-based and a variant
%              of the path-based approach which involves lifted matrices.
%      .opt  - LMIlab solver options
%      .L    - Path length in the path-based approach
% ----- Output ---------------------------------------------------------
%   t        - If this scalar is negative the solver has found a solution
%              to the involved LMI problem and the impulsive system is
%              stable.
function [t, X] = ana_basic(A, AJ, T, opt)
    % Some sanity checks
    arguments
        A  (:, :) {mustBeNumeric}
        AJ (:, :) {mustBeNumeric}
        T  (1, 2) {mustBeNonnegative}
        opt.type  char {mustBeMember(opt.type, ...
                                   {'Lifting','Clock','Path'})} = 'Lifting'
        opt.opt (1, 5) double = [1e-3, 400, 1e9, 50, 1]
        opt.L   (1, 1) {mustBeInteger, mustBeNonnegative} = max(T(T<Inf))+1
    end
    %% Abbreviations
    
    lx = size(A, 1);    % State dimension
    L  = max(T(T<Inf)); % This is to handle minimum and range dwell-time
                        % conditions at the same time
    
    % Generate admissible and postadmissible paths for the path approach.
    if strcmp(opt.type, 'Path')
        [aP, paP] = paths(T, opt.L);
    end

    %% Define variables

    setlmis([]);
    
    % Lyapunov certificates
    if strcmp(opt.type, 'Lifting')
        [X{1}, ~,  sX] = lmivar(1, [lx, 1]);    
    elseif strcmp(opt.type, 'Clock')
        for i = 1 : L+1
            [X{i}, ~,  sX{i}] = lmivar(1, [lx, 1]);    
        end
    elseif strcmp(opt.type, 'Path')
        for i = 1 : length(aP)
            [X{i}, ~,  sX{i}] = lmivar(1, [lx, 1]);    
        end
    end
    
    %% Constraints
  
    % *Positivity*
    for i = 1 : length(X)
        k = newlmi;
        lmiterm([-k, 1, 1, X{i}], 1, 1);
    end

    % *Remaining LMIs*
    if strcmp(opt.type, 'Lifting')
        IX = lmivar(3, blkdiag(sX, -sX)); % Inner term
        for i = T(1):L
            OX = [AJ * A^i; eye(lx)]; % Outer factor
            k  = newlmi;
            lmiterm([k, 1, 1, IX], OX', OX);
        end
        if T(2) == Inf % Minimum dwell-time special case
            OX = [A; eye(lx)]; % Outer factor
            k  = newlmi;
            lmiterm([k, 1, 1, IX], OX', OX);
        end
    elseif strcmp(opt.type, 'Clock')
        % Outer factors
        OXF = [ A; eye(lx)];
        OXJ = [AJ; eye(lx)];
        % Flow LMIs
        for i = 1 : L
            IX = lmivar(3, blkdiag(sX{i+1}, -sX{i})); % Inner term
            k  = newlmi;
            lmiterm([k, 1, 1, IX], OXF', OXF);
        end
        if T(2) == Inf % Minimum dwell-time special case
            IX = lmivar(3, blkdiag(sX{end}, -sX{end})); % Inner term
            k  = newlmi;
            lmiterm([k, 1, 1, IX], OXF', OXF);
        end
        % Jump LMIs
        for i = T(1) : L
            IX = lmivar(3, blkdiag(sX{1}, -sX{i+1})); % Inner term
            k  = newlmi;
            lmiterm([k, 1, 1, IX], OXJ', OXJ);
        end
    elseif strcmp(opt.type, 'Path')
        for i = 1 : length(aP)
            OX = [Av(aP{i}); eye(lx)]; % Outer factor 
            for j = 1 : length(paP{i})
                IX = lmivar(3, blkdiag(sX{paP{i}{j}}, -sX{i}));
                k  = newlmi;
                lmiterm([k, 1, 1, IX], OX', OX);
            end
        end
    end
    
    %% Try to solve LMI problem

    lmis = getlmis;
    
    [t, ~] = feasp(lmis, opt.opt);

    %% Auxilliary functions

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Function : Av
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % This returns the product of the describing system matrices along the
    % path p.
    function tA = Av(p)
       tA = eye(lx);
       for ii = 1 : length(p)
           if p(ii) == 0 % Flow
               tA = A * tA;
           else % Jump
               tA = AJ * tA;
           end
       end
    end     
end
    
    
    
    
