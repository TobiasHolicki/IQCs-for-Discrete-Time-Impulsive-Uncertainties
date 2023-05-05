%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% File    : basis_filter.m                                                %
%                                                                         %
% Author  : Tobias Holicki                                                %
% Version : 03                                                            %
% Date    : 01.07.2022                                                    %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% For a given length and pole and input dimension function returns the 
% state space description of the basis transfer-matrix
%         (1, (s-a)/(s+a), ..., ((s-a)/(s+a))^\length) kron I_inpdim
% or of
%         (1, 1/(s+a), ..., 1/(s+a)^\length)^T kron I_inpdim.
% Adding an extra input argument results in the transposed system which is
% typically used for dual problems.
%
% ----- Input ---------------------------------------------------------- 
%       len - Basis length
%    inpdim - Dimension of the identity matrix
%      rest - Struct involving optional inputs
%          pole - Pole location (should be negative for continuous-time or 
%                 in norm smaller than for discrete-time 
%          type - Specifies if the first (1) or second (2) basis is used
% sampling_time - For returning a discrete-time filter if desired
%    primaldual - For returning the corresponding dual filter
% ----- Output ---------------------------------------------------------
%       sys - System as described above 
function [ sys ] = basis_filter(len, inpdim, rest)
    arguments
       len (1, 1) {mustBeInteger}
       inpdim (1, 1) {mustBeInteger}
       rest.pole (1, 1) double = -0.5
       rest.type (1, 1) {mustBeInteger} = 1
       rest.sampling_time (1, 1) double = 0 
       rest.primaldual {mustBeMember(rest.primaldual,["primal","dual"])}...
                        = 'primal'
    end
    
    % Abbreviation
    p = rest.pole;

    % Build state-space representations
    if rest.type == 1
        T = fliplr(eye(len));
        A = p * (eye(len) + 2* triu(ones(len), 1));
        B =  - sqrt(-2 * p) * ones(len, 1);
        C = [zeros(1, len); sqrt(-2 * p) * T * triu(ones(len), 0)];
        D = ones(len + 1, 1);
    elseif rest.type == 2
        A =  p * eye(len) + diag(ones(len-1,1), -1);
        B = double(1:len == 1)';
        C = [zeros(1, len); eye(len)];
        D = [1; zeros(len, 1)];
    else
        error('Unknown type')
    end
    
    % Incorporate repetitions
    A = kron(A, eye(inpdim));
    B = kron(B, eye(inpdim));
    C = kron(C, eye(inpdim));
    D = kron(D, eye(inpdim));

    if strcmp(rest.primaldual, 'dual') 
        sys = ss(A', C', B', D', rest.sampling_time);
    else
        sys = ss(A, B, C, D, rest.sampling_time);
    end

end