%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% File     : lifted_matrices.m                                            %
%                                                                         %
% Author   : Tobias Holicki                                               %
% Version  : 01                                                           %
% Date     : 14.07.2022                                                   %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% This function generates a lifted version hsys of the discrete-time system
% sys whose input and output signals have the dimensions given in the 
% vectors inp and out. 
% E.g., if sys has the inputs d and u, then hsys will have the lifted 
% inputs hd and hu.
%
% The lifting procedure is described, e.g., in 
% [Chen, T. and Francis, B.A. (1995). Optimal Sampled-Data Control Systems]
%
function [hsys, hinp, hout] = lifted_system(sys, inp, out, h)
    % Some sanity checks
    arguments
        sys {mustBeA(sys, 'ss')}
        inp (1, :) {mustBeInteger}
        out (1, :) {mustBeInteger}
        h   (1, 1) {mustBeInteger}
    end
    if sum(inp) ~= size(sys, 2)
        error('Number of inputs and input partition do not match')
    elseif sum(out) ~= size(sys, 1)
        error('Number of outputs and output partition do not match')
    end

    % Some abbreviations
    n  = size(sys.a, 1);         % Size of system A matrix
    li = length(inp);            % Number of inputs
    lo = length(out);            % Number of outputs
    J  = diag(ones(h-1, 1), -1); % Appearing Jordan block
    
    % System data
    [A, B, C, D] = sssdata(sys, inp, out);
    
    % Matrix of powers of A. This is actually the same as:
    %     inv(eye(h*n) - kron(J, A));
    AA = cell(h, h);
    for i = 0 : h-1
        for j = 0 : h-1
            if j > i
                AA{i+1, j+1} = zeros(n);
            else
                AA{i+1, j+1} = A^(i-j);
            end
        end
    end
    AA = cell2mat(AA);
    
    % Lifted matrices
    hA = A^h;
    hB = cell(1, li);
    hC = cell(lo, 1);
    hD = cell(lo, li);
    
    for i = 1 : li
        hB{i} = AA(end-n+1:end, :)*kron(eye(h), B{i}); 
    end
    for j = 1 : lo
        hC{j} = kron(eye(h), C{j}) * AA(:, 1:n); 
    end
    for i = 1 : li
        for j = 1 : lo
            hD{j, i} = kron(eye(h), D{j, i}) + ...
                       kron(J, C{j}) * AA * kron(eye(h), B{i}); 
        end
    end
    
    hsys = ss(hA, cell2mat(hB), cell2mat(hC), cell2mat(hD));
    hinp = h * inp;
    hout = h * out;
end

