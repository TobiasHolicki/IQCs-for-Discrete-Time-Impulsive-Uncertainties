%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% File     : lifted_matrices.m                                            %
%                                                                         %
% Author   : Tobias Holicki                                               %
% Version  : 01                                                           %
% Date     : 14.07.2022                                                   %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% This function generates the lifted system  matrices corresponding to the 
% standard discrete-time lifting procedure of an LTI system with describing 
% matrices A, B, C and D and for some horizon h.
%
% See for example 
% [Chen, T. and Francis, B.A. (1995). Optimal Sampled-Data Control Systems]
%
function [hA, hB, hC, hD] = lifted_matrices(A, B, C, D, h)
    % Some sanity checks
    arguments
        A double
        B double
        C double
        D double
        h (1, 1) {mustBeInteger}
    end

    n = size(A, 1);             % Size of A matrix
    J = diag(ones(h-1, 1), -1); % Appearing Jordan block
    
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
    hB = AA(end-n+1:end, :)*kron(eye(h), B);
    hC = kron(eye(h), C) * AA(:, 1:n);
    hD = kron(eye(h), D) + kron(J, C) * AA * kron(eye(h), B);
end

