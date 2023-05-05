%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% File    : disturb.m                                                     %
%                                                                         %
% Author  : Tobias Holicki                                                %
% Version : 03                                                            %
% Date    : 18.03.2020                                                    %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Replace all singular values below 'border' of a matrix 'A' with 'new'. 
%
% ----- Input ---------------------------------------------------------- 
%       A - Given matrix
%  border - All singular values below border are considered
%     new - New value for all considered singular values
% ----- Output ---------------------------------------------------------
%       B - Resulting new matrix
%
function B = disturb(A, border, new)

% Decompose A with SVD
[U, S, V] = svd(A);

% Find singular values that are less then border
S      = diag(S);
ind    = find(S < border); 
% Replace them with new
S(ind) = new * ones(size(ind));
S      = diag(S);
% Build new matrix 
B      = U * S * V';

end

