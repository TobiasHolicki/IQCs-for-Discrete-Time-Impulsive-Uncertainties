%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% File     : paths.m                                                      %
%                                                                         %
% Author   : Tobias Holicki                                               %
% Date     : 19.07.2022                                                   %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% This function determines the set of admissible impulse paths of length L
% and the corresponding set of postadmissble impulse paths for sequences of
% impulse instants satisfying a minimum or range dwell-time condition.
% Both sets are encoded as cells.
% ----- Input ---------------------------------------------------------- 
%   T        - Assumed dwell-time conditions, i.e., 
%                       t(k+1) - t(k) - 1 \in [T(1), T(2)] 
%              where T(2) can be Inf.
%   L        - Considered path length.
% ----- Output ---------------------------------------------------------
%   aP       - The set of admissible paths of length L for dwell-time T.
%   paP      - The corresponding set of postadmissble paths.
function [aP, paP] = paths(T, L)
    % Some sanity checks
    arguments
        T (1, 2) {mustBeNonnegative}
        L (1, 1) {mustBeInteger, mustBeNonnegative} = max(T(T<Inf)) + 1
    end
    
    %% Generate all admissible paths of length L for arbitrary dwell-time
    % This is from https://de.mathworks.com/matlabcentral/answers/
    %                       384426-all-possible-combinations-of-0-s-and-1-s
    % and is somewhat expensive and of course more than we actually need.

    C      = cell(1, L);
    [C{:}] = ndgrid(0:1);
    C      = reshape(cat(L+1, C{:}), [], L)';
    
    %% Generate all admissible paths of length L for dwell-time T
    % We do one and no jumps seperately.
    
    e   = @(i) double(1:L == i)'; %  Standard unit vectors 
    ind = 1; % Some index
    
    % No jumps
    if T(2) >= L % Only happens if L is not large enough or for minimum DT
        aP{ind} = zeros(L, 1);
        ind     = ind+1;
    end
            
    % One jump
    for i = 1 : L
        % We omit those that have to many zeros above or below the 1 entry
        if i-1 <= T(2) && L-i <= T(2)
            aP{ind} = e(i);
            ind     = ind+1;
        end
    end
    
    % Two or more jumps
    for i = 1 : length(C)
         ii = find(C(:, i) == 1); % Find impulse instants
         if length(ii) >= 2
             % Determine wether instants satisfy the dwell-time condition
             dtt = (ii(2:end) - ii(1:end-1) - 1 <= T(2)) .* ...
                   (ii(2:end) - ii(1:end-1) - 1 >= T(1));
             if sum(dtt) == length(ii)-1 && ii(1)-1 <= T(2) && ...
                L-ii(end) <= T(2)
                 aP{ind} = C(:, i);
                 ind     = ind+1;
             end
         end
    end         

    %% Find postadmissible paths

    for i = 1 : length(aP)
        ind = 1;
        for j = 1 : length(aP)
            % Append candidate path vertically and from below
            vact = [aP{i}; aP{j}];
            % Find impulse instants
            ii = find(vact == 1);
            % Determine wether instants satisfy the dwell-time condition
            dtt = (ii(2:end) - ii(1:end-1) - 1 <= T(2)) .* ...
                  (ii(2:end) - ii(1:end-1) - 1 >= T(1));
            if sum(dtt) == length(ii)-1 && ii(1)-1 <= T(2) && ...
                L-ii(end) <= T(2)
                paP{i}{ind} = j;
                ind         = ind+1;
            end
        end
    end
end
