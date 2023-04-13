%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% File    : test_stab.m                                                   %
%                                                                         %
% Author  : Tobias Holicki                                                %
% Date    : 13.04.2023                                                    %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% This script realizes the example in Section 2.5.1 of [1].
% Its intention is to compare several stability tests for discrete-time 
% impulsive systems for a simple numerical example. These tests are based 
% on:
% - The classical discrete-time lifting procedure
% - A variation of the clock based approach from [2].
% - A variation of the path based approach from [3].
% - A variation of the IQC theorem from [4] with an IQC inspired by the
%   lifting approach.
% - A variation of the IQC theorem from [4] with an IQC inspired by the
%   clock-based approach
%
% [1] T. Holicki, C. W. Scherer, IQC Based Analysis and Estimator Design
%     for Discrete-Time Systems Affected by Impulsive Uncertainties, 2023
% [2] C. Briat, Convex conditions for robust stability analysis and 
%     stabilization of linear aperiodic impulsive and sampled-data systems 
%     under dwell-time constraints, 2013
% [3] W. Xiang, H.-D. Tran, T. T. Johnson, Nonconservative lifted convex 
%     conditions for stability of discrete-time switched systems under 
%     minimum dwell-time constraint, 2019
% [4] C. W. Scherer, J. Veenman, Stability analysis by dynamic dissipation
%     inequalities: On merging frequency-domain techniques with time domain
%     conditions, 2018


% Clean up
clc
clear

% Addpath for auxiliary functions 
addpath(genpath('AuxiliaryFunctions'));

%% System data

% Describing matrices of the considered impulsive system. Actually, this is
% the discretization of a nice example by [W. P. Dayawansa, C. F. Martin, 
% Dynamical systems which undergo switching 1999].
AJ  = [0, -10; 0.1, 0];
A   = @(beta) eye(2) + beta/100 * [-1, -1; 1, -1];

% Linear part of representation in terms of feedback interconnection
sys = @(beta) ss(A(beta), AJ - A(beta), eye(2), zeros(2), 0.5);


% Dwell-times
T = {[5, 5], [5, 6], [5, 7], [5, 8], [5, 9], [6, 6], [6, 7], [6, 8], ...
     [6, 9], [7, 7], [7, 8], [7, 9]};

% Sampling-time parameter range
beta = linspace(0.1, 5, 40);

% Initialization
e = zeros(length(beta), length(T), 4);

%% Perform several stability tests

% Classical lifting
tic
for i = 1 : length(beta)
    parfor j = 1 : length(T)
        e(i, j, 1) = ana_basic(A(beta(i)), AJ, T{j}, type='Lifting');
    end
end
toc

% Clock based
tic
for i = 1 : length(beta)
    parfor j = 1 : length(T)
        e(i, j, 2) = ana_basic(A(beta(i)), AJ, T{j}, type='Clock');
    end
end
toc

% Path based
tic
for i = 1 : length(beta)
    parfor j = 1 : length(T)
        e(i, j, 3) = ana_basic(A(beta(i)), AJ, T{j}, type='Path', L=11);
    end
end
toc

% IQC and lifting
tic
for i = 1 : length(beta)
    parfor j = 1 : length(T)
        e(i, j, 4) = ana_iqc(sys(beta(i)), T{j}, type='Lifting');
    end
end
toc

% IQC and clock
tic
for i = 1 : length(beta)
    parfor j = 1 : length(T)
        e(i, j, 5) = ana_iqc(sys(beta(i)), T{j}, type='Clock');
    end
end
toc

% IQC and path
% Uncomment the following only if you have a lot of time (4.5 hours on my 
% computer) ... The involved LMIs also seem to be cause numerical issues.
% tic
% for i = 1 : length(beta)
%     parfor j = 1 : length(T)
%         e(i, j, 6) = ana_iqc(sys(beta(i)), T{j}, type='Path', L=10);
%     end
% end
% toc

%% Plot stuff


% Find parameter and dwell-time combinations for which stability could be
% guaranteed
eps = 1e-9;
for l = 1 : length(e(1, 1, :))
    x{l} = [];
    y{l} = [];
    for j = 1 : length(T)
        f = find(e(:, j, l) < -eps);
        y{l} = [y{l}, beta(f)];
        x{l} = [x{l}, j * ones(1, length(f))];
    end
end

% Comparison with path based approach: The latter approach yields indeed 
% the least conservatism
figure
hold on; box on
scatter(x{1}, y{1}, 'o', 'LineWidth', 1); 
scatter(x{2}, y{2}, 'o') % Essentially identical to the above
scatter(x{4}, y{4}, 'o') % Essentially identical to the above
scatter(x{3}, y{3}, '*')
legend('Theorem 2.3, 2.6 & 2.16', 'Theorem 2.5')

xticks(1:length(T));
tostr = @(T) sprintf('(%i, %i)', T(1), T(2));
xticklabels(cellfun(tostr, T, 'UniformOutput', false))
xlim([0.5, length(T)+0.5])
xlabel('(T_{min}, T_{max})')
ylabel('\beta')

% Comparison with IQC and clock based approach: The latter approach seems 
% to involve some additional conservatism
figure
hold on; box on
scatter(x{1}, y{1}, 'o', 'LineWidth', 1); 
scatter(x{2}, y{2}, 'o') % Essentially identical to the above
scatter(x{4}, y{4}, 'o') % Essentially identical to the above
scatter(x{5}, y{5}, '*')
legend('Theorem 2.3, 2.6 & 2.16', 'Theorem 2.13')

xticks(1:length(T));
tostr = @(T) sprintf('(%i, %i)', T(1), T(2));
xticklabels(cellfun(tostr, T, 'UniformOutput', false))
xlim([0.5, length(T)+0.5])
xlabel('(T_{min}, T_{max})')
ylabel('\beta')

