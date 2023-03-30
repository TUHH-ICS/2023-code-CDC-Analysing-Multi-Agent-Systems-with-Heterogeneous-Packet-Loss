%---------------------------------------------------------------------------------------------------
% For Paper
% "A Scalable Approach for Analysing Multi-Agent Systems with Heterogeneous Stochastic Packet Loss"
% by C. Hespe and H. Werner
% Copyright (c) Institute of Control Systems, Hamburg University of Technology. All rights reserved.
% Licensed under the GPLv3. See LICENSE in the project root for license information.
% Author(s): Christian Hespe
%---------------------------------------------------------------------------------------------------

clear

addpath('analysis', 'graphs', 'util')

%% Settings
iter  = 10;     % Iterations for estimating the computation time
rho_l = 0.4;    % Lower probability bound
rho_u = 0.6;    % Upper probability bound

h = (2:2:142)'; % Base length of the triangle graph
N = h.*(h+1)/2;

%% Set up the problem
m  = 1;   % Model mass [kg]
b  = 0.9; % Coefficient of friction [kg/s]
dT = 1;   % Sampling time [s]

kappa = 0.01; % Loop gain

% Discretized state-space model of a mass with friction
A = [ 1     dT    ;
      0  1-dT*b/m ];
B = [ 0; dT/m ];
C = [ 1    0  ];

%% Prepare for system analysis

% Set up generalized plant
sysD = ss(A,          B,          zeros(1,2), 0, 1);
sysC = ss(-kappa*B*C, zeros(2,1), zeros(1,2), 0, 1);
sysP = ss(zeros(2),   zeros(2,1), C,          0, 1);

% Prepare evaluation grid
[H, ~] = ndgrid(h, 1:iter);

%% Evaluate performance for different uncertainty ranges
% Pre-allocate storage
H2_rbst = NaN(size(H));
H2_sing = NaN(size(H));
stats_rbst(length(h), iter) = struct('prep', NaN, 'solver', NaN, 'total', NaN, 'yalmip', NaN);
stats_sing(length(h), iter) = struct('prep', NaN, 'solver', NaN, 'total', NaN, 'yalmip', NaN);

% Create progress meter
meter = ParforProgressMeter(numel(H));
meter.start()

tic
parfor i = 1:numel(H)
    % Calculate Laplacian
    G0     = triangle_graph(H(i));
    L0     = full(laplacian(G0));
    lambda = eig(L0);

    % Improve numerical conditioning for very small eigenvalues
    lambda(lambda < eps) = 0;

    % Evaluate performance
    [H2_rbst(i), ~, stats_rbst(i)] = h2norm_decomposed_robust(sysD, sysC, sysP, lambda, [rho_l, rho_u]);
    [H2_sing(i), ~, stats_sing(i)] = h2norm_decomposed_robust(sysD, sysC, sysP, lambda, (rho_l+rho_u)/2);

    % Update progress meter
    meter.notify(i);
end
disp(['Scalability test completed in ' format_duration(toc) ' on ' char(datetime())])

% Calculate average timings
time_rbst = zeros(length(h),1);
time_sing = zeros(length(h),1);
for i = 1:length(h)
    time_rbst(i) = mean([stats_rbst(i,:).total]);
    time_sing(i) = mean([stats_sing(i,:).total]);
end

H2_rbst = H2_rbst(:,1);
H2_sing = H2_sing(:,1);

%% Visualize
figure()
yyaxis left
loglog(N, H2_rbst, N, H2_sing, '--')
xlabel("Agent Count")
ylabel("H2-Performance")

yyaxis right
loglog(N, time_rbst, N, time_sing, '--')
ylabel("Computation time")

%% CSV Export
data = table(h, N, H2_rbst, H2_sing, time_rbst, time_sing);
name = sprintf('scaling_%d.csv', uint32(posixtime(datetime())));
writetable(data, name)
