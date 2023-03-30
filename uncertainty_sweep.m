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
p_swp = [0.002:0.001:0.01, 0.015:0.005:0.1, 0.11:0.01:0.4, 0.42:0.02:1]';
N_swp = [4,6];

%% Set up the problem
m  = 1;   % Model mass [kg]
b  = 0.9; % Coefficient of friction [kg/s]
dT = 1;   % Sampling time [s]

kappa = 0.15; % Loop gain

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
[P, N] = ndgrid(p_swp, N_swp);

%% Evaluate performance for different uncertainty ranges
H2_dec = NaN(size(P));
H2_enm = NaN(size(P));

% Create progress meter
meter = ParforProgressMeter(numel(P), interval=0.01);
meter.start()

tic
parfor i = 1:numel(P)
    % Calculate Laplacian
    G0     = circular_graph(N(i));
    L0     = full(laplacian(G0));
    lambda = eig(L0);

    % Improve numerical conditioning for very small eigenvalues
    lambda(lambda < eps) = 0;

    % Evaluate performance
    H2_dec(i) = h2norm_decomposed_robust(sysD, sysC, sysP, lambda, [P(i), 1]);
    H2_enm(i) = h2norm_enumerated_robust(sysD, sysC, sysP, L0,     [P(i), 1], points=2);

    % Update progress meter
    meter.notify(i);
end
disp(['Uncertainty sweep completed in ' format_duration(toc) ' on ' char(datetime())])

%% Visualize
figure()
plot(p_swp, H2_dec(:,1), p_swp, H2_enm(:,1))
hold on
plot(p_swp, H2_dec(:,2), '--', p_swp, H2_enm(:,2), '--')
hold off
xlabel("Loss Probability")
ylabel("H2-Performance")

ylim padded
legend('Decomposed', 'Enumerated')

%% Export results
name = sprintf('uncertainty_sweep_%d.csv', uint32(posixtime(datetime())));
tbl  = table(p_swp, H2_dec, H2_enm);
writetable(tbl, name)
