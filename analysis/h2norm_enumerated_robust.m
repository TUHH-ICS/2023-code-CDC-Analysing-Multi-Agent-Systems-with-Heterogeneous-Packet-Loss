%---------------------------------------------------------------------------------------------------
% For Paper
% "A Scalable Approach for Analysing Multi-Agent Systems with Heterogeneous Stochastic Packet Loss"
% by C. Hespe and H. Werner
% Copyright (c) Institute of Control Systems, Hamburg University of Technology. All rights reserved.
% Licensed under the GPLv3. See LICENSE in the project root for license information.
% Author(s): Christian Hespe
%---------------------------------------------------------------------------------------------------

function [H2, X, solver_stats] = h2norm_enumerated_robust(sysD, sysC, sysP, L0, rho, options)
%H2NORM_ENUMERATED_ROBUST Calculates the H2-norm of the given decomposable jump system with
%Bernoulli packet loss.
%   This function calculates the H2-norm of the given decomposable jump system by application of
%   Theorem 2 from the accompanying paper. To do so, it needs to enumerate all possible modes, the
%   number of which grows exponentially with the number of agents in the system.
%
%   Arguments:
%       sysD      -> Decoupled part of the plant
%       sysC      -> Coupled part of the plant
%       sysP      -> Deterministic interconnected part of the plant
%       L0        -> Laplacian describing the nominal graph of the system
%       rho       -> Uncertainty interval for successfull package transmission
%       symmetric -> [Default=false] Handle symmetric or asymmetric loss
%   Returns:
%       H2           -> H2-norm of the system
%       X            -> Storage function matrix of the solution
%       solver_stats -> Timing information about the algorithm

arguments
    sysD ss
    sysC ss {mustBeEqualSize(sysC, sysD), mustBeEqualOrder(sysC, sysD)}
    sysP ss {mustBeEqualSize(sysP, sysD), mustBeEqualOrder(sysP, sysD)}
    L0 (:,:) {mustBeNumeric, mustBeFinite}
    rho (1,2) {mustBeNonnegative, mustBeLessThanOrEqual(rho,1)}

    options.offset (1,1) {mustBePositive, mustBeFinite} = 1e-8
    options.points (1,1) {mustBeInteger, mustBePositive} = 10;
    options.symmetric (1,1) logical = false
    options.verbose (1,1) {mustBeMember(options.verbose, [0,1,2])} = 0
end

% Setup the SDP solver
opts = sdpsettings('verbose', options.verbose);

% Allocate storage for time measurements
solver_stats        = struct;
solver_stats.prep   = NaN;
solver_stats.solver = NaN;
solver_stats.total  = NaN;
solver_stats.yalmip = NaN;

% Initialize return values with reasonable defaults
H2    = inf;
timer = tic;

%% Prepare the system description
[Ad, Bd, Cd, Dd] = ssdata(sysD);
[Ac, Bc, Cc, Dc] = ssdata(sysC);
[Ap, Bp, Cp, Dp] = ssdata(sysP);

nx = size(Ad,1);
nw = size(Bd,2);
N  = size(L0,1);

% Calculate adjacency matrix and truncated diagonalizing transformation
A0 = diag(diag(L0)) - L0;
[T, ~] = eig(L0);
T  = T(:,2:end);
Lr = T'*L0*T;

%% Split network into separete channels
if options.symmetric
    A0 = tril(A0);
end

channels = find(A0);
m = length(channels);

%% Prepare the SDP from Theorem 2
% To generate the required SDP, we need to enumerate all modes of the jump
% system. This is done by splitting the communication into channels and
% disableing them individually.

X = sdpvar(nx*(N-1));
Z = sdpvar(nw*(N-1));
cnstr = X >= options.offset * eye(nx*(N-1));

% Find all edges of the graph
edges = find(tril(A0));
pvec  = linspace(rho(1), rho(2), options.points);

for j = 1:options.points^length(edges)
    % Calculate corresponding probability matrix
    P = zeros(size(A0));
    digits = j - 1;
    for edge = edges'
        P(edge) = pvec(mod(digits, options.points) + 1);
        digits  = floor(digits / options.points);
    end
    P = P + P';

    % Initialize LMIs
    LMI = -X;
    TRC = -Z;

    for i = 1:2^m
        % Calculate new adjacency matrix with certain channels disabled
        A       = A0;
        mask    = channels(dec2bin(i-1, m) ~= '0');
        A(mask) = 0;

        % Calculate the probability for this case
        chance = prod(P(A == 1)) * prod(1 - P(A0 ~= A));

        % Calculate current Laplacian
        if options.symmetric
            A = A + A';
        end
        L = T'*(diag(sum(A, 2)) - A)*T;

        % Assemble system for that loss pattern
        Acl = kron(eye(N-1), Ad) + kron(L, Ac) + kron(Lr, Ap);
        Bcl = kron(eye(N-1), Bd) + kron(L, Bc) + kron(Lr, Bp);
        Ccl = kron(eye(N-1), Cd) + kron(L, Cc) + kron(Lr, Cp);
        Dcl = kron(eye(N-1), Dd) + kron(L, Dc) + kron(Lr, Dp);

        % Assemble the LMI
        LMI = LMI + chance * (Acl'*X*Acl + Ccl'*Ccl);
        TRC = TRC + chance * (Bcl'*X*Bcl + Dcl'*Dcl);
    end

    % For certain LMIs that are obviously infeasible, Yalmip will refuse to
    % construct the SDP and issue an error. This try catch will convert that
    % error into a warning so that we can successfully finish running the
    % remainder of the script.
    try
        cnstr = [ cnstr                                   ;
                  LMI <= -options.offset * eye(nx*(N-1)) ;
                  TRC <= -options.offset * eye(nw*(N-1)) ];
    catch ME
        warning(ME.message)
        X = [];
        solver_stats.total = toc(timer);
        return
    end
end

%% Solve the SDP
solver_stats.prep = toc(timer);
sol = optimize(cnstr, trace(Z), opts);

if sol.problem ~= 0
    if options.verbose > 0
        warning('YALMIP return an error: %s', sol.info)
    end
    X = [];
else
    % We calculate gamma^2 with the LMI constraints, so we need to take the
    % square root here.
    H2 = sqrt(value(trace(Z)));

    X = value(X);
    solver_stats.yalmip = sol.yalmiptime;
    solver_stats.solver = sol.solvertime;
end

solver_stats.total = toc(timer);
