%---------------------------------------------------------------------------------------------------
% For Paper
% "A Scalable Approach for Analysing Multi-Agent Systems with Heterogeneous Stochastic Packet Loss"
% by C. Hespe and H. Werner
% Copyright (c) Institute of Control Systems, Hamburg University of Technology. All rights reserved.
% Licensed under the GPLv3. See LICENSE in the project root for license information.
% Author(s): Christian Hespe
%---------------------------------------------------------------------------------------------------

function [H2, Y, solver_stats] = h2norm_decomposed(sysD, sysC, sysP, lambda, rho, options)
%H2NORM_DECOMPOSED Calculate an upper bound on the H2-norm of a decomposable jump system in a
%scalable manner
%   The LMIs in this function are formulated in terms of the observability Gramian.
%
%   Arguments:
%       sysD   -> Decoupled part of the system
%       sysC   -> Stochastically coupled part of the system
%       sysP   -> Deterministically coupled part of the system
%       lambda -> Eigenvalues of the nominal graph Laplacian of the system
%       p      -> Probability of a successful package transmission
%   Returns:
%       H2           -> Upper bound on the H2-norm of the system
%       Y            -> Storage function matrix of the solution
%       solver_stats -> Timing information about the algorithm

arguments
    sysD ss
    sysC ss {mustBeEqualSize(sysC, sysD), mustBeEqualOrder(sysC, sysD)}
    sysP ss {mustBeEqualSize(sysP, sysD), mustBeEqualOrder(sysP, sysD)}
    lambda (1,:) {mustBeNonnegative, mustBeFinite}
    rho (1,1) {mustBeNonnegative, mustBeLessThanOrEqual(rho,1)}

    options.offset (1,1) {mustBePositive, mustBeFinite} = 1e-8
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
N  = length(lambda);

lambda = sort(lambda);

%% Solve the SDP from Theorem 2
Y = sdpvar(nx);
Z = sdpvar(nw, nw, N-1);

cost  = 0;
cnstr = Y >= options.offset * eye(nx);

% The eigenvalues of L0 are sorted by Matlab. By iteration only over 2:N,
% we ignore lambda_1 = 0 and thus the uncontrollable and marginally stable
% modal subsystem.
for i = 2:N
    li   = lambda(i);
    lit  = rho*(1-rho)*li;
    Zi   = Z(:,:,i-1);
    cost = cost + trace(Zi);

    Acl = Ad + rho*li*Ac + li*Ap;
    Bcl = Bd + rho*li*Bc + li*Bp;
    Ccl = Cd + rho*li*Cc + li*Cp;
    Dcl = Dd + rho*li*Dc + li*Dp;

    LMI = Acl'*Y*Acl + Ccl'*Ccl - Y  + 2*lit * (Ac'*Y*Ac + Cc'*Cc);
    TRC = Bcl'*Y*Bcl + Dcl'*Dcl - Zi + 2*lit * (Bc'*Y*Bc + Dc'*Dc);

    % For certain LMIs that are obviously infeasible, Yalmip will refuse to
    % construct the SDP and issue an error. This try catch will convert
    % that error into a warning so that we can successfully finish running
    % the remainder of the script.
    try
        cnstr = [ cnstr                            ;
                  LMI <= -options.offset * eye(nx) ;
                  TRC <= -options.offset * eye(nw) ];
    catch ME
        warning(ME.message)
        Y = [];
        solver_stats.total = toc(timer);
        return
    end
end

%% Solve the SDP
solver_stats.prep = toc(timer);
sol = optimize(cnstr, cost, opts);

if sol.problem ~= 0
    if options.verbose > 0
        warning('YALMIP return an error: %s', sol.info)
    end
    Y = [];
else
    % We calculate gamma^2 with the LMI constraints, so we need to take the
    % square root here.
    H2 = sqrt(value(cost));

    Y  = value(Y);
    solver_stats.yalmip = sol.yalmiptime;
    solver_stats.solver = sol.solvertime;
end

solver_stats.total = toc(timer);
