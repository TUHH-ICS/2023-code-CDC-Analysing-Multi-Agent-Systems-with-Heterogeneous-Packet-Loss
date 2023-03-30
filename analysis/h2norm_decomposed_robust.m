%---------------------------------------------------------------------------------------------------
% For Paper
% "A Scalable Approach for Analysing Multi-Agent Systems with Heterogeneous Stochastic Packet Loss"
% by C. Hespe and H. Werner
% Copyright (c) Institute of Control Systems, Hamburg University of Technology. All rights reserved.
% Licensed under the GPLv3. See LICENSE in the project root for license information.
% Author(s): Christian Hespe
%---------------------------------------------------------------------------------------------------

function [H2, Y, solver_stats] = h2norm_decomposed_robust(sysD, sysC, sysP, lambda, rho, options)
%H2NORM_DECOMPOSED_ROBUST A function to compute an upper bound on the H2-norm of an MAS with
%heterogeneous stochastic packet loss in a scalable manner.
%   The norm is evaluated in terms of the observability Gramian by applying the full block
%   S-procedure. The transition rates of the Markov chain are considered uncertain but lying inside
%   the interval specified by rho.
%
%   Arguments:
%       sysD   -> Decoupled part of the system
%       sysC   -> Stochastically coupled part of the system
%       sysP   -> Deterministically coupled part of the system
%       lambda -> Eigenvalues of the nominal graph Laplacian of the system
%       rho    -> Uncertainty interval for successfull package transmission
%   Returns:
%       H2           -> Upper bound on the H2-norm of the system
%       Y            -> Storage function matrix of the solution
%       solver_stats -> Timing information about the algorithm

arguments
    sysD ss
    sysC ss {mustBeEqualSize(sysC, sysD), mustBeEqualOrder(sysC, sysD)}
    sysP ss {mustBeEqualSize(sysP, sysD), mustBeEqualOrder(sysP, sysD)}
    lambda (1,:) {mustBeNonnegative, mustBeFinite}
    rho (1,2) {mustBeNonnegative, mustBeLessThanOrEqual(rho,1)}

    options.offset (1,1) {mustBePositive, mustBeFinite} = 1e-8
    options.points (1,1) {mustBeInteger, mustBePositive} = 100;
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
nz = size(Cd,1);
N  = length(lambda);

lambda = sort(lambda);

%% Formulate the SDP from Theorem 6
Y  = sdpvar(nx);
Z  = sdpvar(nw, nw, N-1);
P1 = sdpvar(5*(nx+nz));
P2 = sdpvar(5*(nx+nz));

cost  = 0;
cnstr = Y >= options.offset * eye(nx);

% The eigenvalues of L0 are sorted by Matlab. By iteration only over 2:N,
% we ignore lambda_1 = 0 and thus the uncontrollable and marginally stable
% modal subsystem.
for i = 2:N
    li   = lambda(i);
    Zi   = Z(:,:,i-1);
    cost = cost + trace(Zi);

    M1 = [ eye(nx)       zeros(nx, 3*(nx+nz))                                           ;
           Ad + li*Ap    sqrt(li)*eye(nx)      zeros(nx,2*nx+3*nz)                      ;
           Cd + li*Cp    zeros(nz,nx)          sqrt(li)*eye(nz)     zeros(nz,2*(nx+nz)) ;
           zeros(nx)     zeros(nx,nx+nz)       eye(nx)              zeros(nx,nx+2*nz)   ;
           zeros(nz,nx)  zeros(nz,2*nx+nz)     eye(nz)              zeros(nz,nx+nz)     ];
    M2 = [ eye(nw)       zeros(nw, 3*(nx+nz))                                           ;
           Bd + li*Bp    sqrt(li)*eye(nx)      zeros(nx,2*nx+3*nz)                      ;
           Dd + li*Dp    zeros(nz,nx)          sqrt(li)*eye(nz)     zeros(nz,2*(nx+nz)) ;
           zeros(nx,nw)  zeros(nx,nx+nz)       eye(nx)              zeros(nx,nx+2*nz)   ;
           zeros(nz,nw)  zeros(nz,2*nx+nz)     eye(nz)              zeros(nz,nx+nz)     ];
    K1 = [ zeros(3*(nx+nz),nx)  eye(3*(nx+nz))                             ;
           zeros(nx)            zeros(nx,2*(nx+nz))  eye(nx)  zeros(nx,nz) ;
           zeros(nz,nx)         zeros(nz,3*nx+2*nz)  eye(nz)               ;
           sqrt(li)*Ac          zeros(nx,3*(nx+nz))                        ;
           sqrt(li)*Cc          zeros(nz,3*(nx+nz))                        ];
    K2 = [ zeros(3*(nx+nz),nw)  eye(3*(nx+nz))                             ;
           zeros(nx,nw)         zeros(nx,2*(nx+nz))  eye(nx)  zeros(nx,nz) ;
           zeros(nz,nw)         zeros(nz,3*nx+2*nz)  eye(nz)               ;
           sqrt(li)*Bc          zeros(nx,3*(nx+nz))                        ;
           sqrt(li)*Dc          zeros(nz,3*(nx+nz))                        ];
    Psi_1 = blkdiag(-Y,  Y, eye(nz), 2*Y, 2*eye(nz));
    Psi_2 = blkdiag(-Zi, Y, eye(nz), 2*Y, 2*eye(nz));

    % For certain LMIs that are obviously infeasible, Yalmip will refuse to
    % construct the SDP and issue an error. This try catch will convert
    % that error into a warning so that we can successfully finish running
    % the remainder of the script.
    try
        cnstr = [ cnstr                                                           ;
                  M1'*Psi_1*M1 + K1'*P1*K1 <= -options.offset * eye(4*nx+3*nz)    ;
                  M2'*Psi_2*M2 + K2'*P2*K2 <= -options.offset * eye(3*(nx+nz)+nw) ];
    catch ME
        warning(ME.message)
        Y = [];
        solver_stats.total = toc(timer);
        return
    end
end

% Grid the uncertainty range
for p = linspace(rho(1), rho(2), options.points)
    Delta = [ sqrt(p)    0       ;
              sqrt(1-p)  0       ;
              0          sqrt(p) ];

    M     = [kron(Delta, eye(nx+nz)); eye(2*(nx+nz))];
    cnstr = [ cnstr                                      ;
              M'*P1*M >= options.offset * eye(2*(nx+nz)) ;
              M'*P2*M >= options.offset * eye(2*(nx+nz)) ];
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
