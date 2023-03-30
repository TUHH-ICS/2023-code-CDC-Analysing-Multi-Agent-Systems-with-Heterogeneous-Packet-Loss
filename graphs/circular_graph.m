%---------------------------------------------------------------------------------------------------
% For Paper
% "A Scalable Approach for Analysing Multi-Agent Systems with Heterogeneous Stochastic Packet Loss"
% by C. Hespe and H. Werner
% Copyright (c) Institute of Control Systems, Hamburg University of Technology. All rights reserved.
% Licensed under the GPLv3. See LICENSE in the project root for license information.
% Author(s): Christian Hespe
%---------------------------------------------------------------------------------------------------

function G = circular_graph(nvert, options)
%CIRCULAR_GRAPH Generates a Matlab graph object that represents a circular graph with a given number
%of vertices.
%   This function genenerates a Matlab graph object that contains a circular graph with the given
%   amount of vertices. It can be configured if the connections should be directed or not and to how
%   many of the next vertices the connection should be established.
%
%   Arguments:
%       nvert    -> Number of vertices
%       nedge    -> [Default=1] Number of forward edges per vertex
%       directed -> [Default=false] Directed edges or not

arguments
    nvert (1,1) {mustBeInteger, mustBePositive}
    options.nedge (1,1) {mustBeInteger, mustBePositive, mustBeLessThanOrEqual(options.nedge, nvert)} = 1
    options.directed (1,1) logical = false
end

% Define closure that keeps the index in the ring [1, nvert]
ring = @(i) mod(i-1, nvert) + 1;

%% Generate graph structure
A = zeros(nvert);
for i = 1:nvert
    A(i, ring(i+(1:options.nedge))) = 1;

    % If not direction, also add the other direction
    if ~options.directed
        A(i, ring(i-(1:options.nedge))) = 1;
    end
end

% G needs only to be a digraph if it is directed
if options.directed
    G = digraph(A);
else
    G = graph(A);
end
