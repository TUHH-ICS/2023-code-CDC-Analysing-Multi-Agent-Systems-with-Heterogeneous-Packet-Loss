%---------------------------------------------------------------------------------------------------
% For Paper
% "A Scalable Approach for Analysing Multi-Agent Systems with Heterogeneous Stochastic Packet Loss"
% by C. Hespe and H. Werner
% Copyright (c) Institute of Control Systems, Hamburg University of Technology. All rights reserved.
% Licensed under the GPLv3. See LICENSE in the project root for license information.
% Author(s): Christian Hespe
%---------------------------------------------------------------------------------------------------

function mustBeEqualOrder(a,b)
%MUSTBEEQUALORDER Verifies that a & b have equal (dynamic) order
%   Implementation of a argument verification function to verify that two dynamic systems are of
%   equal dynamic order.

if ~isequal(order(a),order(b))
    eid = 'Order:notEqual';
    msg = 'Inputs must have equal order.';
    throwAsCaller(MException(eid,msg))
end
