%---------------------------------------------------------------------------------------------------
% For Paper
% "A Scalable Approach for Analysing Multi-Agent Systems with Heterogeneous Stochastic Packet Loss"
% by C. Hespe and H. Werner
% Copyright (c) Institute of Control Systems, Hamburg University of Technology. All rights reserved.
% Licensed under the GPLv3. See LICENSE in the project root for license information.
% Author(s): Christian Hespe
%---------------------------------------------------------------------------------------------------

function mustBeEqualSize(a,b)
%MUSTBEEQUALSIZE Verifies that a & b have equal size
%   Implementation of a argument verification function to constrain matrix and system arguments to
%   be of same size.
%   The implementation is based on:
%   https://uk.mathworks.com/help/matlab/matlab_prog/function-argument-validation-1.html

if ~isequal(size(a),size(b))
    eid = 'Size:notEqual';
    msg = 'Inputs must have equal size.';
    throwAsCaller(MException(eid,msg))
end
