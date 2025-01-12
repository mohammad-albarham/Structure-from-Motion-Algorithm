function v_h = toHomogeneous(v)
% TOHOMOGENEOUS - Converts a vector to homogeneous coordinates.
%
% Usage:
%   v_h = toHomogeneous(v)
%
% Inputs:
%   v - Input vector (NxM) where N is the dimensionality, and M is the number of points.
%
% Outputs:
%   v_h - Output vector in homogeneous coordinates ((N+1)xM).
%
% Example:
%   v = [1; 2; 3]; % 3D point
%   v_h = toHomogeneous(v); % Convert to homogeneous -> [1; 2; 3; 1]

% Validate input
if nargin < 1
    error('Input vector is required.');
end

% Add a row of ones to make it homogeneous
v_h = [v; ones(1, size(v, 2))];

end
