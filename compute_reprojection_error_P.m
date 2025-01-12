function errors = compute_reprojection_error_P(P, X, x)
% COMPUTE_REPROJECTION_ERROR - Compute reprojection error for 3D points
% Inputs:
%   P - Projection matrix (3x4)
%   X - 3D points (4xN homogeneous)
%   x - 2D points (3xN homogeneous)
% Outputs:
%   errors - Reprojection errors (1xN)

    % Project 3D points to 2D
    x_proj = P * X;
    x_proj = x_proj ./ x_proj(3, :); % Normalize to homogeneous

    % Compute Euclidean distance error
    errors = sqrt(sum((x_proj(1:2, :) - x(1:2, :)).^2, 1)); % Pixel errors
end