% Function: Compute Reprojection Errors
function errors = compute_reprojection_errors(xs, Xs, R, T)
    % Project 3D points
    projected = (R * Xs(1:3, :) + T);
    projected = projected(1:2, :) ./ projected(3, :); % Normalize

    % Compute errors
    errors = sqrt(sum((xs(1:2, :) - projected).^2, 1)); % Euclidean distance
end
