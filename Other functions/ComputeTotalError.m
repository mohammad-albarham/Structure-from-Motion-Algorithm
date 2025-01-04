function total_error = ComputeTotalError(P1, P2, X, x1, x2)
    % Compute the total reprojection error for all 3D points
    %
    % Inputs:
    %   - P1, P2: Camera projection matrices
    %   - X: 4xN matrix of 3D points in homogeneous coordinates
    %   - x1, x2: 3xN matrices of 2D image points in homogeneous coordinates
    %
    % Outputs:
    %   - total_error: The total reprojection error across all points

    total_error = 0; % Initialize total error
    for j = 1:size(X, 2)
        [err, ~] = ComputeReprojectionError(P1, P2, X(:, j), x1(:, j), x2(:, j));
        total_error = total_error + err;
    end
end