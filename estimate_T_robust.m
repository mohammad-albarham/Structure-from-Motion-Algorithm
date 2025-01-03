% Function: Robust Translation Estimation
function T_best = estimate_T_robust(xs, Xs, R, inlier_threshold, T_init)

    % Initialization
    max_inliers = 0;
    N = size(xs, 2);
    iterations = 5000; % RANSAC iterations

    for i = 1:iterations
        % Sample minimal points (2 points for translation)
        sample_indices = randperm(N, 2);
        xs_sample = xs(:, sample_indices);
        Xs_sample = Xs(:, sample_indices);

        % Estimate candidate T using DLT
        T_candidate = estimate_translation_DLT(xs_sample, Xs_sample, R);

        % Compute reprojection errors
        errors = compute_reprojection_errors(xs, Xs, R, T_candidate);

        % Count inliers
        inliers_candidate = errors.^2 < inlier_threshold.^2;
        num_inliers = sum(inliers_candidate);

        % Update best model
        if num_inliers > max_inliers
            max_inliers = num_inliers;
            T_best = T_candidate;
            inliers = inliers_candidate;
        end
    end

    % Refine with inliers
    if max_inliers > 0
        T_best = T_best; % estimate_translation_DLT(xs(:, inliers), Xs(:, inliers), R);
    else
        warning('No inliers found. Returning initial translation.');
        T_best = T_init;
    end
end


% Function: Translation Estimation using DLT
function T = estimate_translation_DLT(x, X, R)
    % Inputs:
    % x - 2D normalized points in homogeneous coordinates (3xN)
    % X - 3D points in homogeneous coordinates (4xN)
    % R - Rotation matrix (3x3)

    % Outputs:
    % T - Estimated translation vector (3x1)

    % Number of points
    N = size(x, 2);

    % Initialize matrices A and b
    A = [];
    b = [];

    % Loop through each point
    for i = 1:N
        % Extract the 2D point (x) and 3D point (X)
        xi = x(:, i);       % 2D point in homogeneous coordinates
        Xi = X(:, i);       % 3D point in homogeneous coordinates

        % Create the linear system based on the equations in the image
        A_i = [
            Xi(4), 0, 0, -xi(1);   % First row
            0, Xi(4), 0, -xi(2);   % Second row
            0, 0, Xi(4), -xi(3)    % Third row
        ];

        % Append to the full matrix A
        A = [A; A_i];

        % Compute b using the rotation matrix R
        b_i = -R * Xi(1:3);  % Projected 3D point
        b = [b; b_i];
    end

    % Solve the linear system A * T = b using least squares
    T_lambda = A \ b;

    % Extract translation vector T (first 3 elements)
    T = T_lambda(1:3);

    % Compute residual error
    % residual = norm(A * T_lambda - b); % ||A * T_lambda - b||
    % disp('Residual Error:');
    % disp(residual);

    % Flip sign if necessary based on direction consistency
    if dot(T, mean(X(1:3, :), 2)) < 0
        T = -T;  % Flip the translation vector
    end
end


% Function: Compute Reprojection Errors
function errors = compute_reprojection_errors(xs, Xs, R, T)
    % Project 3D points
    projected = (R * Xs(1:3, :) + T);
    projected = projected(1:2, :) ./ projected(3, :); % Normalize

    % Compute errors
    errors = sqrt(sum((xs(1:2, :) - projected).^2, 1)); % Euclidean distance
end
