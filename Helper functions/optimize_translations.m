function [T1_refined, T2_refined, accumulated_error] = optimize_translations(K, R1, R2, T1, T2, X, x1, x2, mu, max_iterations, tolerance)
    % OPTIMIZE_TRANSLATIONS - Refines translation vectors (T1, T2) using 
    % Levenberg-Marquardt optimization.

    % Initialize variables
    T1_temp = T1; % Translation for camera 1
    T2_temp = T2; % Translation for camera 2
    accumulated_error = zeros(1, max_iterations); % Error log for plotting

    % Optimization loop
    for iter = 1:max_iterations
        % Compute Total Error Before Update
        total_error_before = 0;
        for j = 1:size(X, 2)
            err = ComputeReprojectionError(K, R1, R2, T1_temp, T2_temp, X(:, j), x1(:, j), x2(:, j));
            total_error_before = total_error_before + err;
        end
        accumulated_error(iter) = total_error_before; % Store error

        % Compute Jacobian and Residual
        [r, J] = ComputeJacobianAndResidual(K, R1, R2, T1_temp, T2_temp, X, x1, x2);

        % Compute LM Update
        delta_T = ComputeUpdate(r, J, mu); % Update for T1 and T2

        % Apply Updates
        T1_test = T1_temp + delta_T(1:3); % Update T1
        T2_test = T2_temp + delta_T(4:6); % Update T2

        % Compute Total Error After Update
        total_error_after = 0;
        for j = 1:size(X, 2)
            err = ComputeReprojectionError(K, R1, R2, T1_test, T2_test, X(:, j), x1(:, j), x2(:, j));
            total_error_after = total_error_after + err;
        end

        % Accept or Reject Update
        if total_error_after < total_error_before
            T1_temp = T1_test;
            T2_temp = T2_test;
            mu = mu / 2; % Decrease damping
        else
            mu = mu * 2; % Increase damping
        end

        % Check Convergence
        if abs(total_error_after - total_error_before) < tolerance
            break; % Converged
        end
    end

    % Return refined translations
    T1_refined = T1_temp;
    T2_refined = T2_temp;
end
%% ----------------------------------------------------------
% Helper Functions
%% ----------------------------------------------------------

% Compute Reprojection Error
function error = ComputeReprojectionError(K, R1, R2, T1, T2, X, x1, x2)
    % Projection matrices
    P1 = K * [R1, T1]; % Camera 1
    P2 = K * [R2, T2]; % Camera 2

    % Ensure homogeneous coordinates
    if size(X, 1) == 3
        X = [X; 1];
    end

    % Project points
    x1_proj = P1 * X;
    x2_proj = P2 * X;

    % Normalize projections
    x1_proj = x1_proj(1:2) / x1_proj(3);
    x2_proj = x2_proj(1:2) / x2_proj(3);

    % Compute error
    err1 = norm(x1_proj - x1(1:2));
    err2 = norm(x2_proj - x2(1:2));
    error = err1 + err2;
end

% Compute Jacobian and Residuals
function [r, J] = ComputeJacobianAndResidual(K, R1, R2, T1, T2, X, x1, x2)
    % Ensure homogeneous coordinates
    if size(X, 1) == 3
        X = [X; 1];
    end

    % Compute projections
    P1 = K * [R1, T1];
    P2 = K * [R2, T2];
    x1_proj = P1 * X;
    x2_proj = P2 * X;
    x1_proj = x1_proj(1:2) / x1_proj(3);
    x2_proj = x2_proj(1:2) / x2_proj(3);

    % Compute residuals
    r = [x1_proj - x1(1:2); x2_proj - x2(1:2)];

    % Compute Jacobian for T1 and T2
    J1 = K * R1; % Partial derivatives w.r.t T1
    J2 = K * R2; % Partial derivatives w.r.t T2

    % Final Jacobian (4x6)
    J = [J1;J2];
end

function delta = ComputeUpdate(r, J, mu)
    % Hessian and update computation
    H = J' * J + mu * eye(size(J, 2));
    delta = -H \ (J' * r);
end