function [X_refined, accumulated_error] = LevenbergMarquardt(P1, P2, X, x1, x2, mu, max_iterations, tolerance)
    % Levenberg-Marquardt optimization for 3D point refinement

    % Temporary variable for updates
    X_temp = X;
    accumulated_error = zeros(1, max_iterations); % For plotting error

    for iter = 1:max_iterations
        % disp(['Iteration Index: ', num2str(iter)]);

        % Compute total error before update
        total_error_before = 0;
        for j = 1:size(X_temp, 2)
            [err, ~] = ComputeReprojectionError(P1, P2, X_temp(:, j), x1(:, j), x2(:, j));
            total_error_before = total_error_before + err;
        end
        accumulated_error(iter) = total_error_before; % Store error for plotting
        % disp(['Total Error Before: ', num2str(total_error_before)]);

        % Update 3D points
        for j = 1:size(X_temp, 2)
            % Compute residual and Jacobian for point Xj
            [r, J] = LinearizeReprojErr(P1, P2, X_temp(:, j), x1(:, j), x2(:, j));

            % Compute LM update
            delta_Xj = ComputeUpdate(r, J, mu);

            % Update 3D point (temporarily)
            X_test = pflat(X_temp(:, j) + delta_Xj);

            % Compute individual errors
            [err_before, ~] = ComputeReprojectionError(P1, P2, X_temp(:, j), x1(:, j), x2(:, j));
            [err_after, ~] = ComputeReprojectionError(P1, P2, X_test, x1(:, j), x2(:, j));

            % Apply update if error improves
            if err_after < err_before
                X_temp(:, j) = X_test;
            end
        end

        % Compute total error after update
        total_error_after = 0;
        for j = 1:size(X_temp, 2)
            [err, ~] = ComputeReprojectionError(P1, P2, X_temp(:, j), x1(:, j), x2(:, j));
            total_error_after = total_error_after + err;
        end
        % disp(['Total Error After: ', num2str(total_error_after)]);

        % Adjust damping factor and check for convergence
        if total_error_after < total_error_before
            mu = mu / 2; % Decrease damping factor
        else
            mu = mu * 2; % Increase damping factor
        end

        if abs(total_error_after - total_error_before) < tolerance
            % disp('Converged!');
            break;
        end
    end

    % Return refined points and accumulated errors
    X_refined = X_temp;
end