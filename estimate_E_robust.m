function [E_best, indices] = estimate_E_robust(x1, x2, eps, K)
    % Robustly estimates the essential matrix E using RANSAC.
    
    % Normalize points
    x1_normalized = inv(K) * x1;
    x2_normalized = inv(K) * x2;

    % RANSAC parameters
    max_iterations = 10000;
    num_points = size(x1, 2);
    best_inlier_count = 0;
    E_best = [];
    inliers = false(1, num_points);

    % RANSAC loop
    for iter = 1:max_iterations
        % Randomly sample 8 points correspondences
        sample_indices = randperm(num_points, 8);
        x1_sample = x1_normalized(:, sample_indices);
        x2_sample = x2_normalized(:, sample_indices);

        % Estimate essential matrix
        E_candidate = estimate_F_DLT(x1_sample, x2_sample);
        E_candidate = enforce_essential(E_candidate);

        % Compute epipolar errors for all points
        distances_l2_x2 = compute_epipolar_errors(E_candidate, x1_normalized, x2_normalized);
        distances_l1_x1 = compute_epipolar_errors(E_candidate', x2_normalized, x1_normalized);

        % Symmetric epipolar error
        errors = (distances_l2_x2.^2 + distances_l1_x1.^2) / 2;

        % Determine inliers
        inliers_current = errors < eps^2;
        num_inliers = sum(inliers_current);

        % Update best model
        if num_inliers > best_inlier_count
            best_inlier_count = num_inliers;
            E_best = E_candidate;
            inliers = inliers_current;
        end
    end

    fprintf('Best inlier count: %d\n', best_inlier_count);
    indices = find(inliers);
end