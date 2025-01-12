function [E_best, indices] = estimate_E_robust(x1_h_n, x2_h_n, epipolar_threshold)
    % Robustly estimates the essential matrix E using RANSAC.

    % x1: 3xN , N: Number of points
    % x2: 3XN , N: Number of points
    
    global enableInfo;


    % RANSAC parameters
    max_iterations = 10000;
    num_points = size(x1_h_n, 2);

    best_inlier_count = 0;
    E_best = [];
    inliers = false(1, num_points);

    % RANSAC loop
    for iter = 1:max_iterations
        % Randomly sample 8 points correspondences
        sample_indices = randperm(num_points, 8);

        x1_h_n_sample = x1_h_n(:, sample_indices);
        x2_h_n_sample = x2_h_n(:, sample_indices);

        % Estimate Essential matrix
        E_candidate = estimate_F_DLT(x1_h_n_sample, x2_h_n_sample); 
        E_candidate = enforce_essential(E_candidate);
        
        %Just for validation
        validate_E(E_candidate)

        % Compute epipolar errors for all points
        distances_l2_x2 = compute_epipolar_errors(E_candidate, x1_h_n, x2_h_n);
        distances_l1_x1 = compute_epipolar_errors(E_candidate', x2_h_n, x1_h_n);

        % Symmetric epipolar error
        errors = (distances_l2_x2.^2 + distances_l1_x1.^2) / 2;

        % Determine inliers
        inliers_current = errors < epipolar_threshold^2;
        num_inliers = sum(inliers_current);

        % Update best model
        if num_inliers > best_inlier_count
            best_inlier_count = num_inliers;
            E_best = E_candidate;
            inliers = inliers_current;
        end
    end

    % fprintf('Best inlier count: %d\n', best_inlier_count);
    % fprintf('Percentage of inliers: %.2f%%\n', 100 * (best_inlier_count) / num_points);

    info("Percentage of inliers: %.2f \n:  ", 100 * (best_inlier_count) / num_points);


    % Return the indices of the inliners
    indices = find(inliers);
end