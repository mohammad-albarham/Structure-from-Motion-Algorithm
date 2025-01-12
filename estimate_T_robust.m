% Function: Robust Translation Estimation
function T_best = estimate_T_robust(xs, Xs, R, inlier_threshold)


    % Initialization
    max_inliers = 0;
    N = size(xs, 2);
    iterations = 10000; % RANSAC iterations

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
        inliers_candidate = errors < inlier_threshold;
        num_inliers = sum(inliers_candidate);

        % Update best model
        if num_inliers > max_inliers
            max_inliers = num_inliers;
            T_best = T_candidate;
            inliers = inliers_candidate;
        end
    end


    % Refine with inliers
    if max_inliers == 0
        warning('No inliers found. Returning initial translation.');
        T_best = zeros(3,1);
    end

    errors = compute_reprojection_errors(xs, Xs, R, T_best);

    % Count inliers
    inliers_candidate = errors < inlier_threshold;
    num_inliers = sum(inliers_candidate);

    info('Percentage of inliers in T_robust: %.2f\n',2,100 * (num_inliers) / N);
end
