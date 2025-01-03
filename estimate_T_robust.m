%%% Previous version 

%%
% Function: Robust Translation Estimation
function T_best = estimate_T_robust(xs, Xs, R, inlier_threshold,T_init, K)

    % Add a dimension to X to make it homogenous 
    % Fix dimensions
    xs = xs';
    Xs = Xs';  

    % Initialize variables
    max_inliers = 0; 
    N = size(xs, 1); % Updated for Nx2
    iterations = 10000; 
    inliers = false(1, N);

    for i = 1:iterations
        % Random sample (2 points) to solve this problem "Proof on the
        % paper"
        sample_indices = randperm(N, 2);
        xs_sample = xs(sample_indices, :);
        Xs_sample = Xs(sample_indices, :); %TODO: Solve the problem here 

        % Estimate candidate translation
        T_candidate = estimate_translation_DLT(xs_sample, Xs_sample, R);
        
        
        % Compute reprojection errors
        errors = compute_reprojection_errors(xs, Xs, R, T_candidate);
    
        errors = errors / K(1,1);

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
    % Final refinement with inliers
    if max_inliers > 0
        T_best = estimate_translation_DLT(xs(inliers, :), Xs(inliers, :), R);
    else
        warning('No inliers found. Returning initial translation.');
        T_best = T_init; % Fallback translation
    end
end


% Function: Translation Estimation using DLT
function T = estimate_translation_DLT(xs, Xs, R)
    xs = xs';
    Xs = Xs';
    N = size(xs, 2); % Number of points
    M = zeros(3*N,3+N); % Initialize the matrix 
    b = zeros(3*N, 1);

    for i = 1:N
        % First row in the block
        M(3*i-2, 1) = Xs(4,i)';
        M(3*i-2, 3+i) = -xs(1,i);

        % Second row in the block
        M(3*i-1, 2) = Xs(4,i)';
        M(3*i-1, 3+i) = -xs(2,i);

        % Third row in the block
        M(3*i, 3) = Xs(4,i)';
        M(3*i, 3+i) = -xs(3,i);

        % the b part in the block 

        b(3*(i-1)+1:3*i, 1) = - R *  Xs(1:3,i);
    end


    % Solve Using Pseudo-Inverse (if b ≠ 0)
    lambda = 1e-6; % Regularization parameter

    x_pinv = (M' * M + lambda * eye(size(M, 2))) \ (M' * (-b));

    T = x_pinv(1:3);

    % Fix sign ambiguity
    if dot(T, mean(Xs(1:3, :), 2)) < 0
        T = -T;
    end
    

end

% Function: Compute Reprojection Errors
function errors = compute_reprojection_errors(xs, Xs, R, T)
    projected = (R * Xs(:,1:3)' + T)'; 
    projected = projected(:, 1:2) ./ projected(:, 3); % Normalize to 2D
    errors = sqrt(sum((xs(1:2) - projected).^2, 2)); % Compute Euclidean distance
end


% % Function: Robust Translation Estimation
% function T_best = estimate_T_robust(x1, x2, Xs, R, translation_threshold,T_init)
% 
%     % Add a dimension to X to make it homogenous 
% 
%     % RANSAC parameters
%     max_iterations = 1000;
%     num_points = size(x1, 2);
%     best_inlier_count = 0;
%     E_best = [];
%     inliers = false(1, num_points);
% 
% 
%     for iter = 1:max_iterations
%         % Random sample (2 points) to solve this problem "Proof on the
%         % paper"
%         sample_indices = randperm(num_points, 2);
%         x1_sample = x1_normalized(:, sample_indices);
%         x2_sample = x2_normalized(:, sample_indices);
%         Xs_sample = Xs(sample_indices, :);
% 
% 
%         % Estimate candidate translation
%         T_candidate1 = estimate_translation_DLT(x1_sample, Xs_sample, R);
% 
%                 % Estimate candidate translation
%         T_candidate2 = estimate_translation_DLT(x1_sample, Xs_sample, R);
% 
%         % Compute reprojection errors
%         projc_error_1 = compute_reprojection_errors(x1_normalized, Xs_sample, R, T_candidate1);
% 
%         % Compute reprojection errors
%         projc_error_2 = compute_reprojection_errors(x1_normalized, Xs_sample, R, T_candidate2);
% 
% 
%         % Symmetric epipolar error
%         errors = (projc_error_1.^2 + projc_error_2.^2) / 2;
% 
% 
%         % Count inliers
%         inliers_candidate = errors < translation_threshold^2;
%         num_inliers = sum(inliers_candidate);
% 
%         % Update best model
%         if num_inliers > max_inliers
%             max_inliers = num_inliers;
%             T_best = T_candidate1;
%             inliers = inliers_candidate;
%         end
%     end
% 
% 
%     % Refine with inliers
%     % Final refinement with inliers
%     if max_inliers > 0
%         T_best = estimate_translation_DLT(xs(inliers, :), Xs(inliers, :), R);
%     else
%         warning('No inliers found. Returning initial translation.');
%         T_best = T_init; % Fallback translation
%     end
% end
% 
% 
% % Function: Translation Estimation using DLT
% function T = estimate_translation_DLT(xs, Xs, R)
%     xs = xs';
%     Xs = Xs';
%     N = size(xs, 1); % Number of points
%     M = zeros(3*N,3+N); % Initialize the matrix 
%     b = zeros(3*N, 1);
% 
%     for i = 1:N
%         % First row in the block
%         M(3*i-2, 1) = Xs(4,i)';
%         M(3*i-2, 3+i) = -xs(1,i);
% 
%         % Second row in the block
%         M(3*i-1, 2) = Xs(4,i)';
%         M(3*i-1, 3+i) = -xs(2,i);
% 
%         % Third row in the block
%         M(3*i, 3) = Xs(4,i)';
%         M(3*i, 3+i) = -xs(3,i);
% 
%         % the b part in the block 
% 
%         b(3*(i-1)+1:3*i, 1) = -R *  Xs(1:3,i);
%     end
% 
% 
%     % Solve Using Pseudo-Inverse (if b ≠ 0)
%     x_pinv = pinv(M) * (-b); % Least-squares solution
% 
%     T = x_pinv(1:3);
% 
% end
% 
% % Function: Compute Reprojection Errors
% function errors = compute_reprojection_errors(xs, Xs, R, T)
%     projected = (R * Xs(:,1:3)' + T)'; % R * Xs' = 3xN, T is 3x1
%     projected = projected(:, 1:2) ./ projected(:, 3); % Normalize to 2D
%     errors = sqrt(sum((xs(1:2) - projected).^2, 2)); % Compute Euclidean distance
% end


%
