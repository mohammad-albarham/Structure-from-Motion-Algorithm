function [X_best, P] = Cheirality_triangulate(x1_h_n_in, x2_h_n_in, E_ransac)
% CHEIRALITY_TRIANGULATE - Triangulate 3D points and resolve E ambiguity
% Inputs:
%   x1, x2      - 2D unnormalized homogeneous points in images 1 and 2 (3xN)
%   inliers     - Indices of inlier matches after RANSAC
%   K           - Intrinsic camera matrix (3x3)
%   E_ransac    - Estimated essential matrix (3x3)
% Outputs:
%   X_best      - Triangulated 3D points (4xN homogeneous coordinates)
%   P           - Correct camera projection matrices {P1, P2}

    % Use the selected inlier points for triangulation
    
    % Initialize best parameters
    P2_correct_ind = 0; % Index for correct P2
    max_points_1 = 0;   % Maximum valid points satisfying chirality
    X_best = [];        % Store the best 3D points

    % Step 2: Define camera matrices
    % First camera at the origin (canonical form)
    P1 = [eye(3), zeros(3, 1)];

    % Step 3: Decompose Essential Matrix (E) to get R, t candidates
    [U, ~, V] = svd(E_ransac);
    W = [0 -1 0; 1 0 0; 0 0 1]; % Pre-defined rotation matrix

    % Four possible solutions for [R, t]
    R1 = U * W * V';
    R2 = U * W.'*V';
    
    t = U(:, 3);

    % Ensure proper rotations (determinant = +1)
    if det(R1) < 0, R1 = -R1; end
    if det(R2) < 0, R2 = -R2; end

    % Step 4: Generate P2 candidates (Four cases)
    P2_candidates = {
        [R1, t],
        [R1, -t],
        [R2, t],
        [R2, -t]
    };

    % Step 5: Evaluate each P2 candidate based on chirality condition
    for i = 1:4
        P2 = P2_candidates{i};

        % Triangulate 3D points using current P2
        X = triangulate_3D_point_DLT(x1_h_n_in, x2_h_n_in, P1, P2);

        % Project back to both views
        x1_proj = P1 * X;
        x2_proj = P2 * X;

        % Count points in front of both cameras (chirality check)
        valid_points = (x1_proj(3,:) > 0) & (x2_proj(3,:) > 0);
        num_valid = sum(valid_points);


        % Keep the candidate with the most valid points
        if num_valid > max_points_1
            max_points_1 = num_valid;
            P2_correct_ind = i;
            X_best = X; % (:, valid_points); % Store the best 3D points % Filter valid 3D points
        end
    end
    
 
    % Step 6: Final output - Correct projection matrices
    P{1} = P1;
    P{2} = P2_candidates{P2_correct_ind};

    info('\n--- Cheirality_triangulate Validation Start ---\n',2);

    % Print the selected solution
    info('Selected P2 Index: %d\n',2,P2_correct_ind);
    info('Number of Valid Points: %d out of %d\n',2, max_points_1, length(x1_h_n_in));

    % Optional: Compute reprojection error
    err1 = compute_reprojection_error_P(P{1}, X_best, x1_h_n_in); % (:, valid_points));
    err2 = compute_reprojection_error_P(P{2}, X_best, x2_h_n_in); % (:, valid_points));
    info('Mean Reprojection Error: %.4f (Image 1), %.4f (Image 2)\n',2, mean(err1), mean(err2));
    info('\n--- Cheirality_triangulate Validation End ---\n',2);
end