function validate_E(E)
% VALIDATE_E - Validates the Essential matrix E based on rank, singular values, 
% epipolar constraint, and determinant.
%
% Inputs:
%   E  - Essential matrix (3x3)
%   x1 - Normalized homogeneous coordinates of points in image 1 (3xN)
%   x2 - Normalized homogeneous coordinates of points in image 2 (3xN)
%
% Example:
% validate_E(E_candidate, x1_h_n_sample, x2_h_n_sample);

    % fprintf('\n--- Essential Matrix Validation Start ---\n');
    info("\n--- Essential Matrix Validation Start ---\n", 2);

    % 1. Rank Constraint
    rank_E = rank(E); % Compute rank
    % fprintf('Rank of E: %d (Expected: 2)\n', rank_E);

    info("Rank of E: %d (Expected: 2)\n:  ", 2, rank_E);



    % 2. Singular Value Constraint
    [U, S, V] = svd(E); % SVD
    singular_values = diag(S);
    % fprintf('Singular values of E: [%.6f, %.6f, %.6f]\n', singular_values);

    info("Singular values of E: [%.6f, %.6f, %.6f]\n", 2, singular_values);


    if abs(S(1,1) - S(2,2)) < 1e-6 && S(3,3) < 1e-6
        % fprintf('Singular value constraint: Satisfied\n');
        info("Singular value constraint: Satisfied\n", 2);

    else
        % fprintf('Singular value constraint: Not satisfied\n');
        info("Singular value constraint: Not satisfied\n", 2);

    end

    % 4. Determinant Check
    det_E = det(E); % Compute determinant
    % fprintf('Determinant of E: %.6f (Expected: Close to 0)\n', det_E);
    info("Determinant of E: %.6f (Expected: Close to 0)\n", 2, det_E);


    % Final Assessment
    if rank_E == 2 && abs(S(1,1) - S(2,2)) < 1e-6 && S(3,3) < 1e-6 && abs(det_E) < 1e-6
        % fprintf('Validation Status: PASS - Essential matrix is valid.\n');
        info("Validation Status: PASS - Essential matrix is valid.\n", 2);

    else
        % fprintf('Validation Status: FAIL - Essential matrix is invalid.\n');
        info("Validation Status: FAIL - Essential matrix is invalid.\n", 2);
    end

    % fprintf('\n--- Essential Matrix Validation End ---\n');
    info("\n--- Essential Matrix Validation End ---\n", 2);

end
