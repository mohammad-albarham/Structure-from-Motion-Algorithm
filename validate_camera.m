function validate_camera(P)
% VALIDATE_CAMERA - Validates the projection matrix of a camera
%
% Inputs:
%   P - 3x4 Camera projection matrix
%   K - 3x3 Intrinsic parameter camera matrix
% Outputs:
%   Prints validation results for R and T:
%   - Orthogonality error of R
%   - Determinant of R
%   - Translation vector (T) and its norm

    info("\n--- validate_camera Validation Start ---\n",2);

    
    % --- Extract Rotation (R) and Translation (T) ---
    R = P(1:3, 1:3); % Rotation matrix
    T = P(1:3, 4);   % Translation vector

    % --- Validate Orthogonality of R ---
    orthogonality_error = norm(R * R' - eye(3)); % Should be close to 0

    info("R Orthogonality Error: %.6e\n",2,orthogonality_error);


    % --- Validate Determinant of R ---
    det_R = det(R); % Should be close to +1
    info("Determinant of R: %.6f\n",2, det_R);

    % --- Display Translation Information ---
    info("Translation Vector T: [%.4f, %.4f, %.4f]\n",2, T(1), T(2), T(3));
    info("Norm of T: %.4f\n",2, norm(T));

    % --- Additional Checks ---
    % Warn if errors exceed acceptable thresholds
    if orthogonality_error > 1e-6
        warning('R is not orthogonal. Check rotation matrix!');
    end
    if abs(det_R - 1) > 1e-3
        warning('Determinant of R is not 1. Check rotation matrix!');
    end
    
    info("\n--- validate_camera Validation End ---\n",2, norm(T));

end
