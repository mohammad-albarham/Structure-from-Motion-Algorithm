function T = estimate_translation_DLT(xs, Xs, R)
    % Inputs:
    % xs - 2D points in homogeneous coordinates (3xN)
    % Xs - 3D points in homogeneous coordinates (4xN)
    % R  - Rotation matrix (3x3)

    % Extract the first two points
    x1 = xs(:,1);
    x2 = xs(:,2);
    X1 = Xs(:,1);
    X2 = Xs(:,2);

    % Construct matrix directly using skew-symmetric form
    M = [
        skew(x1) * R * X1(1:3), skew(x1);
        skew(x2) * R * X2(1:3), skew(x2)
    ];

    % Solve using SVD
    [~, ~, V] = svd(M);
    T = V(2:end, end) ./ V(1,4); % Extract translation vector (last column, rows 4-6)
end

function S = skew(x)
    % Constructs the skew-symmetric matrix for cross product
    S = [  0   -x(3)  x(2);
          x(3)   0   -x(1);
         -x(2)  x(1)   0 ];
end