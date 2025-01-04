% Generate synthetic data
N = 100; % Number of 3D points
X_gt = [rand(3, N) * 10; ones(1, N)]; % 3D points in homogeneous coordinates

% Ground Truth Translation (T_gt) and Rotation (R_gt)
T_gt = [1; 2; 3];         % Ground truth translation
R_gt = eye(3);            % Identity rotation (no rotation)

% Intrinsics
K = eye(3);               % Assume identity intrinsic matrix for simplicity

% Project 3D points into 2D
x = K * [R_gt, T_gt] * X_gt;
x = x ./ x(3, :);         % Normalize

% Normalize points
x_norm = K \ x;

sample_indices = randperm(N, 2);
xs_sample = x(:, sample_indices);
Xs_sample = X_gt(:, sample_indices);

% Estimate Translation using DLT
T_est = estimate_translation_DLT(xs_sample, Xs_sample, R_gt);

% Display results
disp('Estimated T:');
disp(T_est);
disp('Ground Truth T:');
disp(T_gt);

% Compute translation error
error = norm(T_gt - T_est);
disp('Translation Error:');
disp(error);



function T = estimate_translation_DLT(x, X, R)
    % Inputs:
    % x - 2D normalized points in homogeneous coordinates (3xN)
    % X - 3D points in homogeneous coordinates (4xN)
    % R - Rotation matrix (3x3)

    % Outputs:
    % T - Estimated translation vector (3x1)

    % Number of points
    N = size(x, 2);

    % Initialize matrices A and b
    A = [];
    b = [];

    % Loop through each point
    for i = 1:N
        % Extract the 2D point (x) and 3D point (X)
        xi = x(:, i);       % 2D point in homogeneous coordinates
        Xi = X(:, i);       % 3D point in homogeneous coordinates

        % Create the linear system based on the equations in the image
        A_i = [
            Xi(4), 0, 0, -xi(1);   % First row
            0, Xi(4), 0, -xi(2);   % Second row
            0, 0, Xi(4), -xi(3)    % Third row
        ];

        % Append to the full matrix A
        A = [A; A_i];

        % Compute b using the rotation matrix R
        b_i = -R * Xi(1:3);  % Projected 3D point
        b = [b; b_i];
    end

    % Solve the linear system A * T = b using least squares
    T_lambda = A \ b;

    % Extract translation vector T (first 3 elements)
    T = T_lambda(1:3);

    % Compute residual error
    residual = norm(A * T_lambda - b); % ||A * T_lambda - b||
    disp('Residual Error:');
    disp(residual);

    % Flip sign if necessary based on direction consistency
    if dot(T, mean(X(1:3, :), 2)) < 0
        T = -T;  % Flip the translation vector
    end
end

%%
% Load example images
img1 = imread("data/2/DSC_0025.JPG");  % Replace with your image
img2 = imread("data/2/DSC_0033.JPG");  % Replace with your image

% Detect features in both images
gray1 = rgb2gray(img1);
gray2 = rgb2gray(img2);
points1 = detectSURFFeatures(gray1);
points2 = detectSURFFeatures(gray2);

% Extract features
[features1, validPoints1] = extractFeatures(gray1, points1);
[features2, validPoints2] = extractFeatures(gray2, points2);

% Match features
indexPairs = matchFeatures(features1, features2);
matchedPoints1 = validPoints1(indexPairs(:, 1));
matchedPoints2 = validPoints2(indexPairs(:, 2));

% Define camera intrinsics
focalLength = [2.3785e+03, 2.3785e+03];    % Replace with actual values
principalPoint = [968, 648];   % Replace with actual values
imageSize = size(gray1);
intrinsics = cameraIntrinsics(focalLength, principalPoint, imageSize(1:2));

% Estimate the Essential matrix
[E, inliers] = estimateEssentialMatrix(...
    matchedPoints1, matchedPoints2, intrinsics);

% Select inlier matches
inlierPoints1 = matchedPoints1(inliers);
inlierPoints2 = matchedPoints2(inliers);

% Recover relative camera pose
[R, T] = relativeCameraPose(E, intrinsics, inlierPoints1, inlierPoints2);

% Display results
disp('Estimated Translation (T):');
disp(T);

% Display results
disp('Estimated Rotation (R):');
disp(R);