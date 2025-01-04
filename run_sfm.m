
function [] = run_sfm(dataset_num, varargin)
% RUN_SFM - Structure from Motion Pipeline
% Inputs:
%   dataset_num - (Required) Dataset number to process.
%   varargin    - (Optional) Key-value pairs:
%       'plot_initial_rel_3D' (default: false) - Flag to plot initial relative 3D points.
%       'plot_initial_rel_3D_pair' (default: false) - Flag to plot initial
%       relative 3D points of the initial pair.
%
% Example Usage:
%   run_sfm(2) % Default, no plot
%   run_sfm(2, 'plot_initial_rel_3D', true) % Enable plotting for all
%   relative pairs
%   run_sfm(2, 'plot_initial_rel_3D_pair', true) % Enable plotting for the
%   relative pair

% --- Parse Optional Arguments ---
p = inputParser; % Create parser object
addRequired(p, 'dataset_num'); % Required input
addParameter(p, 'plot_initial_rel_3D', false, @(x) islogical(x)); % Optional parameter
addParameter(p, 'plot_initial_rel_3D_pair', false, @(x) islogical(x)); % Optional parameter

parse(p, dataset_num, varargin{:}); % Parse inputs

% Extract parameters
dataset_num = p.Results.dataset_num;
plot_initial_rel_3D = p.Results.plot_initial_rel_3D; % Optional flag
plot_initial_rel_3D_pair = p.Results.plot_initial_rel_3D_pair; % Optional flag

% reset and load data
% clc; close all;
% Add the vl_toolbox to be accessible in the current directory
run('C:\Users\pain\Desktop\CV_project\vlfeat-0.9.21-bin\vlfeat-0.9.21\toolbox\vl_setup')
% Get dataset info by providing dataset input

% Control random number generator - to get consistant resutls
rng(42)

%% Load the dataset data information
disp("Step 0: Load the dataset data information")

% Extract dataset information
[K, img_names, init_pair, pixel_threshold] = get_dataset_info(dataset_num);

%% Step 1: Calculate the relative Orientation
% This step calculates the relative rotations (R_rel)
% between consecutive image pairs in the dataset.
% It also performs feature extraction, matches features, estimates the essential
% matrix using RANSAC, and triangulates 3D points.

disp("Step 1: Calculate the relative Rotation between image pairs")

% Number of images in the dataset
N = length(img_names);

% Initialize containers(cells) for storing relative rotations
R_rel = cell(1, N-1); % Relative rotations for each image pair

% Containers for descriptors and features
desc_im = cell(1,N);
feat_im = cell(1,N);

% Container(cell) for inlier indices after RANSAC
inliers_indices = cell(1,N-1);

% Loop through image pairs to compute relative rotations
for n=1:N-1
    %%% Feature Extraction using SIFT %%%

    % Load images for the current and next view
    im1 = imread(img_names{n});
    im2 = imread(img_names{n+1});

    % Extract features and descriptors using SIFT for both images
    [x1, x2, ~, fA, dA, fB, dB] = feature_extraction(im1, im2);

    % Assign features and descriptors only once for each image
    if n == 1
        % Store the first image's features and descriptors (processed once)
        feat_im{n} = fA;
        desc_im{n} = dA;
    end

    % Store features and descriptors for the second image in the pair
    feat_im{n+1} = fB;
    desc_im{n+1} = dB;
    

    %%% Estimate the Essential Matrix using RANSAC %%%
  
    % Define the pixel threhold
    eps = pixel_threshold * 2 / (K(1,1) + K(2,2)); % Normalize threshold

    % e.g: 0.000420430520853354

    % Estimate Essential Matrix robustly using RANSAC
    [E_ransac, indices] = estimate_E_robust(x1, x2, eps, K);

    % Store inlier indices for later refinement
    inliers_indices{n} = indices;

    % Triangulate the points 3D points using the relative pose and inliers
    % and Check for Cheirality Condition to get P2

    [X,P] = Cheirality_triangulate(x1,x2, indices,K,E_ransac);

    % Save relative rotation

    R_rel{n} = K\P{2}(1:3, 1:3);


    if plot_initial_rel_3D
        % Optional: Step 5: Visualization
        figure;

        plot3(X(1, :), X(2, :), X(3, :), '.', 'MarkerSize', 10); % Plot 3D points

        hold on;
        plotcams(P); % Plot cameras
        xlabel('x'); ylabel('y'); zlabel('z');
        title('Relative plot of 3D points');
        axis equal;  % Equal scaling for axes
        grid on;
        hold off;
    end
end


%% Step 2: Upgrade to Absolute Rotations
disp("Step 2: Upgrade to Absolute Rotations")

% Get the total number of relative rotations
num_rel_rotations = length(R_rel); % Total relative rotations (N-1)

% Initialize cell array for storing absolute rotations
R_abs_i = cell(1, N);            % Absolute rotations for N cameras
R_abs_i{1} = eye(3);             % First camera has identity rotation (reference)

% Compute absolute rotations from relative rotations
for k = 1:num_rel_rotations
    % Compute absolute rotation by chaining relative rotations
    R_abs_i{k+1} = R_rel{k} * R_abs_i{k}; % Accumulate rotations

    % Display progress for debugging
    fprintf('Computed Absolute Rotation (Camera %d):\n', k+1);
    disp(R_abs_i{k+1});
end

% Final output
disp("Absolute Rotations Computed for all Cameras!");

% Verify orthogonality (R * R' = I)
fprintf('\nVerifying Orthogonality of Rotations:\n');
for k = 1:N 
    error_orthogonality = norm(R_abs_i{k} * R_abs_i{k}' - eye(3)); % Should be close to 0 using Euclidean norm
    fprintf('Camera %d Orthogonality Error: %.6f\n', k, error_orthogonality);
end

%% Step 3: Reconstruct initial 3D points from initial image pair
disp("Step 3: Reconstruct initial 3D points from initial image pair")

% get the 3D points from the suggested pair
% Load the images
im1 = imread(img_names{init_pair(1)});
im2 = imread(img_names{init_pair(2)});

% Pipeline for extracting the featuers,

[x1, x2, desc_X, ~,~,~,~] = feature_extraction(im1,im2);

% Step: Estimate E robust using RANSAC

% Define the pixel threhold
eps = pixel_threshold * 2 / (K(1,1) + K(2,2)); % Normalize threshold

% Get the estimate_E_robust
[E_ransac, indices] = estimate_E_robust(x1, x2, eps, K);

% Save descriptors for future use "Only inliers"
desc_X = desc_X(:,indices);
% save('desc_X.mat', 'desc_X'); % Save to a mat file

% Triangulate points using relative pose

[X_0,P] = Cheirality_triangulate(x1,x2, indices,K,E_ransac);


% Rotate the initial 3D points to the world coordinate using the first
% absolute roation matrix in the initial pair
X_wc = R_abs_i{init_pair(1)}' * X_0(1:3,:);

% Step 1: Compute Center of Gravity
X_mean = mean(X_wc, 2); % 3x1 mean position

% Step 2: Compute Distances
distances = vecnorm(X_wc - X_mean, 2, 1); % Euclidean distances

% Step 3: Compute Threshold
q90 = quantile(distances, 0.9); % 90th percentile
threshold = 5 * q90;           % 5 times the 90% quantile

% Step 4: Filter Outliers
inliers = distances <= threshold;           % Logical indices for inliers
X_wc_filtered = X_wc(:, inliers);           % Filtered 3D points

if plot_initial_rel_3D_pair
    % Step 5: Visualization for initial pair
    figure;
    plot3(X_wc_filtered(1, :), X_wc_filtered(2, :), X_wc_filtered(3, :), '.', 'MarkerSize', 10);
    hold on;
    plotcams(P); % Plot cameras
    xlabel('x'); ylabel('y'); zlabel('z');
    title('Filtered 3D Points and Cameras');
    axis equal;
    grid on;
    hold off;
end

%%
% Step 4:
disp("Step 4: T robust estimation")

% Establish correspondences between i and 3D points (using desc X),
% Number of images
% Initialize translation estimates
T_best = cell(1, N);
T_best{1} = zeros(3, 1);

for n = 1:N-1
    % Step 1: Establish correspondences using descriptor matching
    [matches_2d_3d, ~] = vl_ubcmatch(desc_im{n}, desc_X);

    unique_indices = unique(matches_2d_3d(2, :));

    disp(['Unique 3D Points: ', num2str(length(unique_indices))]);
    disp(['Total Matches: ', num2str(size(matches_2d_3d, 2))]);
    disp(['Percentage%:  ', 100 * num2str(length(unique_indices)) / num2str(size(matches_2d_3d, 2))] )

    [~, ia] = unique(matches_2d_3d(2, :), 'stable');
    filtered_matches = matches_2d_3d(:, ia);

    % Step 2: Extract 2D points and convert to homogeneous coordinates
    xA = feat_im{n}(1:2, filtered_matches(1, :));
    x1 = [xA; ones(1, length(xA))];
    xsn1 = K \ x1; % Normalize 2D points

    % Step 3: Extract 3D points and convert to homogeneous
    Xs = X_wc(:, filtered_matches(2, :)); % Use pre-computed 3D points
    Xsn = [Xs; ones(1, length(Xs))];

    % Step 4: Robust estimation of translation
    translation_threshold = 3 * pixel_threshold / K(1, 1);

    disp("Determinant of R: ")
    disp(det(R_abs_i{n}))

    % Estimate T using Robust RANSAC + DLT
    %TODO: should we prodive the initial as zero ? zeros(3, 1)
    T_best{n+1} = estimate_T_robust(xsn1, Xsn, R_abs_i{n}, translation_threshold, zeros(3, 1));
end

%%
disp("Step 5: Plot the cameras")

P_all = cell(1,N);

for i=1:N
    P_all{i} = K * [R_abs_i{i} T_best{i}];
end

figure();
plotcams(P_all)
title('Camera Poses Visualization');

disp("Cameras plotted successfully!");
%%
% Step 6: Triangulate points for all pairs (i, i + 1) and visualize 3D points + cameras

disp("% Step 6: Triangulate points for all pairs (i, i + 1) and visualize 3D points + cameras")
X_all = [];  % To store all 3D points

for i=1:N-1

    % Compute the matches
    matches = vl_ubcmatch(desc_im{n},desc_im{n+1});

    % Get the 2D points correspondences on each images' pair
    xA = feat_im{n}(1:2, matches(1 ,:));
    xB = feat_im{n+1}(1:2 , matches(2 ,:));

    % Convert the vectors to be homogenous for this pair of images
    % Add a homogenous dimension to the points
    x1 = [xA; ones(1,length(xA))];
    x2 = [xB; ones(1,length(xB))];

    % % Normalize the points
    x1_n = K\x1;
    x2_n = K\x2;

    % Keep only the inliers
    x1_n_in = x1_n(:, inliers_indices{n});
    x2_n_in = x2_n(:, inliers_indices{n});

    % Triangulate the 3D points using the cameras and the 2D points

    X = triangulate_3D_point_DLT(x1_n_in, x2_n_in, P_all{i},  P_all{i+1});
    
    figure();
    % Plot 3D points
    plot3(X(1, :), X(2, :), X(3, :), '.', 'MarkerSize', 10);
    hold on;

    % % Plot cameras
    plotcams(P_all);

    % Labels and axes
    xlabel('X'); ylabel('Y'); zlabel('Z');
    title('3D Points and Cameras');
    axis equal;
    grid on;
    hold off;

    % % Levenberg Marquardt Method
    % 
    % % Compute initial reprojection error
    % total_error_before = ComputeTotalError(P_all{i},  P_all{i+1}, X, x1_n_in, x2_n_in);
    % 
    % disp(['Initial Total sum Error: ', num2str(sum(total_error_before))]);
    % disp(['Initial Total median Error: ', num2str(median(total_error_before / size(X, 2)))]);
    % 
    % % Parameters for LM optimization
    % mu = 1e-1; % Damping factor
    % max_iterations = 10; % Maximum LM iterations
    % tolerance = 1e-6; % Convergence tolerance
    % 
    % % Run Levenberg-Marquardt optimization
    % [X_refined, accumulated_error] = LevenbergMarquardt(P_all{i},  P_all{i+1}, X, x1_n_in, x2_n_in, mu, max_iterations, tolerance);
    % 
    % % Compute the final reprojection error
    % total_error_after = ComputeTotalError(P_all{i},  P_all{i+1}, X_refined, x1_n_in, x2_n_in);
    % 
    % disp(['After Total sum Error: ', num2str(sum(total_error_after))]);
    % disp(['After Total median Error: ', num2str(median(total_error_after / size(X_refined, 2)))]);

end


end


%% Feature Extraction Function %%
function [x1,x2, desc_X, fA, dA, fB, dB] = feature_extraction(im1,im2)
        
        % Return a homogenous projections

        % Extract the features using SIFT
        % 1 specifies a threshold for rejecting low-contrast keypoints detected by SIFT.

        [fA dA] = vl_sift(single(rgb2gray(im1)), 'PeakThresh', 1);
        [fB dB] = vl_sift(single(rgb2gray(im2)), 'PeakThresh', 1);
        
        % Compute the matches
        matches = vl_ubcmatch(dA, dB);
        
        % Get the 2D points correspondences on each images' pair
        xA = fA(1:2 , matches(1 ,:));
        xB = fB(1:2 , matches(2 ,:));
        
        % Convert the vectors to be homogenous for this pair of images 
        % Add a homogenous dimension to the points
        x1 = [xA; ones(1,length(xA))];
        x2 = [xB; ones(1,length(xB))];
    
        % Save descriptors for matched points (3D points will correspond to these)
        % Extract descriptors only for matched points in image 1
        desc_X = dA(:, matches(1, :)); 
       
end

%% estimate_E_robust

function [E_best, indices] = estimate_E_robust(x1, x2, eps, K)
    % Robustly estimates the essential matrix E using RANSAC.

    % x1: 3xN , N: Number of points
    % x2: 3XN , N: Number of points
    
    % Normalize points
    x1_normalized = K\x1;
    x2_normalized = K\x2;

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
    fprintf('Percentage of inliers: %.2f%%\n', 100 * (best_inlier_count) / num_points);
    indices = find(inliers);
end
%%

function F_norm_esit = estimate_F_DLT(x1s, x2s)
% Function to compute the normalized fundamental matrix using the 8-point algorithm
% The output F will be normalized by frobenius norm

% Get the dimenion of the M matrix
n = length(x1s);

% Initialize M
M = zeros(n,9);

% Loop over the points 
for i=1:n
    mul = x2s(:,i) * x1s(:,i)';
    mul = reshape(mul.', [], 1);
    mul = mul';
    M(i,:) = mul;
end

% Calculate the S, V, U of Singular value decomposition
[~,~,V] = svd(M);

% v is the last column of V
v = V(:, end); % 9 elements

% disp('||Mv||:'); disp(Mv_abs); 

F_norm_esit = reshape(v,[3, 3])';

% F_norm_esit = F_norm_esit / norm(F_norm_esit, 'fro'); % Ask about this later?

end
%%

function E = enforce_essential(E_approx)
% ENFORCE_ESSENTIAL Apply the enforce on E

% Apply the SVD on E approximate 
[U,~,V] = svd( E_approx );

if det(U*V') < 0 
    V = -V; 
end

% Creates a valid essential matrix from an approximate solution 

E = U * diag ([1 1 0])* V';

end

%%
function distances = compute_epipolar_errors(F, x1s, x2s)
    % Compute epipolar errors for a given fundamental matrix F and point correspondences
    %
    % Inputs:
    %   F - Fundamental matrix (3x3)
    %   x1s - Points in the first image (3xN, homogeneous coordinates)
    %   x2s - Points in the second image (3xN, homogeneous coordinates)
    %
    % Output:
    %   distances - Vector of epipolar distances for each point correspondence
    % Compute epipolar lines for points in the second image
    l = F * x1s; % Epipolar lines in the second image

    % Normalize the epipolar lines to have unit norm
    l = l ./ sqrt(l(1, :).^2 + l(2, :).^2);

    % Compute distances from points in the second image to their epipolar lines
    distances = abs(sum(l .* x2s, 1));
   
end
%%

function [X_best,P] = Cheirality_triangulate(x1,x2, inliers,K,E_ransac)

  % Include only the inliers on the triangulation
        x1_inliers = x1(:,inliers);
        x2_inliers = x2(:,inliers);
        
        % Correct camera P2
        P2_correct_ind = 0 ;
        max_points_1 = 0;
        
        % X best 
        X_best = [];
        
        % Define the cameras & assume first camera is at the origin
        P1 = K* [eye(3), zeros(3, 1)];
        
        % Get the U and V canditates
        [U, ~, V] = svd(E_ransac);
        W = [0 -1 0; 1 0 0; 0 0 1]; % Rotation matrix for decomposition
        
        % Four possible solutions for P2
        R1 = U * W * V';
        R2 = U * W' * V';
        t = U(:, 3);
        
        % Ensure proper rotations (determinant = +1)
        if det(R1) < 0, R1 = -R1; end
        if det(R2) < 0, R2 = -R2; end
        
        % Four possible camera matrices
        P2_candidates = {
            K * [R1, t],
            K * [R1, -t],
            K * [R2, t],
            K * [R2, -t]
        };
                
        % Choose the correct P2 by testing chirality (points must be in front of both cameras)
        for i = 1:4
            P2 = P2_candidates{i};
        
            X = triangulate_3D_point_DLT(x1_inliers, x2_inliers,P1, P2);
            x1_proj = P1 * X; 
            x2_proj = P2 * X; 
        
            % This is a general condition for the position of the points 
            cond = sum( (x1_proj(3,:) > 0) & ((x2_proj(3,:) > 0)) );
            
            if cond > max_points_1
                max_points_1 = cond;
                P2_correct_ind = i;
                X_best = X;
            end
        end
        
        P{1} = P1; 
        P{2} = P2_candidates{P2_correct_ind};
end

%%
function X = triangulate_3D_point_DLT(x1, x2, P1, P2)
    % INPUT:
    % x1 - 2D points in the first image (3xN homogeneous coordinates)
    % x2 - 2D points in the second image (3xN homogeneous coordinates)
    % P1 - Camera projection matrix for the first image (3x4)
    % P2 - Camera projection matrix for the second image (3x4)
    % OUTPUT:
    % X - Triangulated 3D points (4xN homogeneous coordinates)

    % Number of points given in x1 and x2 
    N = size(x1, 2);

    X = zeros(4, N); 
    
    % Loop through all points to get all 3D points
    for i = 1:N
        % Construct the A matrix for triangulation
        A = [
            x1(1,i) * P1(3,:) - P1(1,:);
            x1(2,i) * P1(3,:) - P1(2,:);

            x2(1,i) * P2(3,:) - P2(1,:);
            x2(2,i) * P2(3,:) - P2(2,:);
        ];
        
        % Solve the homogeneous system using SVD as before and we need only
        % V to extract v 

        [~, ~, V] = svd(A);
        X(:,i) = V(:,end); % Last column of V gives the solution ( last eigenvector)
    end
    
    % Normalize the points to make them homogeneous (4th coordinate = 1)
    X = X(1:4, :) ./ X(4, :);
end

%%

% Function: Robust Translation Estimation
function T_best = estimate_T_robust(xs, Xs, R, inlier_threshold, T_init)

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
        inliers_candidate = errors.^2 < inlier_threshold.^2;
        num_inliers = sum(inliers_candidate);

        % Update best model
        if num_inliers > max_inliers
            max_inliers = num_inliers;
            T_best = T_candidate;
            inliers = inliers_candidate;
        end
    end

    % Refine with inliers
    if max_inliers > 0
        T_best = estimate_translation_DLT(xs(:, inliers), Xs(:, inliers), R);
    else
        warning('No inliers found. Returning initial translation.');
        T_best = T_init;
    end
end


% Function: Translation Estimation using DLT
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
    % residual = norm(A * T_lambda - b); % ||A * T_lambda - b||
    % disp('Residual Error:');
    % disp(residual);

    % Flip sign if necessary based on direction consistency
    if dot(T, mean(X(1:3, :), 2)) < 0
        T = -T;  % Flip the translation vector
    end
end


% Function: Compute Reprojection Errors
function errors = compute_reprojection_errors(xs, Xs, R, T)
    % Project 3D points
    projected = (R * Xs(1:3, :) + T);
    projected = projected(1:2, :) ./ projected(3, :); % Normalize

    % Compute errors
    errors = sqrt(sum((xs(1:2, :) - projected).^2, 1)); % Euclidean distance
end



%% Other helper functions


function [V] = pflat(V)
% Normalize the n points in P^(m-1)

V = V(:,:) ./ V(end, :);

end
%%
function plotcams(P)
    % Function to plot camera positions and viewing directions
    % Inputs:
    %   P - Cell array of camera projection matrices
    
    % Initialize arrays to store camera centers and directions
    c = zeros(4, length(P)); % Camera centers
    v = zeros(3, length(P)); % Viewing directions
    
    for i = 1:length(P)
        c(:,i) = null(P{i});          % Compute camera center (homogeneous coordinates)
        v(:,i) = P{i}(3, 1:3);        % Viewing direction (3rd row of P)
    end
    
    % Normalize homogeneous coordinates
    c = c ./ repmat(c(4,:), [4 1]);   % Convert to 3D by dividing by the last element

    % Plot cameras using quiver3
    quiver3(c(1,:), c(2,:), c(3,:), v(1,:), v(2,:), v(3,:), 'r-'); 
    hold on;
    

    %Label all cameras
    for i = 1:length(P)
        label = sprintf('Camera %d', i);
        text(c(1,i), c(2,i), c(3,i), label, 'FontSize', 1, 'Color', 'b');
    end
    
    % Enhance the plot for clarity
    xlabel('X');
    ylabel('Y');
    zlabel('Z');
    grid on;
    axis equal;
    hold off;
end
%%