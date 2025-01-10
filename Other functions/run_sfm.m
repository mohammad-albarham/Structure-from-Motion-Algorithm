
% function [] = run_sfm(dataset_num, varargin)
% RUN_SFM - Structure from Motion Pipeline
%% 
% Inputs:
%   dataset_num - (Required) Dataset number to process.
%   varargin    - (Optional) Key-value pairs:
%       'plot_initial_rel_3D' (default: false) - Flag to plot initial relative 3D points.
%       'plot_initial_rel_3D_pair' (default: false) - Flag to plot initial
%       relative 3D points of the initial pair.
%
% Example Usage:
%   run_sfm(2) % Default, no plot++
%   run_sfm(2, 'plot_initial_rel_3D', true) % Enable plotting for all
%   relative pairs
%   run_sfm(2, 'plot_initial_rel_3D_pair', true) % Enable plotting for the
%   relative pair

% --- Parse Optional Arguments ---
% p = inputParser; % Create parser object
% addRequired(p, 'dataset_num'); % Required input
% addParameter(p, 'plot_initial_rel_3D', false, @(x) islogical(x)); % Optional parameter
% addParameter(p, 'plot_initial_rel_3D_pair', false, @(x) islogical(x)); % Optional parameter
% 
% parse(p, dataset_num, varargin{:}); % Parse inputs

% Extract parameters
% dataset_num = p.Results.dataset_num;
% plot_initial_rel_3D = p.Results.plot_initial_rel_3D; % Optional flag
% plot_initial_rel_3D_pair = p.Results.plot_initial_rel_3D_pair; % Optional flag
close all; clear; clc
dataset_num = 15;
plot_initial_rel_3D = false; 
plot_initial_rel_3D_pair = false;

global enableInfo;
global info_level;

enableInfo = true; % Enable or disable info globally
info_level = 1;
%% 

% reset and load data
% clc; close all;
% Add the vl_toolbox to be accessible in the current directory
run('.\vlfeat-0.9.21-bin\vlfeat-0.9.21\toolbox\vl_setup')
% Get dataset info by providing dataset input

% Control random number generator - to get consistant resutls
rng(42)

%% Load the dataset data information
% disp("Step 0: Load the dataset data information")
info("Step 0: Load the dataset data information")


% Extract dataset information
[K, img_names, init_pair, pixel_threshold] = get_dataset_info(dataset_num);

%% Step 1: Calculate the relative Orientation
% This step calculates the relative rotations (R_rel)
% between consecutive image pairs in the dataset.
% It also performs feature extraction, matches features, estimates the essential
% matrix using RANSAC, and triangulates 3D points.

% disp("Step 1: Calculate the relative Rotation between image pairs")


info("Step 1: Calculate the relative Rotation between image pairs")

% Number of images in the dataset
N = length(img_names);

% Initialize containers(cells) for storing relative rotations
R_rel = cell(1, N-1); % Relative rotations for each image pair

% Containers for descriptors and features
desc_im = cell(1,N);
feat_im = cell(1,N);
matches_im = cell(1,N-1);

% Container(cell) for inlier indices after RANSAC
inliers_indices = cell(1,N-1);

% Loop through image pairs to compute relative rotations
for n=1:N-1
    %%% Feature Extraction using SIFT %%%

    % Load images for the current and next view
    im1 = imread(img_names{n});
    im2 = imread(img_names{n+1});

    % Extract features and descriptors using SIFT for both images
    [x1, x2, ~, fA, dA, fB, dB, matches_im{n}] = feature_extraction(im1, im2);

    
    % vectors with dim as
    % 3 x N; 
    % N-> Number of the points
    x1_h = toHomogeneous(x1); 
    x2_h = toHomogeneous(x2);

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
    epipolar_threshold = pixel_threshold / K(1,1); % Normalize threshold

    % e.g: 0.000420430520853354

    x1_h_n = K\x1_h;
    x2_h_n = K\x2_h;

    % Estimate Essential Matrix robustly using RANSAC
    [E_ransac, indices] = estimate_E_robust(x1_h_n, x2_h_n, epipolar_threshold);

    % Store inlier indices for later refinement
    inliers_indices{n} = indices;

    % Triangulate the points 3D points using the relative pose and inliers
    % and Check for Cheirality Condition to get P2
    
    x1_h_n_in = x1_h_n(:,inliers_indices{n});
    x2_h_n_in = x2_h_n(:,inliers_indices{n});
    
    [X,P_pair] = Cheirality_triangulate(x1_h_n_in,x2_h_n_in,E_ransac);

    validate_camera(P_pair{1})

    validate_camera(P_pair{2})

    
    % Save relative rotation

    R_rel{n} = P_pair{2}(1:3, 1:3);


    if plot_initial_rel_3D
        % Optional: Step 5: Visualization
        figure;

        plot3(X(1, :), X(2, :), X(3, :), '.', 'MarkerSize', 5); % Plot 3D points

        hold on;

        P_pair_un = CalibratedToUncalibrated(P_pair, K);
        plotcams(P_pair_un); % Plot cameras
        xlabel('x'); ylabel('y'); zlabel('z');
        title('Relative plot of 3D points');
        axis equal;  % Equal scaling for axes
        grid on;
        hold off;
    end
end


%% Step 2: Upgrade to Absolute Rotations
% disp("Step 2: Upgrade to Absolute Rotations")

info("Step 2: Upgrade to Absolute Rotations")


% Initialize cell array for storing absolute rotations
R_abs_i = cell(1, N);            % Absolute rotations for N cameras
R_abs_i{1} = eye(3);             % First camera has identity rotation (reference)

% Compute absolute rotations from relative rotations
for k = 1:N-1
    % Compute absolute rotation by chaining relative rotations
    R_abs_i{k+1} = R_rel{k} * R_abs_i{k}; % Accumulate rotations

    % Display progress for debugging
    % fprintf('Computed Absolute Rotation (Camera %d):\n', k+1);
    info("Computed Absolute Rotation (Camera %d):\n",k+1);
    disp(R_abs_i{k+1});
end

% Final output
% disp("Absolute Rotations Computed for all Cameras!");
info("Absolute Rotations Computed for all Cameras!");


% Verify orthogonality (R * R' = I)
% fprintf('\nVerifying Orthogonality of Rotations:\n');
info("\nVerifying Orthogonality of Rotations:\n");

for k = 1:N 
    error_orthogonality = norm(R_abs_i{k} * R_abs_i{k}' - eye(3)); % Should be close to 0 using Euclidean norm
    info("Camera %d Orthogonality Error: %.6f\n", 1,k, error_orthogonality);
end

%% Step 3: Reconstruct initial 3D points from initial image pair
% disp("Step 3: Reconstruct initial 3D points from initial image pair")

info("Step 3: Reconstruct initial 3D points from initial image pair:\n");


% get the 3D points from the suggested pair
% Load the images
im1 = imread(img_names{init_pair(1)});
im2 = imread(img_names{init_pair(2)});

% Pipeline for extracting the featuers,

[x1, x2, desc_X, ~,~,~,~,~] = feature_extraction(im1, im2);


% vectors with dim as
% 3 x N; 
% N-> Number of the points
x1_h = toHomogeneous(x1); 
x2_h = toHomogeneous(x2);

x1_h_n = K\x1_h;
x2_h_n = K\x2_h;

% Get the estimate_E_robust
[E_ransac, indices] = estimate_E_robust(x1_h_n, x2_h_n, epipolar_threshold);


% save('desc_X.mat', 'desc_X'); % Save to a mat file

% Triangulate points using relative pose

x1_h_n_in = x1_h_n(:,indices);
x2_h_n_in = x2_h_n(:,indices);

[X_init_pair,P_init_pair] = Cheirality_triangulate(x1_h_n_in,x2_h_n_in,E_ransac);

validate_camera(P_init_pair{1})

validate_camera(P_init_pair{2})

% X_init_pair = removePointsExcessivelyFarAway(X_init_pair);
% 
% 
% new_indices = find(indices < length(X_init_pair));
% 
% % Save descriptors for future use "Only inliers"
desc_X_init_pair = desc_X(:,indices);


% Move the 3D pionts to the 3D world coordinate 
X_init_pair_wc = R_abs_i{init_pair(1)}' * X_init_pair(1:3, :);
X_init_pair_wc_h = toHomogeneous(X_init_pair_wc);

%%
if true
    % Step 5: Visualization for initial pair
    figure;
    plot3(X_init_pair_wc_h(1, :), X_init_pair_wc_h(2, :), X_init_pair_wc_h(3, :), ".");
    hold on;
    xlabel('x'); ylabel('y'); zlabel('z');
    title('Filtered 3D Points');
    xlim([-20, 20]);
    ylim([-20, 20]);
    zlim([-20, 20]);
    hold off;
end

%%
% Step 4:
% disp("Step 4: T robust estimation")
info("Step 4: T robust estimation:\n");


% Establish correspondences between i and 3D points (using desc X),
% Number of images
% Initialize translation estimates
T_best = cell(1, N);
% T_best{1} = zeros(3, 1);

for n = 1:N
    % Step 1: Establish correspondences using descriptor matching
    [matches_2d_3d, ~] = vl_ubcmatch(desc_X_init_pair, desc_im{n});

    % unique_indices = unique(matches_2d_3d(2, :));
    % disp(['Unique 3D Points: ', num2str(length(unique_indices))]);
    % disp(['Total Matches: ', num2str(size(matches_2d_3d, 2))]);
    % info("Percentage of unique matches : %.2f%%\n:  ", 100 * length(unique_indices) / size(matches_2d_3d, 2));


    % [~, ia] = unique(matches_2d_3d(2, :), 'stable');
    % filtered_matches = matches_2d_3d(:, ia);

    % Step 2: Extract 2D points and convert to homogeneous coordinates
    xA = feat_im{n}(1:2, matches_2d_3d(2, :));

    x1_h = toHomogeneous(xA);

    x1_h_n = K \ x1_h; % Normalize 2D points

    % Step 3: Extract 3D points and convert to homogeneous
    Xs_h = X_init_pair_wc_h(:, matches_2d_3d(1, :)); % Use pre-computed 3D points
    
    % Step 4: Robust estimation of translation
    translation_threshold = 3 * pixel_threshold / K(1, 1);

    disp("Determinant of R: ")
    disp(det(R_abs_i{n}))

    % Estimate T using Robust RANSAC + DLT
    %TODO: should we prodive the initial as zero ? zeros(3, 1)
    T_best{n} = estimate_T_robust(x1_h_n, Xs_h, R_abs_i{n}, translation_threshold);


    % Normalize it to avoid scale issues
    % T_best{n+1} = T_best{n+1} / norm(T_best{n+1});
end

%%
% disp("Step 5: Plot the cameras")

info("Step 5: Plot the cameras: \n:  ");


P_all = cell(1,N);

for i=1:N
    P_all{i} = [R_abs_i{i} T_best{i}]; 
end

figure();

P_all_un = CalibratedToUncalibrated(P_all, K);
plotcams(P_all_un);

title('Camera Poses Visualization');

% disp("Cameras plotted successfully!");
% info("Cameras plotted successfully!: \n:  ", 100 * length(unique_indices) / size(matches_2d_3d, 2));

%%

% Step 6: Triangulate points for all pairs (i, i + 1) and visualize 3D points + cameras

% disp("% Step 6: Triangulate points for all pairs (i, i + 1) and visualize 3D points + cameras")
info("Step 6: Triangulate points for all pairs (i, i + 1) and visualize 3D points + cameras: \n:  ");

X_all = [];  % To store all 3D points

% Generate a colormap with N-1 distinct colors
colormap = lines(N-1);

% Initialize a cell array to store 3D points and colors
points_with_colors = cell(N-1, 1);

for n = 1:N-1
    % Compute the matches
    % matches = vl_ubcmatch(desc_im{n}, desc_im{n+1});
    
    % Get the 2D points correspondences on each images' pair
    x1 = feat_im{n}(1:2, matches_im{n}(1, :));
    x2 = feat_im{n+1}(1:2, matches_im{n}(2, :));

    % Convert the vectors to homogeneous for this pair of images
    x1_h = toHomogeneous(x1);
    x2_h = toHomogeneous(x2);

    % Normalize the points
    x1_h_n = K \ x1_h;
    x2_h_n = K \ x2_h;

    % Estimate the Essential Matrix and get inliers
    [~, indices] = estimate_E_robust(x1_h_n, x2_h_n, epipolar_threshold);

    % Keep only inliers
    x1_h_n_in = x1_h_n(:, indices);
    x2_h_n_in = x2_h_n(:, indices);

    % Triangulate the 3D points using cameras and the 2D points
    X = triangulate_3D_point_DLT(x1_h_n_in, x2_h_n_in, P_all{n}, P_all{n+1});

    % Remove excessively far away points for clarity
    X = removePointsExcessivelyFarAway(X);

    % Save points and corresponding color in the cell array
    points_with_colors{n} = struct('points', X, 'color', colormap(n, :));
end

%%
% Visualize all points with their colors

% Visualize all points with their colors
figure();
for n = 1:N-1
    % Extract points and color
    X = points_with_colors{n}.points;
    color = points_with_colors{n}.color;

    % Plot the points with the corresponding color
    plot3(X(1, :), X(2, :), X(3, :), ".", 'Color',color);
    hold on;
end


% Plot cameras
plotcams(P_all_un);

% Labels and axes
xlabel('X');
ylabel('Y');
zlabel('Z');
title('3D Points and Cameras - All Cameras');
axis equal;
grid off;
hold off;



% Save the structured data
save('points_with_colors.mat', 'points_with_colors');

% Save the required data
filename = sprintf('sfm_data_%d.mat', dataset_num);
save(filename, 'X', 'P_all', 'R_abs_i', 'T_best', 'inliers_indices');
imagename = sprintf('dataset_image_%d.png', dataset_num);
saveas(gcf, imagename);  % Save as PNG


%% Feature Extraction Function %%
function [x1,x2, desc_X, fA, dA, fB, dB, matches] = feature_extraction(im1,im2)
    % FEATURE_EXTRACTION - Extracts keypoints, descriptors, and matches between two images.
    %
    % Inputs:
    %   im1, im2 - Input images (grayscale or RGB).
    %
    % Outputs:
    %   x1, x2   - 2D coordinates of matched keypoints in both images (2xN).
    %   desc_X   - Descriptors of matched keypoints in the first image (128xN).
    %   fA, dA   - Keypoints and descriptors for the first image.
    %   fB, dB   - Keypoints and descriptors for the second image.

        
        % Return a 2D corresponces and other ....

        % Extract the features using SIFT
        % 1 specifies a threshold for rejecting low-contrast keypoints detected by SIFT.

        [fA, dA] = vl_sift(single(rgb2gray(im1)),'PeakThresh', 2);
        [fB, dB] = vl_sift(single(rgb2gray(im2)),'PeakThresh', 2);
        
        % Compute the matches
        matches = vl_ubcmatch(dA, dB);
        
        % Get the 2D points correspondences on each images' pair
        x1 = fA(1:2 , matches(1 ,:));
        x2 = fB(1:2 , matches(2 ,:));
  
    
        % Save descriptors for matched points (3D points will correspond to these)
        % Extract descriptors only for matched points in image 1
        desc_X = dA(:, matches(1, :)); 
       
end

%% estimate_E_robust

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
%%

function F_norm_esit = estimate_F_DLT(x1_h_n, x2_h_n)
% Function to compute the normalized fundamental matrix using the 8-point algorithm
% The output F will be normalized by frobenius norm

% ESTIMATE_F_DLT - Computes the Fundamental Matrix using the 8-point algorithm.
%
% Inputs:
%   x1s - 2D normalized points in the first image (3xN homogeneous coordinates).
%   x2s - 2D normalized points in the second image (3xN homogeneous coordinates).
%
% Outputs:
%   F_norm_esit - Estimated fundamental matrix (3x3) normalized by Frobenius norm.

% Get the dimenion of the M matrix
n = length(x1_h_n);

% Initialize M
M = zeros(n,9);

% Loop over the points 
for i=1:n
    mul = x2_h_n(:,i) * x1_h_n(:,i)';
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
end
%%

function E = enforce_essential(E_approx)
% ENFORCE_ESSENTIAL Apply the enforce on E

% Apply the SVD on E approximate 
[U,~,V] = svd( E_approx );

% if det(U*V') < 0 
%     V = -V; 
% end

% Check the rank of F
rank_E_approx = rank(E_approx);

% Display the rank
% disp(['Rank of E_approx: ', num2str(rank_E_approx)]);

info("Rank of E_approx: %.2f \n:  ", 2, rank_E_approx);


% Creates a valid essential matrix from an approximate solution 

E = U * diag ([1 1 0])* V';

% Check the rank of F
rank_E = rank(E);

% Display the rank
% disp(['Rank of E: ', num2str(rank_E)]);
info("Rank of E: %.2f \n:  ", 2, rank_E);


end

%%
function distances = compute_epipolar_errors(F, x1_h_n, x2_h_n)
    % Compute epipolar errors for a given fundamental matrix F and point correspondences
    %
    % Inputs:
    %   F - Fundamental matrix (3x3)
    %   x1_h_n - Points in the first image (3xN, homogeneous coordinates) and
    %   normalized
    %   x2_h_n - Points in the second image (3xN, homogeneous coordinates) and
    %   normalized
    %
    % Output:
    %   distances - Vector of epipolar distances for each point correspondence
    % Compute epipolar lines for points in the second image
    l = F * x1_h_n; % Epipolar lines in the second image

    % Normalize the epipolar lines to have unit norm
    l = l ./ sqrt(l(1, :).^2 + l(2, :).^2);

    % Compute distances from points in the second image to their epipolar lines
    distances = abs(sum(l .* x2_h_n, 1));
   
end
%%

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

    fprintf('\n--- Cheirality_triangulate Validation Start ---\n');

    % Print the selected solution
    fprintf('Selected P2 Index: %d\n', P2_correct_ind);
    fprintf('Number of Valid Points: %d out of %d\n', max_points_1, length(x1_h_n_in));

    % Optional: Compute reprojection error
    err1 = compute_reprojection_error_P(P{1}, X_best, x1_h_n_in); % (:, valid_points));
    err2 = compute_reprojection_error_P(P{2}, X_best, x2_h_n_in); % (:, valid_points));
    fprintf('Mean Reprojection Error: %.4f (Image 1), %.4f (Image 2)\n', mean(err1), mean(err2));

    fprintf('\n--- Cheirality_triangulate Validation End ---\n');
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
    
    fprintf('\n--- triangulate_3D_point_DLT Validation Start ---\n');
    % Optional: Compute reprojection error
    err1 = compute_reprojection_error_P(P1, X, x1);
    err2 = compute_reprojection_error_P(P2, X, x2);
    fprintf('Mean Reprojection Error: %.4f (Image 1), %.4f (Image 2)\n', mean(err1), mean(err2));

    fprintf('\n--- triangulate_3D_point_DLT Validation End ---\n');

    % Normalize the points to make them homogeneous (4th coordinate = 1)
    X = X(1:4, :) ./ X(4, :);
end

%%

% Function: Robust Translation Estimation
function T_best = estimate_T_robust(xs, Xs, R, inlier_threshold)

    % T_init
    T_init = zeros(3,1);

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
        
        %T_candidate = compute_translation(xs_sample, Xs_sample, R);

        % Compute reprojection errors
        % % errors = compute_reprojection_error_P([R,T_candidate], Xs, xs);
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

    fprintf('Percentage of inliers in T_robust: %.2f%%\n', 100 * (max_inliers) / N);


    % Refine with inliers
    if max_inliers > 0
        T_best = T_best;
    else
        warning('No inliers found. Returning initial translation.');
        T_best = T_init;
    end

    errors = compute_reprojection_errors(xs, Xs, R, T_best);

    % Count inliers
    inliers_candidate = errors.^2 < inlier_threshold.^2;
    num_inliers = sum(inliers_candidate);

    fprintf('Percentage of inliers in T_robust best: %.2f%%\n', 100 * (num_inliers) / N);

end




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


% Function: Compute Reprojection Errors
function errors = compute_reprojection_errors(xs, Xs, R, T)
    % Project 3D points
    projected = (R * Xs(1:3, :) + T);
    projected = projected(1:2, :) ./ projected(3, :); % Normalize

    % Compute errors
    errors = sqrt(sum((xs(1:2, :) - projected).^2, 1)); % Euclidean distance
end


function errors = compute_reprojection_error_P(P, X, x)
% COMPUTE_REPROJECTION_ERROR - Compute reprojection error for 3D points
% Inputs:
%   P - Projection matrix (3x4)
%   X - 3D points (4xN homogeneous)
%   x - 2D points (3xN homogeneous)
% Outputs:
%   errors - Reprojection errors (1xN)

    % Project 3D points to 2D
    x_proj = P * X;
    x_proj = x_proj ./ x_proj(3, :); % Normalize to homogeneous

    % Compute Euclidean distance error
    errors = sqrt(sum((x_proj(1:2, :) - x(1:2, :)).^2, 1)); % Pixel errors
end


%% LM Helper Functions
function [X_refined, accumulated_error] = LevenbergMarquardt(P1, P2, X, x1, x2, mu, max_iterations, tolerance)
    % Levenberg-Marquardt optimization for 3D point refinement

    % Temporary variable for updates
    X_temp = X;
    accumulated_error = zeros(1, max_iterations); % For plotting error

    for iter = 1:max_iterations
        % disp(['Iteration Index: ', num2str(iter)]);

        % Compute total error before update
        total_error_before = 0;
        for j = 1:size(X_temp, 2)
            [err, ~] = ComputeReprojectionError(P1, P2, X_temp(:, j), x1(:, j), x2(:, j));
            total_error_before = total_error_before + err;
        end
        accumulated_error(iter) = total_error_before; % Store error for plotting
        % disp(['Total Error Before: ', num2str(total_error_before)]);

        % Update 3D points
        for j = 1:size(X_temp, 2)
            % Compute residual and Jacobian for point Xj
            [r, J] = LinearizeReprojErr(P1, P2, X_temp(:, j), x1(:, j), x2(:, j));

            % Compute LM update
            delta_Xj = ComputeUpdate(r, J, mu);

            % Update 3D point (temporarily)
            X_test = pflat(X_temp(:, j) + delta_Xj);

            % Compute individual errors
            [err_before, ~] = ComputeReprojectionError(P1, P2, X_temp(:, j), x1(:, j), x2(:, j));
            [err_after, ~] = ComputeReprojectionError(P1, P2, X_test, x1(:, j), x2(:, j));

            % Apply update if error improves
            if err_after < err_before
                X_temp(:, j) = X_test;
            end
        end

        % Compute total error after update
        total_error_after = 0;
        for j = 1:size(X_temp, 2)
            [err, ~] = ComputeReprojectionError(P1, P2, X_temp(:, j), x1(:, j), x2(:, j));
            total_error_after = total_error_after + err;
        end
        % disp(['Total Error After: ', num2str(total_error_after)]);

        % Adjust damping factor and check for convergence
        if total_error_after < total_error_before
            mu = mu / 2; % Decrease damping factor
        else
            mu = mu * 2; % Increase damping factor
        end

        if abs(total_error_after - total_error_before) < tolerance
            % disp('Converged!');
            break;
        end
    end

    % Return refined points and accumulated errors
    X_refined = X_temp;
end


function J = ComputeJacobian(P, Z)
    % Decompose rows of the camera P
    P1 = P(1, :);
    P2 = P(2, :);
    P3 = P(3, :);

    % Compute Jacobian (2x4)
    J = [
        (Z(1) * P3  - Z(3) * P1 ) / Z(3)^2;
        (Z(2) * P3  - Z(3) * P2 ) / Z(3)^2;
    ];
end


function [err, res] = ComputeReprojectionError(P1, P2, Xj, x1j, x2j)
    % Project 3D point into both cameras
    Z1 = P1 * Xj;
    Z2 = P2 * Xj;

    % Compute residuals
    r1 = x1j(1:2) - Z1(1:2) ./ Z1(3); % 2x1
    r2 = x2j(1:2) - Z2(1:2) ./ Z2(3); % 2x1
    res = [r1; r2]; % 4x1

    % Compute total reprojection error
    err = sum(res.^2); % 1x1
end


%%% Confirmed 


function total_error = ComputeTotalError(P1, P2, X, x1, x2)
    % Compute the total reprojection error for all 3D points
    %
    % Inputs:
    %   - P1, P2: Camera projection matrices
    %   - X: 4xN matrix of 3D points in homogeneous coordinates
    %   - x1, x2: 3xN matrices of 2D image points in homogeneous coordinates
    %
    % Outputs:
    %   - total_error: The total reprojection error across all points

    total_error = 0; % Initialize total error
    for j = 1:size(X, 2)
        [err, ~] = ComputeReprojectionError(P1, P2, X(:, j), x1(:, j), x2(:, j));
        total_error = total_error + err;
    end
end

function delta_Xj = ComputeUpdate(r, J, mu)
    % Compute normal equations for the update
    H = J.' * J + mu * eye(size(J, 2)); % Hessian approximation (4x4)
    g = J.' * r; % Gradient (4x1)

    % Compute the update step
    delta_Xj = -H \ g; % inv(H)
% Computes the LM update .
end


function [r, J] = LinearizeReprojErr(P1, P2, Xj, x1j, x2j)
    % Project 3D point into both cameras
    Z1 = P1 * Xj;
    Z2 = P2 * Xj;

    % Compute residuals
    r1 = x1j(1:2) - Z1(1:2) ./ Z1(3); % 2x1
    r2 = x2j(1:2) - Z2(1:2) ./ Z2(3); % 2x1
    r = [r1; r2]; % 4x1

    % Compute Jacobians using chain rule
    J1 = ComputeJacobian(P1, Z1); % 2x4
    J2 = ComputeJacobian(P2, Z2); % 2x4

    % Combine Jacobians
    J = [J1; J2]; % 4x4
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
    

    % %Label all cameras
    % for i = 1:length(P)
    %     label = sprintf('Camera %d', i);
    %     text(c(1,i), c(2,i), c(3,i), label, 'FontSize', 1, 'Color', 'b');
    % end
    
    % Enhance the plot for clarity
    xlabel('X');
    ylabel('Y');
    zlabel('Z');
    grid on;
    axis equal;
    hold off;
end
%%

function v_h = toHomogeneous(v)
% TOHOMOGENEOUS - Converts a vector to homogeneous coordinates.
%
% Usage:
%   v_h = toHomogeneous(v)
%
% Inputs:
%   v - Input vector (NxM) where N is the dimensionality, and M is the number of points.
%
% Outputs:
%   v_h - Output vector in homogeneous coordinates ((N+1)xM).
%
% Example:
%   v = [1; 2; 3]; % 3D point
%   v_h = toHomogeneous(v); % Convert to homogeneous -> [1; 2; 3; 1]
%

% Validate input
if nargin < 1
    error('Input vector is required.');
end

% Add a row of ones to make it homogeneous
v_h = [v; ones(1, size(v, 2))];

end

%%
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

%%
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

    fprintf('\n--- validate_camera Validation Start ---\n');

    
    % --- Extract Rotation (R) and Translation (T) ---
    R = P(1:3, 1:3); % Rotation matrix
    T = P(1:3, 4);   % Translation vector

    % --- Validate Orthogonality of R ---
    orthogonality_error = norm(R * R' - eye(3)); % Should be close to 0
    fprintf('R Orthogonality Error: %.6e\n', orthogonality_error);

    % --- Validate Determinant of R ---
    det_R = det(R); % Should be close to +1
    fprintf('Determinant of R: %.6f\n', det_R);

    % --- Display Translation Information ---
    fprintf('Translation Vector T: [%.4f, %.4f, %.4f]\n', T(1), T(2), T(3));
    fprintf('Norm of T: %.4f\n', norm(T));

    % --- Additional Checks ---
    % Warn if errors exceed acceptable thresholds
    if orthogonality_error > 1e-6
        warning('R is not orthogonal. Check rotation matrix!');
    end
    if abs(det_R - 1) > 1e-3
        warning('Determinant of R is not 1. Check rotation matrix!');
    end
    
    fprintf('\n--- validate_camera Validation End ---\n');

end

%%

function info(message, varargin)
% INFO - Custom function to display hierarchical information messages.
% Prints messages based on the global info_level and enableInfo settings.
%
% Usage:
%   info('Message')                        % Defaults to Level 1
%   info('Message', 1)                      % Explicit Level 1
%   info('Message with %d', 1, 10)          % Level 1 with variables
%   info('Message with %f', 2, 3.14159)     % Level 2 with variables
%
% Global Settings:
%   - enableInfo: Toggle info messages (true/false).
%   - info_level: Controls verbosity (1 = only Level 1, 2 = Level 1 and 2).

% Access global variables
global enableInfo;
global info_level;

% Set default level to 1 if not explicitly provided
if nargin < 2 || ~isnumeric(varargin{1}) % Check if level is missing
    level = 1;               % Default level
    args = varargin;         % Treat everything else as variables
else
    level = varargin{1};     % Extract level
    args = varargin(2:end);  % Remaining are variables
end

% Print messages based on levels
if enableInfo
    % Print if the current level is less than or equal to the global level
    if level <= info_level
        fprintf('[INFO][Level %d]: ', level);
        fprintf(message, args{:}); % Pass additional arguments
        fprintf('\n');
    end
end
end

%%
function X_clean = removePointsExcessivelyFarAway(X_dirty)
    % remove points behind the camera
    X_in_front = X_dirty(:,find(X_dirty(3,:)>0));

    % remove points far away from center of gravity
    cent = mean(X_in_front(1:3,:), 2);
    dist = sqrt(sum((X_in_front(1:3,:) - cent).^2, 1));
    threshold = 4 * quantile(dist, 0.9);
    ind_clean = dist < threshold;
    X_clean = X_in_front(:, ind_clean);

end

%%
function P_un = CalibratedToUncalibrated(P_ca, K)

N = length(P_ca);

for i=1:N
    P_un{i} = K * P_ca{i};
end

end


% function T = compute_translation(xs, Xs, R)    
%     x1 = xs(:,1);
%     x2 = xs(:,2);
% 
%     X1 = Xs(:,1);
%     X2 = Xs(:,2);
% 
%     x1_cross = [0, -x1(3), x1(2); x1(3), 0, -x1(1); -x1(2), x1(1), 0];
%     x2_cross = [0, -x2(3), x2(2); x2(3), 0, -x2(1); -x2(2), x2(1), 0];
% 
%     M = [x1_cross*R*X1(1:3), x1_cross; x2_cross*R*X2(1:3), x2_cross];
% 
%     [~, ~, V] = svd(M);
%     T = V(2:end, end) ./ V(1,4);
% end