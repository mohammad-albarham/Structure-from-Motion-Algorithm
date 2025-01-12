function [] = run_sfm(dataset_num, varargin)
% RUN_SFM - Structure from Motion (SfM) Pipeline
%
% This function processes a dataset for Structure from Motion, with optional
% visualization and logging control through key-value parameters.
%
% Inputs:
%   dataset_num - (Required) Dataset number to process (integer).
%   varargin    - (Optional) Key-value pairs:
%       'plot_initial_rel_3D' (default: false) - Enable plotting of initial
%           relative 3D points for all pairs.
%       'plot_initial_rel_3D_pair' (default: false) - Enable plotting of initial
%           relative 3D points for the first pair.
%       'enableInfo' (default: true) - Enable or disable informational messages.
%       'info_level' (default: 1) - Set informational detail level:
%           1 - Basic informational messages.
%           2 - Detailed execution logs.
%
% Example Usage:
%   run_sfm(2); % Default settings
%   run_sfm(2, 'plot_initial_rel_3D', true); % Enable 3D plotting for all pairs
%   run_sfm(2, 'enableInfo', false); % Disable informational messages
%   run_sfm(2, 'info_level', 2); % Enable detailed execution logs
%
% Notes:
%   - This function uses global variables 'enableInfo' and 'info_level'.
%   - 'info_level' accepts values of 1 or 2 only.

% --- Parse Optional Arguments ---
p = inputParser; % Create input parser
addRequired(p, 'dataset_num', @(x) isnumeric(x) && isscalar(x)); % Dataset number
addParameter(p, 'plot_initial_rel_3D', false, @(x) islogical(x)); % Plotting flag
addParameter(p, 'plot_initial_rel_3D_pair', false, @(x) islogical(x)); % Initial pair plotting flag
addParameter(p, 'enableInfo', false, @(x) islogical(x)); % Enable or disable logs
addParameter(p, 'info_level', 1, @(x) isnumeric(x) && ismember(x, [1, 2])); % Log detail level (1 or 2)

parse(p, dataset_num, varargin{:}); % Parse inputs

% --- Set Global Variables ---
global enableInfo;
global info_level;


% --- Extract Parameters ---
dataset_num = p.Results.dataset_num;
plot_initial_rel_3D = p.Results.plot_initial_rel_3D;
plot_initial_rel_3D_pair = p.Results.plot_initial_rel_3D_pair;
enableInfo = p.Results.enableInfo;
info_level = p.Results.info_level;


% enableInfo = true; % Enable or disable info globally
% info_level = 1;

% reset and load data

% Add the vl_toolbox to be accessible in the current directory
run('.\vlfeat-0.9.21-bin\vlfeat-0.9.21\toolbox\vl_setup')
% Get dataset info by providing dataset input

% Control random number generator - to get consistant resutls
rng(42)

% Load the dataset data information
info("Step 0: Load the dataset data information")

% Extract dataset information
[K, img_names, init_pair, pixel_threshold] = get_dataset_info(dataset_num);

% Step 1: Calculate the relative Orientation
% This step calculates the relative rotations (R_rel)
% between consecutive image pairs in the dataset.
% It also performs feature extraction, matches features, estimates the essential
% matrix using RANSAC, and triangulates 3D points.

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

    validate_camera(P_pair{1});
    validate_camera(P_pair{2});

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


% Step 2: Upgrade to Absolute Rotations

info("Step 2: Upgrade to Absolute Rotations")

% Initialize cell array for storing absolute rotations
R_abs_i = cell(1, N);            % Absolute rotations for N cameras
R_abs_i{1} = eye(3);             % First camera has identity rotation (reference)

% Compute absolute rotations from relative rotations
for k = 1:N-1
    % Compute absolute rotation by chaining relative rotations
    R_abs_i{k+1} = R_rel{k} * R_abs_i{k}; % Accumulate rotations

    % Display progress for debugging
    info("Computed Absolute Rotation (Camera %d):\n",1,k+1);
    info("Absolute camera",1,R_abs_i{k+1});
end

% Final output
info("Absolute Rotations Computed for all Cameras!");


% Verify orthogonality (R * R' = I)
info("\nVerifying Orthogonality of Rotations:\n");

for k = 1:N
    error_orthogonality = norm(R_abs_i{k} * R_abs_i{k}' - eye(3)); % Should be close to 0 using Euclidean norm
    info("Camera %d Orthogonality Error: %.6f\n", 1,k, error_orthogonality);
end

% Step 3: Reconstruct initial 3D points from initial image pair

info("Step 3: Reconstruct initial 3D points from initial image pair:\n");

% Get the 3D points from the suggested pair

% Compute the matches
matches = vl_ubcmatch(desc_im{init_pair(1)}, desc_im{init_pair(2)});

% Get the 2D points correspondences on each images' pair
x1 = feat_im{init_pair(1)}(1:2 , matches(1 ,:));
x2 = feat_im{init_pair(2)}(1:2 , matches(2 ,:));

% Save descriptors for matched points (3D points will correspond to these)
% Extract descriptors only for matched points in image 1
desc_X = desc_im{init_pair(1)}(:, matches(1, :)); 

% vectors with dim as
% 3 x N;
% N-> Number of the points
x1_h = toHomogeneous(x1);
x2_h = toHomogeneous(x2);

x1_h_n = K\x1_h;
x2_h_n = K\x2_h;

% Get the estimate_E_robust
[E_ransac, indices] = estimate_E_robust(x1_h_n, x2_h_n, epipolar_threshold);

% Triangulate points using relative pose

x1_h_n_in = x1_h_n(:,indices);
x2_h_n_in = x2_h_n(:,indices);

[X_init_pair,P_init_pair] = Cheirality_triangulate(x1_h_n_in,x2_h_n_in,E_ransac);

validate_camera(P_init_pair{1})
validate_camera(P_init_pair{2})

% % Save descriptors for future use "Only inliers"
desc_X_init_pair = desc_X(:,indices);

% Move the 3D pionts to the 3D world coordinate
X_init_pair_wc = R_abs_i{init_pair(1)}' * X_init_pair(1:3, :);
X_init_pair_wc_h = toHomogeneous(X_init_pair_wc);


if plot_initial_rel_3D_pair
    % Step 5: Visualization for initial pair
    figure();
    plot3(X_init_pair_wc_h(1, :), X_init_pair_wc_h(2, :), X_init_pair_wc_h(3, :), ".");
    hold on;
    xlabel('x'); ylabel('y'); zlabel('z');
    title('Filtered 3D Points');
    hold off;
end

% Step 4:
info("Step 4: T robust estimation:\n");

% Establish correspondences between i and 3D points (using desc X),
% Number of images
% Initialize translation estimates
T_best = cell(1, N);
% T_best{1} = zeros(3, 1);

for n = 1:N
    % Step 1: Establish correspondences using descriptor matching
    [matches_2d_3d, ~] = vl_ubcmatch(desc_X_init_pair, desc_im{n});

    % Step 2: Extract 2D points and convert to homogeneous coordinates
    xA = feat_im{n}(1:2, matches_2d_3d(2, :));

    x1_h = toHomogeneous(xA);

    x1_h_n = K \ x1_h; % Normalize 2D points

    % Step 3: Extract 3D points and convert to homogeneous
    Xs_h = X_init_pair_wc_h(:, matches_2d_3d(1, :)); % Use pre-computed 3D points

    % Step 4: Robust estimation of translation
    translation_threshold = 3 * pixel_threshold / K(1, 1);

    info("Determinant of R: ",2,det(R_abs_i{n}));

    % Estimate T using Robust RANSAC + DLT
    T_best{n} = estimate_T_robust(x1_h_n, Xs_h, R_abs_i{n}, translation_threshold);
end

% Step 5: Plot the cameras:

info("Step 5: Plot the cameras: \n:  ");

P_all = cell(1,N);

for i=1:N
    P_all{i} = [R_abs_i{i} T_best{i}];
end

figure();

P_all_un = CalibratedToUncalibrated(P_all, K);
plotcams(P_all_un);

title('Camera Poses Visualization');

info("Cameras plotted successfully!  ");

% Step 6: Triangulate points for all pairs (i, i + 1) and visualize 3D points + cameras
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

    % Keep good points
    X = keep_good_points(points_in_front(X));

    % Save points and corresponding color in the cell array
    points_with_colors{n} = struct('points', X, 'color', colormap(n, :));
end

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

info("The SFM Algorithm is Finished");

end