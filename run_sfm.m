
% RUN_SFM - Main script for Structure from Motion (SfM)
% Usage: run_sfm(<dataset number>)

% reset and load data
% clc; close all;
run('C:\Users\pain\Desktop\CV_project\vlfeat-0.9.21-bin\vlfeat-0.9.21\toolbox\vl_setup')
% Get dataset info by providing dataset input

% Control random number generator - to get consistant resutls
rng(42)

%%
% Extract dataset info
[K, img_names, init_pair, pixel_threshold] = get_dataset_info(2);

%%
% Step 1: Calculate the relative Orientation

disp("Step 1: Calculate the relative Orientation")

% Number of images
N = length(img_names);

% Initialize rotations and translations for relative poses
R_rel = cell(1, N-1); % Relative rotations, normalised
T_rel = cell(1, N-1); % Relative translations, normalised

desc_im = cell(1,N); 
feat_im = cell(1,N);

inliers_indices = cell(1,N-1);

for n=1:N-1
    % Step 2: Compute the features using SIFT, extract the matches and the corresponding points
    % Load the images 
    im1 = imread(img_names{n});
    im2 = imread(img_names{n+1});
    
    % Pipeline for extracting the featuers 
    [x1, x2, ~, fA, dA, fB, dB] = feature_extraction(im1,im2);
    
    feat_im{n} = fA; 
    desc_im{n} = dA;

    feat_im{n+1} = fB; 
    desc_im{n+1} = dB;

    % Step3: Estimate E robust using RANSAC
    % Define the pixel threhold 
    eps = pixel_threshold * 2 / (K(1,1) + K(2,2)); % Normalize threshold

    % e.g: 0.000420430520853354
    
    % Get the estimate_E_robust
    [E_ransac, indices] = estimate_E_robust(x1, x2, eps, K);
    inliers_indices{n} = indices;

    % Step 4: Triangulate the points and Check for Cheirality Condition using Triangulation:
    [X,P] = Cheirality_triangulate(x1,x2, indices,K,E_ransac);
     
    % Save relative pose
    % Extract rotation from P2
    R_rel{n} = K\P{2}(1:3, 1:3); 

    % Extract translation vector
    T_rel{n} = K\P{2}(1:3, 4); % Extract relative translation

    % Step 5: Visualization
    figure;
    plot3(X(1, :), X(2, :), X(3, :), '.', 'MarkerSize', 10);
    
    hold on;
    plotcams(P); % Plot cameras
    xlabel('x'); ylabel('y'); zlabel('z');
    title('Relative plot of 3D points');
    axis equal;
    grid on;
    hold off;

end

%%
% Step 2: Upgrade to absolute Rotations
disp("Step 2: Upgrade to absolute Rotations")


tot = length(R_rel);

% P{1} = K* [eye(3), zeros(3, 1)];

R_abs_i = cell(1, N);
R_abs_i{1} = eye(3,3);

for k=1:tot
    % Calculate the rotation translation 
    R_abs_i{k+1} = R_rel{1,k} * R_abs_i{k};
    fprintf('Absolute Rotation R_abs{%d}:\n', i);
    disp(R_abs_i{k});
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
save('desc_X.mat', 'desc_X'); % Save to a mat file

% Triangulate points using relative pose

[X_0,P] = Cheirality_triangulate(x1,x2, indices,K,E_ransac);


P_X_0 = P{2};


% Rotate the initial 3D points to the world coordinate
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

% Step 5: Visualization
figure;
plot3(X_wc_filtered(1, :), X_wc_filtered(2, :), X_wc_filtered(3, :), '.', 'MarkerSize', 10);
hold on;
plotcams(P); % Plot cameras
xlabel('x'); ylabel('y'); zlabel('z');
title('Filtered 3D Points and Cameras');
axis equal;
grid on;
hold off;



% ------------------------------------------------------------------------
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

    [unique_indices, ia] = unique(matches_2d_3d(2, :), 'stable');
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
%%
% Optional: Step 5: Optimize the translation vectors

% % Optimization parameters
% mu = 1e-3;
% max_iterations = 100;
% tolerance = 1e-6;

% % Optimize translations
% [T1_refined, T2_refined, error] = optimize_translations(K, R_abs_i{1}, R_abs_i{2}, T_best{1}, T_best{2}, [X_wc; ones(1,length(X_wc))], x1, x2, mu, max_iterations, tolerance);


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

    % Normalize the points
    x1_n = K\x1;
    x2_n = K\x2;

    % Keep only the inliers 
    x1_n_in = x1_n(:, inliers_indices{n});
    x2_n_in = x2_n(:, inliers_indices{n});

     % Triangulate the 3D points using the cameras and the 2D points 

    X = triangulate_3D_point_DLT(x1_n_in, x2_n_in, inv(K)* P_all{i},  inv(K) * P_all{i+1});
    % X = K * X(1:3,:);
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

end

%%