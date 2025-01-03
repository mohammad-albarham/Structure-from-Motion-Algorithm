% reset and load data
clc; close all; clear;
% run('C:\Users\pain\Desktop\CV_project\vlfeat-0.9.21-bin\vlfeat-0.9.21\toolbox\vl_setup')
% Get dataset info by providing dataset input

% Control random number generator - to get consistant resutls

rng(42)

% Pre-settings 
output_folder = 'plot_folder';

% Step 2: Create the folder if it does not exist
if ~exist(output_folder, 'dir')
    mkdir(output_folder);
end

%% Setp 1: Extract dataset info

[K, img_names, init_pair, pixel_threshold] = get_dataset_info(2);

%% 

% Loop over all images 

img_n = length(img_names);

for n=1:img_n-1
    
    % Step 2: Compute the features using SIFT, extract the matches and the corresponding points
    
    % Load the images 
    im1 = imread(img_names{n});
    im2 = imread(img_names{n+1});
    
    % Extract the features using SIFT
    [fA dA] = vl_sift(single(rgb2gray(im1)) );
    [fB dB] = vl_sift(single(rgb2gray(im2)) );
    
    % Compute the matches
    matches = vl_ubcmatch(dA ,dB );
    
    % Get the 2D points correspondences on each images' pair
    xA = fA(1:2 , matches(1 ,:));
    xB = fB(1:2 , matches(2 ,:));
    
    % Convert the vectors to be homogenous for this pair of images 
    
    % Add a homogenous dimension to the points
    x1 = [xA; ones(1,length(xA))];
    x2 = [xB; ones(1,length(xB))];
    
    % Normalized points using K (camera intrinsic matrix)
    x1_norm = inv(K) * x1; % Points in the first image
    x2_norm = inv(K) * x2; % Points in the second image
    
    % Step3: Estimate E robust using RANSAC
    
    % Define the pixel threhold 
    inlier_threshold_px = pixel_threshold;
    eps = inlier_threshold_px * 2 / (K(1,1) + K(2,2)); % Normalize threshold
    
    % Get the estimate_E_robust
    [E_ransac, inliers] = estimate_E_robust(x1, x2, eps, K);
    
    % Display results
    disp('Estimated Essential Matrix (RANSAC):');
    disp(E_ransac);
    disp('Number of inliers:');
    disp(sum(inliers));
    
    
    % Step 4: Triangulate the points and Check for Cheirality Condition using Triangulation:
    
    % Include only the inliers on the triangulation
    x1_inliers = x1(:,inliers);
    x2_inliers = x2(:,inliers);
    
    % Correct camera P2
    P2_correct_ind = 0 ;
    max_points_1 = 0;
    max_points_2 = 0;
    
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
    
    correct_index = -1; 
    
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
    
    X = X_best;
    
    % Plot the best camera with best 3D points
    
    figure();
    plot3(X(1,:), X(2,:), X(3,:), '.', 'MarkerSize', 10);
    hold on;
    P{1} = P1; 
    P{2} = P2_candidates{P2_correct_ind};
    plotcams(P);
    
    xlabel('x'); ylabel('y'); zlabel('z');
    title("Best camera, 3D Reconstruction")
    axis equal;
    hold off;
    
    
    % Get the robust T using DLT
    
    % disply T before 
    disp('Translation T before:');
    disp(inv(K) * P{2}(1:end,end));
    
    inlier_threshold = 1;
    
    % Use inlier correspondences for 2D points
    xn = inv(K) * x2_inliers;
    xsn = xn(1:2, :);  % 2D points (Nx2)
    
    % Use triangulated 3D points
    Xs = inv(K) * X(1:3, :);           % 3D points (Nx3)
    Xsn = Xs;
    
    R = inv(K) * P{2}(1:3, 1:3); % Extract rotation from P2
    
    translation_threshold = 3 * pixel_threshold / K(1,1);
    
    % The normalized coordinates is used here 
        T_best = estimate_T_robust(xsn, Xsn, R, translation_threshold, inv(K) * P{2}(1:end,end))
        
        disp('Estimated Robust Translation T:');
        disp(T_best);

end