% Define 3D points in the world coordinate system
n_points = 10; % Number of points
P_world = rand(3, n_points) * 10; % Random 3D points

% Camera 1 (reference camera)
R1 = eye(3); % Identity rotation
t1 = [0; 0; 0]; % Zero translation

% Camera 2 (freely chosen pose)
theta = pi/4; % 45 degrees rotation around z-axis
R2 = [cos(theta), -sin(theta), 0; sin(theta), cos(theta), 0; 0, 0, 1]; % Rotation
t2 = [5; 2; 3]; % Translation

% Compute relative rotation and translation
R12 = R2; % Relative rotation is the rotation of Camera 2
t12 = -R2 * t2; % Relative translation

% Transform points
P_cam1 = R1 * (P_world - t1); % Points in Camera 1 frame (same as P_world)
P_cam2 = R2 * (P_world - t2); % Points in Camera 2 frame

% Validate transformation
P_cam2_reconstructed = R12 * P_cam1 + t12; % Reconstructed points in Camera 2
error = norm(P_cam2 - P_cam2_reconstructed, 'fro'); % Validation error

% Display results
disp('Relative Rotation Matrix (R12):');
disp(R12);
disp('Relative Translation Vector (t12):');
disp(t12);
disp('Reconstruction Error:');
disp(error);

% Visualization
figure;
plot3(P_world(1, :), P_world(2, :), P_world(3, :), 'ko'); hold on;
plot3(P_cam1(1, :), P_cam1(2, :), P_cam1(3, :), 'r*');
plot3(P_cam2(1, :), P_cam2(2, :), P_cam2(3, :), 'b+');
legend('World Points', 'Camera 1 Points', 'Camera 2 Points');
xlabel('X'); ylabel('Y'); zlabel('Z');
title('3D Points in World, Camera 1, and Camera 2 Frames');
grid on;