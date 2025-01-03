function P = estimate_camera_DLT(x,X)
% Calculate the camera matrix using DLT
% Inputs: x, X: 3x1

% Add a dimension to X to make it homogenous 
X = [X; ones(1,length(X))]; % 4x1 

% Define the matrix M 
% Number of points 
N = size(X, 2);
M = zeros(3*N,12+N); % Initialize the matrix 

% Fill the values as blocks 
for i= 1:N
    % First row in the block
    M(3*i-2, 1:4) = X(:,i)';
    M(3*i-2, 12+i) = -x(1,i);

    % Second row in the block
    M(3*i-1, 5:8) = X(:,i)';
    M(3*i-1, 12+i) = -x(2,i);

    % Third row in the block
    M(3*i, 9:12) = X(:,i)';
    M(3*i, 12+i) = -x(3,i);
    
end

% Calculate the S, V, U of Singular value decomposition
[U,S,V] = svd(M);

% Check the smallest singular value
diagonal_numbers = diag(S); % Smallest value in the diagonal of S

smallest_singular_value = diagonal_numbers(end);

disp('Smallest Singular Value:');
disp(smallest_singular_value);

% v is the last column of V
v = V(:, end); % N+12
P = reshape(v(1:12), [4 3])'; % P is 3x4

% Check ||Mv||
Mv_abs = norm(M * v);
disp('||Mv||:');
disp(Mv_abs); % 0.0150653751103407


% Check the determinant of the rotation matrix
A = P(1:3, 1:3);

if det(A)<0
    P = -P;
end


end