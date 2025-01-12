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