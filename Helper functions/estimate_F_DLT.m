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
