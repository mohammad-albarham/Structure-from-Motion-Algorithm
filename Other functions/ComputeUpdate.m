function delta_Xj = ComputeUpdate(r, J, mu)
    % Compute normal equations for the update
    H = J.' * J + mu * eye(size(J, 2)); % Hessian approximation (4x4)
    g = J.' * r; % Gradient (4x1)

    % Compute the update step
    delta_Xj = -H \ g; % inv(H)
% Computes the LM update .
end