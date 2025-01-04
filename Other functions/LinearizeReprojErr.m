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