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
