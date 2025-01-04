function [err, res] = ComputeReprojectionError(P1, P2, Xj, x1j, x2j)
    % Project 3D point into both cameras
    Z1 = P1 * Xj;
    Z2 = P2 * Xj;

    % Compute residuals
    r1 = x1j(1:2) - Z1(1:2) ./ Z1(3); % 2x1
    r2 = x2j(1:2) - Z2(1:2) ./ Z2(3); % 2x1
    res = [r1; r2]; % 4x1

    % Compute total reprojection error
    err = sum(res.^2); % 1x1
end


%%% Confirmed 


