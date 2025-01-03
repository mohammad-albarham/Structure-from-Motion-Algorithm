function distances = compute_epipolar_errors(F, x1s, x2s)
    % Compute epipolar errors for a given fundamental matrix F and point correspondences
    %
    % Inputs:
    %   F - Fundamental matrix (3x3)
    %   x1s - Points in the first image (3xN, homogeneous coordinates)
    %   x2s - Points in the second image (3xN, homogeneous coordinates)
    %
    % Output:
    %   distances - Vector of epipolar distances for each point correspondence
    % Compute epipolar lines for points in the second image
    l = F * x1s; % Epipolar lines in the second image

    % Normalize the epipolar lines to have unit norm
    l = l ./ sqrt(l(1, :).^2 + l(2, :).^2);

    % Compute distances from points in the second image to their epipolar lines
    distances = abs(sum(l .* x2s, 1));
   
end

