function distances = compute_epipolar_errors(F, x1_h_n, x2_h_n)
    % Compute epipolar errors for a given fundamental matrix F and point correspondences
    %
    % Inputs:
    %   F - Fundamental matrix (3x3)
    %   x1_h_n - Points in the first image (3xN, homogeneous coordinates) and
    %   normalized
    %   x2_h_n - Points in the second image (3xN, homogeneous coordinates) and
    %   normalized
    %
    % Output:
    %   distances - Vector of epipolar distances for each point correspondence
    % Compute epipolar lines for points in the second image
    l = F * x1_h_n; % Epipolar lines in the second image

    % Normalize the epipolar lines to have unit norm
    l = l ./ sqrt(l(1, :).^2 + l(2, :).^2);

    % Compute distances from points in the second image to their epipolar lines
    distances = abs(sum(l .* x2_h_n, 1));
   
end