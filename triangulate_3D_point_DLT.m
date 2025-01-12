function X = triangulate_3D_point_DLT(x1, x2, P1, P2)
    % INPUT:
    % x1 - 2D points in the first image (3xN homogeneous coordinates)
    % x2 - 2D points in the second image (3xN homogeneous coordinates)
    % P1 - Camera projection matrix for the first image (3x4)
    % P2 - Camera projection matrix for the second image (3x4)
    % OUTPUT:
    % X - Triangulated 3D points (4xN homogeneous coordinates)

    % Number of points given in x1 and x2 
    N = size(x1, 2);

    X = zeros(4, N); 
    
    % Loop through all points to get all 3D points
    for i = 1:N
        % Construct the A matrix for triangulation
        A = [
            x1(1,i) * P1(3,:) - P1(1,:);
            x1(2,i) * P1(3,:) - P1(2,:);

            x2(1,i) * P2(3,:) - P2(1,:);
            x2(2,i) * P2(3,:) - P2(2,:);
        ];
        
        % Solve the homogeneous system using SVD as before and we need only
        % V to extract v 

        [~, ~, V] = svd(A);
        X(:,i) = V(:,end); % Last column of V gives the solution ( last eigenvector)
    end
    
    info("\n--- triangulate_3D_point_DLT Validation Start ---\n");
    % Optional: Compute reprojection error
    err1 = compute_reprojection_error_P(P1, X, x1);
    err2 = compute_reprojection_error_P(P2, X, x2);

    info('Mean Reprojection Error: %.4f (Image 1), %.4f (Image 2)\n', 1, mean(err1), mean(err2));

    info('\n--- triangulate_3D_point_DLT Validation End ---\n');

    % Normalize the points to make them homogeneous (4th coordinate = 1)
    X = X(1:4, :) ./ X(4, :);
end