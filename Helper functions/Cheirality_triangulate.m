function [X_best,P] = Cheirality_triangulate(x1,x2, inliers,K,E_ransac)

  % Include only the inliers on the triangulation
        x1_inliers = x1(:,inliers);
        x2_inliers = x2(:,inliers);
        
        % Correct camera P2
        P2_correct_ind = 0 ;
        max_points_1 = 0;
        
        % X best 
        X_best = [];
        
        % Define the cameras & assume first camera is at the origin
        P1 = K* [eye(3), zeros(3, 1)];
        
        % Get the U and V canditates
        [U, ~, V] = svd(E_ransac);
        W = [0 -1 0; 1 0 0; 0 0 1]; % Rotation matrix for decomposition
        
        % Four possible solutions for P2
        R1 = U * W * V';
        R2 = U * W' * V';
        t = U(:, 3);
        
        % Ensure proper rotations (determinant = +1)
        if det(R1) < 0, R1 = -R1; end
        if det(R2) < 0, R2 = -R2; end
        
        % Four possible camera matrices
        P2_candidates = {
            K * [R1, t],
            K * [R1, -t],
            K * [R2, t],
            K * [R2, -t]
        };
                
        % Choose the correct P2 by testing chirality (points must be in front of both cameras)
        for i = 1:4
            P2 = P2_candidates{i};
        
            X = triangulate_3D_point_DLT(x1_inliers, x2_inliers,P1, P2);
            x1_proj = P1 * X; 
            x2_proj = P2 * X; 
        
            % This is a general condition for the position of the points 
            cond = sum( (x1_proj(3,:) > 0) & ((x2_proj(3,:) > 0)) );
            
            if cond > max_points_1
                max_points_1 = cond;
                P2_correct_ind = i;
                X_best = X;
            end
        end
        
        P{1} = P1; 
        P{2} = P2_candidates{P2_correct_ind};
end

