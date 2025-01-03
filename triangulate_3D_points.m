function [X] = triangulate_3D_points(x1,x2, inliers, P1, P2)

  % Include only the inliers on the triangulation
        x1_inliers = x1(:,inliers);
        x2_inliers = x2(:,inliers);
        
        % Correct camera P2
        P2_correct_ind = 0 ;
        max_points_1 = 0;
        max_points_2 = 0;
        
        % X best 
        X_best = [];
       
        X = triangulate_3D_point_DLT(x1_inliers, x2_inliers,P1, P2);
        x1_proj = P1 * X; 
        x2_proj = P2 * X; 
    
        % This is a general condition for the position of the points 
        cond = sum((x1_proj(3,:) > 0) & ((x2_proj(3,:) > 0)) );
       
        disp("Points in front of the camera")
        disp(cond)
        
        % Plot the best camera with best 3D points
        % figure();
        % plot3(X(1,:), X(2,:), X(3,:), '.', 'MarkerSize', 10);
        % hold on;
        % P{1} = P1; 
        % P{2} = P2;
        % plotcams(P);
        % 
        % xlabel('x'); ylabel('y'); zlabel('z');
        % title("Initial pair, 3D Reconstruction")
        % axis equal;
        % hold off;
end

