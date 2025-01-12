function X_front= points_in_front(X_points)
    % remove points behind the camera
    X_front = X_points(:,X_points(3,:)>0);
end 