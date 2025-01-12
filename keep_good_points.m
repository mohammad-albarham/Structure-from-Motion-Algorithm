function X_good = keep_good_points(X_front)
    % remove points far away from center of gravity
    center_points = mean(X_front(1:3,:), 2);
    distance_from_center = sqrt(sum((X_front(1:3,:) - center_points).^2, 1));
    threshold = 4 * quantile(distance_from_center, 0.9);
    ind_clean = distance_from_center < threshold;
    X_good = X_front(:, ind_clean);
end
