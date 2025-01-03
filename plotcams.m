function plotcams(P)
    % Function to plot camera positions and viewing directions
    % Inputs:
    %   P - Cell array of camera projection matrices
    
    % Initialize arrays to store camera centers and directions
    c = zeros(4, length(P)); % Camera centers
    v = zeros(3, length(P)); % Viewing directions
    
    for i = 1:length(P)
        c(:,i) = null(P{i});          % Compute camera center (homogeneous coordinates)
        v(:,i) = P{i}(3, 1:3);        % Viewing direction (3rd row of P)
    end
    
    % Normalize homogeneous coordinates
    c = c ./ repmat(c(4,:), [4 1]);   % Convert to 3D by dividing by the last element

    % Plot cameras using quiver3
    quiver3(c(1,:), c(2,:), c(3,:), v(1,:), v(2,:), v(3,:), 'r-'); 
    hold on;
    

    %Label all cameras
    for i = 1:length(P)
        label = sprintf('Camera %d', i);
        text(c(1,i), c(2,i), c(3,i), label, 'FontSize', 10, 'Color', 'b');
    end
    
    % Enhance the plot for clarity
    xlabel('X');
    ylabel('Y');
    zlabel('Z');
    grid on;
    axis equal;
    hold off;
end