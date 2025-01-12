function [V] = pflat(V)
% Normalize the n points in P^(m-1)

V = V(:,:) ./ V(end, :);

end