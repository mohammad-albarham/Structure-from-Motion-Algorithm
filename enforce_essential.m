function E = enforce_essential(E_approx)
% ENFORCE_ESSENTIAL Apply the enforce on E

% Apply the SVD on E approximate 
[U,~,V] = svd( E_approx );

% if det(U*V') < 0 
%     V = -V; 
% end

% Check the rank of F
rank_E_approx = rank(E_approx);

% Display the rank
% disp(['Rank of E_approx: ', num2str(rank_E_approx)]);

info("Rank of E_approx: %.2f \n:  ", 2, rank_E_approx);


% Creates a valid essential matrix from an approximate solution 

E = U * diag ([1 1 0])* V';

% Check the rank of F
rank_E = rank(E);

% Display the rank
% disp(['Rank of E: ', num2str(rank_E)]);
info("Rank of E: %.2f \n:  ", 2, rank_E);


end