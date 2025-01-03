function E = enforce_essential(E_approx)
% ENFORCE_ESSENTIAL Apply the enforce on E

% Apply the SVD on E approximate 
[U,S,V] = svd( E_approx );

if det(U*V') < 0 
    V = -V; 
end

% Creates a valid essential matrix from an approximate solution 

E = U * diag ([1 1 0])* V';

end