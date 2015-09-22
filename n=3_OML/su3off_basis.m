% This function returns the basis for su(3), minus its diagonal subspace. There
% are six such basis elements

function[basis] = su3off_basis()

basis = zeros(3,3,6);

% 1-2 elements
basis(1,2,1) = 1i;
basis(2,1,1) = 1i;
basis(1,2,2) = 1;
basis(2,1,2) = -1;

% 2-3 elements
basis(2,3,3) = 1i;
basis(3,2,3) = 1i;
basis(2,3,4) = 1;
basis(3,2,4) = -1;

% 1-3 elements
basis(3,1,5) = 1i;
basis(1,3,5) = 1i;
basis(3,1,6) = 1;
basis(1,3,6) = -1;

end

