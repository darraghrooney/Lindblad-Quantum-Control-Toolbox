% This function returns the generalized Pauli matrices, which form a basis
% on su(3)

function[gP] = gPaulis()

gP = zeros(3,3,8);

gP(:,:,1) = [[0,1,0];[1,0,0];[0,0,0]]/sqrt(2);
gP(:,:,2) = [[0,-1i,0];[1i,0,0];[0,0,0]]/sqrt(2);
gP(:,:,3) = [[1,0,0];[0,-1,0];[0,0,0]]/sqrt(2);
gP(:,:,4) = [[0,0,1];[0,0,0];[1,0,0]]/sqrt(2);
gP(:,:,5) = [[0,0,-1i];[0,0,0];[1i,0,0]]/sqrt(2);
gP(:,:,6) = [[0,0,0];[0,0,1];[0,1,0]]/sqrt(2);
gP(:,:,7) = [[0,0,0];[0,0,-1i];[0,1i,0]]/sqrt(2);
gP(:,:,8) = [[1,0,0];[0,1,0];[0,0,-2]]/sqrt(6);


end
