% This function computes the flag derivative of A(w) and b(w) for n=3, given
% a collection of Lindblad operators and a flag.

function[dA, db] = DF_Ab(Ls, flag)

% Compute the flag derivative of the wij's
dw = DF_w(Ls,flag);

% Initialize
db = zeros(2,6);
dA = zeros(2,2,6);

% Sweep over flag-tangent basis
for j = 1:6
	% Form a square matrix out of the six DF_w's 
	dw_square = zeros(3,3);
	dw_square(1,2) = dw(1,j);
	dw_square(1,3) = dw(2,j);
	dw_square(2,1) = dw(3,j);
	dw_square(2,3) = dw(4,j);
	dw_square(3,1) = dw(5,j);
	dw_square(3,2) = dw(6,j);		

	% Compute dA and db for this basis component
	[dAj,dbj,trash] = A_b_rinf(dw_square);
	dA(:,:,j) = dAj;
	db(:,j) = dbj;
end

end
