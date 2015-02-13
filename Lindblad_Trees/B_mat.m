% The eigenvalues of a density matrix under Lindblad dynamics obey the ODE:
% 	d(lambda)/dt = B*lambda
% This function computes the matrix B from the Lindblad rates w(i,j)

function [B] = B_mat(w)

	% The rates w(i,j) are the off-diagonal elements. The diagonal
	% elements are such that each column sums to zero (so that the
	% sum of eigenvalues is stationary at one.
	
	B = w - diag(sum(w));

end
