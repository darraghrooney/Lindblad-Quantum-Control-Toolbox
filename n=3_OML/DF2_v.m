% This function computes the second flag derivative of the tangent vector b - A*x
% for n=3, given a collection of Lindblad operators, a flag, and a position x.

function[d2v] = DF2_v(Ls, flag, x)

% Fetch the second derivative of the Lindblad rates
d2w = DF2_w(Ls,flag);

% Initialize
d2v = zeros(2,6,6);

% Sweep over pairs of basis components
for j = 1:6
	for k = 1:6

		% Put the 2nd deriv of w's into square form
		dw2_square = zeros(3,3);
		dw2_square(1,2) = d2w(j,k,1);
		dw2_square(1,3) = d2w(j,k,2);
		dw2_square(2,1) = d2w(j,k,3);
		dw2_square(2,3) = d2w(j,k,4);
		dw2_square(3,1) = d2w(j,k,5);
		dw2_square(3,2) = d2w(j,k,6);
	
		% Fetch the 2nd deriv of A and b
		[d2A,d2b] = A_b_rinf(dw2_square);

		% Calculate 2nd deriv of tangent vector
		d2v(:,j,k) = d2b - d2A*x;
	end
end

end
