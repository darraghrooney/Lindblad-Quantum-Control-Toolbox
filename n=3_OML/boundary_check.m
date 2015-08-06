% This function checks that the boundary of B03 is closed. If it is not, it returns
% is_closed = 0, a vector of ints between 1 and 6 that tell you which vector fields
% point outwards, and theta, between 0 and pi/3, which tells you the angular position
% on the sector I boundary where closure fails. 

function[is_closed,lc,theta] = boundary_check(w, boundaryE)

% Find the six b-vectors and A-matrices

[A1,b1,rinf] = A_b_rinf(w);
S3g = S3gen;
As = zeros(2,2,6);
bs = zeros(2,6);
for j = 1:6
	As(:,:,j) = S3g(:,:,j)*A1*S3g(:,:,j)';
	bs(:,j) = S3g(:,:,j)*b1;
end

% Initialization. Presumes the boundary is closed.

is_closed = 1;
lc = [];
theta = NaN;
grid_size = size(boundaryE,2)/6-1;

% Sweep over boundary points in sector I

for j = 1:grid_size+1

	% Get info about boundary position and the current vector field
	x = boundaryE(1:2,j);
	k = boundaryE(3,j);
	v = bs(:,k) - As(:,:,k)*x;

	% Determine if the trajectory is clockwise or counter
	is_cw = sign(x'*[[0,-1];[1,0]]*v);
	
	% Sweep over the vector fields
	for l = 1:6

		% Don't bother with trajectory vector field
		if k != l

			% Determine if vf points in or out
			vl = bs(:,l) - As(:,:,l)*x;
			if is_cw*(vl'*[[0,-1];[1,0]]*v)/sqrt(norm(vl)*norm(v)) > 1e-5 
				is_closed = 0;
				lc = [lc,l];
			end
		end
	end

	% Calculate theta
	if !is_closed
		theta = (j-1)/grid_size*pi/3;
		break;		% Don't bother finding all closure failures
	end		
end

end
