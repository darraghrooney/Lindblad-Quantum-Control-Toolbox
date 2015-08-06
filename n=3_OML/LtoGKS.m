% This function takes a collection of Lindblad operators and returns the
% GKS matrix. Ls must be 3x3xK where K>=1

function[GKS] = LtoGKS(Ls)

% Initialize
K = size(Ls,3);
GKS = zeros(8);

% Sweep through the Lindblad operators
for k = 1:K

	% Load basis matrices
	gPs = gPaulis;

	% Find components of L op
	comps = zeros(8,1);
	comps(1) = trace(Ls(:,:,k)*gPs(:,:,1));
	comps(2) = trace(Ls(:,:,k)*gPs(:,:,2));
	comps(3) = trace(Ls(:,:,k)*gPs(:,:,3));
	comps(4) = trace(Ls(:,:,k)*gPs(:,:,4));
	comps(5) = trace(Ls(:,:,k)*gPs(:,:,5));
	comps(6) = trace(Ls(:,:,k)*gPs(:,:,6));
	comps(7) = trace(Ls(:,:,k)*gPs(:,:,7));
	comps(8) = trace(Ls(:,:,k)*gPs(:,:,8));
	
	% Add contribution from this L op to the GKS matrix
	GKS += comps*comps';
end

end
