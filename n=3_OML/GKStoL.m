% This function takes a GKS matrix and returns the collection of 
% Lindblad operators. GKS must be 8x8 and positive semi-definite

function[Ls,rnk] = GKStoL(GKS)

% Initialize
Ls = zeros(3,3,8);
rnk = 0;

% Retrieve eigeninformation
[evecs,evals] = eig(GKS);

% Load generalized Pauli matrices
gPs = gPaulis;

% Sweep through the eigenvalues
for k = 1:8

	% Don't bother with zero eigenvalues
	if (1)%abs(evals(k,k)) > 1E-10)
		rnk += 1; 	
			
		% Sweep through basis matrices
		for l = 1:8
			Ls(:,:,rnk) += sqrt(evals(k,k))*evecs(l,k)*gPs(:,:,l);
		end
		Ls(:,:,rnk) *= exp(-1i*arg(Ls(1,1,rnk)));
	end
end

% Truncate trivial L ops
Ls = Ls(:,:,1:rnk);

end
