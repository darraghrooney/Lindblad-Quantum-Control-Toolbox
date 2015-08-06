% This function takes a GKS matrix and produces the Lindblad rates w(i,j)

function[w] = GKStow(GKS,flag)


% Initialize
w = zeros(3);

% Load the basis matrices
gPs = gPaulis;

% Sweep over Lindblad rates
for i=1:3
	for j=1:3

		% Skip self-rates
		if (i==j) continue; 
		end

		% Sweep over pairs of basis matrices
		for k = 1:8
			for l = 1:8
		
				% Add contribution to rate
				w(i,j) += GKS(k,l)*flag(:,i)'*gPs(:,:,k)*flag(:,j)*...
					flag(:,j)'*gPs(:,:,l)*flag(:,i);		
			end
		end
	end
end

end
