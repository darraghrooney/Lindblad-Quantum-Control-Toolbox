% This function takes a collection of Lindblad operators and flag
% calculates the Lindblad rates
 
function[w] = Ltow(Ls,flag)

% Find operator number
K = size(Ls,3);

% Initialize
w = zeros(4);

% Sweep over Lindblad rtes
for i=1:4
	for j=1:4

		% Skip self=rates
		if (i==j) continue; 
		end

		% Compute rate
		for k = 1:K
			sqw = flag(:,i)'*Ls(:,:,k)*flag(:,j);
			w(i,j) += sqw'*sqw;
		end
	end
end

end
