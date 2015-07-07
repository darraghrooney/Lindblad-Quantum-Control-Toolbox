% This function computes the flag derivative of the Lindblad rates for n=3, given
% a collection of Lindblad operators and a flag.

function[deriv] = DF_w(Ls, flag)

% Get the basis for the tangent space

gens = su3off_basis();

% Initialize derivative

deriv = zeros(6,6);
count = 0;

% Sweep over Lindblad rates
for j=1:3
for k=1:3
if j != k
	

	count += 1;
	% Sweep over basis elements	
	for l = 1:6
		% Sweep over Lindblad ops			
		for m = 1:size(Ls,3)
			deriv(count,l) += flag(:,j)'*Ls(:,:,m)*gens(:,:,l)*flag(:,k)...
						*flag(:,k)'*Ls(:,:,m)'*flag(:,j);
			deriv(count,l) -= flag(:,j)'*gens(:,:,l)*Ls(:,:,m)*flag(:,k)...
						*flag(:,k)'*Ls(:,:,m)'*flag(:,j);
			deriv(count,l) += flag(:,j)'*Ls(:,:,m)*flag(:,k)...
						*flag(:,k)'*Ls(:,:,m)'*gens(:,:,l)*flag(:,j);
			deriv(count,l) -= flag(:,j)'*Ls(:,:,m)*flag(:,k)...
						*flag(:,k)'*gens(:,:,l)*Ls(:,:,m)'*flag(:,j);
		end
	end

end
end
end

% The derivative should be real, so get rid of imag machine error
deriv = real(deriv);

end
