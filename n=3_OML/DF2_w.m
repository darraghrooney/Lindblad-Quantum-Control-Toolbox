% This function computes the second flag derivative of the Lindblad rates for n=3, given
% a collection of Lindblad operators and a flag.

function[dd] = DF2_w(Ls, flag)

% Fetch the basis matrices for the tangent space

gens = su3off_basis();
for j=1:6
  gens(:,:,j) = flag*gens(:,:,j)*flag';
end

% Initialize

dd = zeros(6,6,6);
c = 0;

% Sweep over wij's
for j=1:3
for k=1:3
if j != k
	c += 1;
	
	% Sweep over component-pairs of the basis
	for l1 = 1:6
	for l2 = 1:6

	% Sweep over Lindblad operators
	for m = 1:size(Ls,3)
		
		% Compute double derivative
		dd(l1,l2,c) += flag(:,j)'*com(com(Ls(:,:,m),gens(:,:,l1)),gens(:,:,l2))*flag(:,k)*flag(:,k)'*Ls(:,:,m)'*flag(:,j);
		dd(l1,l2,c) += flag(:,j)'*Ls(:,:,m)*flag(:,k)*flag(:,k)'*com(com(Ls(:,:,m)',gens(:,:,l1)),gens(:,:,l2))*flag(:,j);
		dd(l1,l2,c) += flag(:,j)'*com(Ls(:,:,m),gens(:,:,l2))*flag(:,k)*flag(:,k)'*com(Ls(:,:,m)',gens(:,:,l1))*flag(:,j);
		dd(l1,l2,c) += flag(:,j)'*com(Ls(:,:,m),gens(:,:,l1))*flag(:,k)*flag(:,k)'*com(Ls(:,:,m)',gens(:,:,l2))*flag(:,j);
	end
	end
	end
end
end
end

% Get rid of imaginary error
dd = real(dd);

end

% Helper commutator function
function[com] = com(A,B)

	com = A*B- B*A;
end
