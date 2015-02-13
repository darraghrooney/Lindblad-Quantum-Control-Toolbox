% For a given collection of Lindblad operators, this function finds a principal flag
% and the corresponding principal rates. The principal flag is an eigenbasis of 
% the nucleus, sum([L_k,L_k']).

function [w,flag] = principal_w(Ls)

% Calculate the nucleus
nucleus = zeros(3);
for j = 1:size(Ls,3)
	nucleus += Ls(:,:,j)*Ls(:,:,j)' - Ls(:,:,j)'*Ls(:,:,j);
end

% Fetch eigeninformation
[vecs,vals] = eig(nucleus);

% Calculate the rates
w = zeros(3);
for i = 1:3
	for j = 1:3
		if (i != j)
			for k = 1:size(Ls,3);
				w(i,j) += vecs(:,i)'*Ls(:,:,k)...
				*vecs(:,j)*vecs(:,j)'*Ls(:,:,k)'*vecs(:,i); 
			end
		end
	end
end

% Discard machine error in the imaginary part
w = real(w);

% Permute the basis 
[w,Js,ind] = order_rates(w);
flag = vecs(:,ind);

end
