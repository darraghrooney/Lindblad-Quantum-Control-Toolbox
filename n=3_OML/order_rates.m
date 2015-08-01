% This function takes rates w_ij for n=3 and permutes the indices so that the 
% input rates J_k are ordered J_1 >= J_2 >= J_3. The input rate is defined to 
% be J_k = sum w(i1,j1)w(i2,j2) where the edges (i,j) form a tree on
% the three vertices, with root k. The sum is done over all such trees.

function [w, Jvec, ind] = order_rates(w);

% Calculate input rates

J1 = 	w(1,2)*w(1,3) + w(1,2)*w(2,3) + w(1,3)*w(3,2);
J2 = 	w(2,1)*w(2,3) + w(2,1)*w(1,3) + w(2,3)*w(3,1);
J3 = 	w(3,1)*w(3,2) + w(3,1)*w(1,2) + w(3,2)*w(2,1);

Jvec = [J1,J2,J3];

% Check to see that the system is non-degenerate

if J1 + J2 + J3 == 0
    	error('degenerate A')

% Check to see that the input rates are ordered

elseif J3 > J2 || J2 > J1

	% Re-order

	holding = w;

	[S,ind] = sort([J1,J2,J3]);

	w(1,2) = holding(ind(3),ind(2));
	w(2,1) = holding(ind(2),ind(3));
	w(1,3) = holding(ind(3),ind(1));
	w(3,1) = holding(ind(1),ind(3));
	w(2,3) = holding(ind(2),ind(1));
	w(3,2) = holding(ind(1),ind(2));

	J3 = S(1);
	J2 = S(2);	
	J1 = S(3);
	Jvec = [J1,J2,J3];
	ind = ind(3:-1:1); 
end

