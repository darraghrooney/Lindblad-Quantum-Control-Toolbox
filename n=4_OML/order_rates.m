% This function takes rates w_ij for n=4 and permutes the indices so that the 
% input rates J_k are ordered J_1 >= ... >= J_4. The input rate is defined to 
% be J_k = sum w(i1,j1)w(i2,j2)w(i3,j3) where the edges (i,j) form a tree on
% the four vertices, with root k. The sum is done over all such trees.

function [w, Jvec] = order_rates(w);

% Calculate input rates

J1 = 	w(1,2)*w(1,3)*w(1,4) + w(1,2)*w(2,3)*w(3,4) + w(1,2)*w(2,4)*w(4,3) + ...
	w(1,3)*w(3,2)*w(2,4) + w(1,3)*w(3,4)*w(4,2) + w(1,4)*w(4,3)*w(3,2) + ...
	w(1,4)*w(4,2)*w(2,3) + w(1,2)*w(1,3)*w(2,4) + w(1,2)*w(1,3)*w(3,4) + ...
	w(1,2)*w(1,4)*w(2,3) + w(1,2)*w(1,4)*w(4,3) + w(1,3)*w(1,4)*w(3,2) + ...
	w(1,3)*w(1,4)*w(4,2) + w(1,2)*w(2,3)*w(2,4) + w(1,3)*w(3,2)*w(3,4) + ...
	w(1,4)*w(4,3)*w(4,2);
J2 = 	w(2,1)*w(2,3)*w(2,4) + w(2,1)*w(1,3)*w(3,4) + w(2,1)*w(1,4)*w(4,3) + ...
	w(2,3)*w(3,1)*w(1,4) + w(2,3)*w(3,4)*w(4,1) + w(2,4)*w(4,3)*w(3,1) + ...
	w(2,4)*w(4,1)*w(1,3) + w(2,1)*w(2,3)*w(1,4) + w(2,1)*w(2,3)*w(3,4) + ...
	w(2,1)*w(2,4)*w(1,3) + w(2,1)*w(2,4)*w(4,3) + w(2,3)*w(2,4)*w(3,1) + ...
	w(2,3)*w(2,4)*w(4,1) + w(2,1)*w(1,3)*w(1,4) + w(2,3)*w(3,1)*w(3,4) + ...
	w(2,4)*w(4,3)*w(4,1);
J3 = 	w(3,1)*w(3,2)*w(3,4) + w(3,1)*w(1,2)*w(2,4) + w(3,1)*w(1,4)*w(4,2) + ...
	w(3,2)*w(2,1)*w(1,4) + w(3,2)*w(2,4)*w(4,1) + w(3,4)*w(4,2)*w(2,1) + ...
	w(3,4)*w(4,1)*w(1,2) + w(3,1)*w(3,2)*w(1,4) + w(3,1)*w(3,2)*w(2,4) + ...
	w(3,1)*w(3,4)*w(1,2) + w(3,1)*w(3,4)*w(4,2) + w(3,2)*w(3,4)*w(2,1) + ...
	w(3,2)*w(3,4)*w(4,1) + w(3,1)*w(1,2)*w(1,4) + w(3,2)*w(2,1)*w(2,4) + ...
	w(3,4)*w(4,2)*w(4,1);
J4 = 	w(4,1)*w(4,2)*w(4,3) + w(4,1)*w(1,2)*w(2,3) + w(4,1)*w(1,3)*w(3,2) + ...
	w(4,2)*w(2,1)*w(1,3) + w(4,2)*w(2,3)*w(3,1) + w(4,3)*w(3,2)*w(2,1) + ...
	w(4,3)*w(3,1)*w(1,2) + w(4,1)*w(4,2)*w(1,3) + w(4,1)*w(4,2)*w(2,3) + ...
	w(4,1)*w(4,3)*w(1,2) + w(4,1)*w(4,3)*w(3,2) + w(4,2)*w(4,3)*w(2,1) + ...
	w(4,2)*w(4,3)*w(3,1) + w(4,1)*w(1,2)*w(1,3) + w(4,2)*w(2,1)*w(2,3) + ...
	w(4,3)*w(3,2)*w(3,1);
Jvec = [J1,J2,J3,J4];

% Check to see that the system is non-degenerate

if J1 + J2 + J3 + J4 == 0
    	error('degenerate A')

% Check to see that the input rates are ordered

elseif J4 > J3 || J3 > J2 || J2 > J1
	disp('Input rates not ordered:');
	disp([J1,J2,J3,J4]);
	disp('Re-ordering w_ij');

	% Re-order

	holding = w;

	[S,I] = sort([J1,J2,J3,J4]);

	w(1,2) = holding(I(4),I(3));
	w(2,1) = holding(I(3),I(4));
	w(1,3) = holding(I(4),I(2));
	w(3,1) = holding(I(2),I(4));
	w(2,3) = holding(I(3),I(2));
	w(3,2) = holding(I(2),I(3));
	w(1,4) = holding(I(4),I(1));
	w(4,1) = holding(I(1),I(4));
	w(2,4) = holding(I(3),I(1));
	w(4,2) = holding(I(1),I(3));
	w(3,4) = holding(I(2),I(1));
	w(4,3) = holding(I(1),I(2));

	J4 = S(1);	
	J3 = S(2);
	J2 = S(3);	
	J1 = S(4);
	Jvec = [J1,J2,J3,J4];
end

