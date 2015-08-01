% This function calculates A(w) and b(w) for an n=4 Lindblad system.
% w must be input as a 4x4 matrix (the diagonal entries are ignored).
% The asymptotic end-state rinf = A^(-1)b is also returned if unique.
% If not unique (i.e. A is not invertible), the pure state is returned.

function[A,b,rinf] = A_b_rinf(w)

% Calculate A(w)

A=zeros(3,3);
Omega = w - diag(sum(w));
Pi = [[1,-1,0,0]/sqrt(2);[1,1,-2,0]/sqrt(6);[1,1,1,-3]/sqrt(12)];
A = -Pi*Omega*Pi';

% Calculate b(w)

b=zeros(3,1);
iota = [1;1;1;1]/4;
b = Pi*Omega*iota;

% Calculate rinf = A^-1 *b

if (det(A) > 0)
	rinf = A\b;
else
	rinf = [1,1/sqrt(3),1/sqrt(6)]/sqrt(2);
end

end
