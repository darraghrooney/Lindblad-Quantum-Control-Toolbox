% This function calculates A(w) and b(w) for an n=3 Lindblad system.
% w must be input as a 3x3 matrix (the diagonal entries are ignored).
% The asymptotic end-state rinf = A^(-1)b is also returned if unique.
% If not unique (i.e. A is not invertible), the pure state is returned.

function[A,b,rinf] = A_b_rinf(w)

% Calculate A(w)

A=zeros(2,2);

A(1,1)=1/2*(2*w(1,2)+2*w(2,1)+w(3,1)+w(3,2)); 
A(2,2)=1/2*(2*w(1,3)+2*w(2,3)+w(3,1)+w(3,2));
A(1,2)=1/2/sqrt(3)*(2*w(2,1)-2*w(1,2)+2*w(1,3)-2*w(2,3)+w(3,1)-w(3,2)); 
A(2,1)=sqrt(3)/2*(w(3,1)-w(3,2));

% Calculate b(w)

b=zeros(2,1);

b(1)=1/3*(2*w(1,2)-2*w(2,1)+w(1,3)-w(3,1)+w(3,2)-w(2,3)); 
b(2)=1/sqrt(3)*(w(1,3)-w(3,1)+w(2,3)-w(3,2));

% Calculate rinf = A^-1 *b

if (det(A) > 0)
	rinf = A\b;
else
	rinf = [1,1/sqrt(3)];
end

end
