% This function calculates A(w) and b(w) for an n=4 Lindblad system.
% w must be input as a 4x4 matrix (the diagonal entries are ignored).
% The asymptotic end-state rinf = A^(-1)b is also returned if unique.
% If not unique (i.e. A is not invertible), the pure state is returned.

function[A,b,rinf] = A_b_rinf(w)

% Calculate A(w)

A=zeros(3,3);

A(1,1)=1/2*(2*w(1,2)+2*w(2,1)+w(3,1)+w(3,2)+w(4,1)+w(4,2)); 
A(2,2)=1/6*(6*w(1,3)+6*w(2,3)+3*w(3,1)+3*w(3,2)+4*w(4,3)+w(4,1)+w(4,2));
A(3,3)=1/3*(w(4,1)+w(4,2)+w(4,3)+3*(w(1,4)+w(2,4)+w(3,4)));
A(1,2)=1/2/sqrt(3)*(2*w(2,1)-2*w(1,2)+2*w(1,3)-2*w(2,3)+w(3,1)-w(3,2)+ ...
	w(4,1)-w(4,2)); 
A(2,1)=1/sqrt(3)/2*(3*w(3,1)-3*w(3,2)+w(4,1)-w(4,2));
A(1,3)=1/2/sqrt(6)*(3*(w(1,4)-w(2,4))+2*(w(2,1)-w(1,2))+w(2,3)-w(1,3)+ ...
	w(3,1)-w(3,2)+w(4,1)-w(4,2));
A(3,1)=2/sqrt(6)*(w(4,1)-w(4,2));
A(2,3)=1/6/sqrt(2)*(3*(w(1,4)-w(1,3)+w(2,4)-w(2,3)+w(3,1)+w(3,2))- ...
	6*w(3,4)+w(4,1)+w(4,2)-2*w(4,3));
A(3,2)=sqrt(2)/3*(w(4,1)+w(4,2)-2*w(4,3));

% Calculate b(w)

b=zeros(3,1);

b(1)=1/4*(2*w(1,2)-2*w(2,1)+w(1,3)-w(3,1)+w(3,2)-w(2,3)+w(1,4)- ...
	w(4,1)+w(4,2)-w(2,4)); 
b(2)=1/4/sqrt(3)*(3*w(1,3)-3*w(3,1)+3*w(2,3)-3*w(3,2)+2*w(4,3)- ...
	2*w(3,4)+w(1,4)-w(4,1)+w(2,4)-w(4,2));
b(3)=1/sqrt(6)*(w(1,4)-w(4,1)+w(2,4)-w(4,2)+w(3,4)-w(4,3));

% Calculate rinf = A^-1 *b

if (det(A) > 0)
	rinf = A\b;
else
	rinf = [1,1/sqrt(3),1/sqrt(6)];
end

end
