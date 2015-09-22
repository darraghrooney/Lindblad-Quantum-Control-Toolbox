% This function creates a random 3x3 matrix that is Hermitian and positive
% semi-definite but with purely imaginary off-diagonal terms. The vector of 
% diagonal terms, a, is decreasing and the first is normalized to 100. The second two
% are uniform on the simplex with vertices (0,100), (100,0) and (0,0).

% The vector of off-diagonal terms, b, is uniform on the ball with radius 
% 4(sqrt(a1 a2)+sqrt(a2 a3)+sqrt(a1 a3)).
 
% The function returns a, b, and A, the matrix.

function [ b, A ] = randb(a)
    b=zeros(3,1);

% Construct b
    bstar = randn(3,1);
    bstar = bstar/norm(bstar);
    r = rand;
    bstar = r^(1/3)*bstar;
    b(1)=2*sqrt(a(2)*a(3))*bstar(1);
    b(2)=2*sqrt(a(3)*a(1))*bstar(2);
    b(3)=2*sqrt(a(1)*a(2))*bstar(3);
    
    
% Construct A
    A=[[a(1),1i*b(3)/2,-1i*b(2)/2];[-1i*b(3)/2,a(2),1i*b(1)/2];[1i*b(2)/2,-1i*b(1)/2,a(3)]];

end

