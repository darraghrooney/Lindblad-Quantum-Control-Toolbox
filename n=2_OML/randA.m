% This function creates a random 3x3 matrix that is Hermitian and positive
% semi-definite but with purely imaginary off-diagonal terms. The vector of 
% diagonal terms, a, is decreasing and the first is normalized to 100. The second two
% are uniform on the simplex with vertices (0,100), (100,0) and (0,0).

% The vector of off-diagonal terms, b, is uniform on the ball with radius 
% 4(sqrt(a1 a2)+sqrt(a2 a3)+sqrt(a1 a3)).
 
% The function returns a, b, and A, the matrix.

% EDIT: if desired, a random rotation is generated and A is rotated
s
function [ a, b, A ] = randA(rot_switch)
    a=zeros(3,1);
    b=zeros(3,1);

% Construct a
    a(1)=100;
    temp=100*[rand() rand()]; 
    a(2)=max(temp);
    a(3)=min(temp);
    
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

    
% If a random rotation is desired:
if (rot_switch)

  % Generate random quaternion
  us = rand(3,1);
  q = [sqrt(1-us(1))*[sin(2*pi*us(2));cos(2*pi*us(2))];sqrt(us(1))*[sin(2*pi*us(3));cos(2*pi*us(3))];];
  % Map quaternion to SO(3)
  rotation = eye(3);
  rotation(1,1) -= 2*q(3)^2 + 2*q(4)^2;
  rotation(2,2) -= 2*q(2)^2 + 2*q(4)^2;
  rotation(3,3) -= 2*q(2)^2 + 2*q(3)^2;
  rotation(1,2) += 2*q(2)*q(3) - 2*q(1)*q(4);
  rotation(2,1) += 2*q(2)*q(3) + 2*q(1)*q(4);
  rotation(3,1) += 2*q(2)*q(4) - 2*q(1)*q(3);
  rotation(1,3) += 2*q(2)*q(4) + 2*q(1)*q(3);
  rotation(2,3) += 2*q(3)*q(4) - 2*q(1)*q(2);
  rotation(3,2) += 2*q(3)*q(4) + 2*q(1)*q(2);
  
  % Rotate
  A = rotation*A*rotation';
  
endif
 
end

