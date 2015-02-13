% This function generates a representation of the symmetric group S3
% onto matrices in R2. Specifically, the representation matrices 
% permute the following vertices of a regular tetrahedron:
%  	v1: (1,1/sqrt(3))	v2: (-1,1/sqrt(3))
%	v3: (0,-2/sqrt(3))
% So for example g12 exchanges v1 and v2 and g123 cycles v1, v2 and v3

function[gens] = S3gen()

% Initialize the 6 matrices. The first is the group identity
gens = zeros(2,2,6);
gens(:,:,1) = eye(2);

% Representations of 2-cycles
g12 = diag([-1,1]);
g23 = diag([1/2,-1/2]);
g23(1,2) = sqrt(3)/2;
g23(2,1) = sqrt(3)/2;
g13 = g12*g23*g12;

gens(:,:,2) = g12;
gens(:,:,3) = g13;
gens(:,:,4) = g23;

% Representations of 3-cycles
g123 = g12*g23;
g132 = g13*g23;

gens(:,:,5) = g123;
gens(:,:,6) = g132;

end
