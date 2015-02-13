% This function generates a representation of the symmetric group S4
% onto matrices in R3. Specifically, the representation matrices 
% permute the following vertices of a regular tetrahedron:
%  	v1: (1,1/sqrt(3),1/sqrt(6))	v2: (-1,1/sqrt(3),1/sqrt(6))
%	v3: (0,-2/sqrt(3),1/sqrt(6))	v3: (0,0,-3/sqrt(6))
% So for example g12 exchanges v1 and v2 and g123 cycles v1, v2 and v3

function[gens] = S4gen()

% Initialize the 24 matrices. The first is the group identity
gens = zeros(3,3,24);
gens(:,:,1) = eye(3,3);

% Representations of 2-cycles
g12 = diag([-1,1,1],0);
g23 = diag([1/2,-1/2,1],0);
g23(1,2) = sqrt(3)/2;
g23(2,1) = sqrt(3)/2;
g34 = diag([1,1/3,-1/3],0);
g34(2,3) = 4/sqrt(18);
g34(3,2) = 4/sqrt(18);
g13 = g12*g23*g12;
g14 = g13*g34*g13;
g24 = g23*g34*g23; 

gens(:,:,2) = g12;
gens(:,:,3) = g13;
gens(:,:,4) = g14;
gens(:,:,5) = g23;
gens(:,:,6) = g24;
gens(:,:,7) = g34;

% Representations of two disjoint 2-cycles
g12g34 = g12*g34; 
g14g23 = g14*g23; 
g13g24 = g13*g24;

gens(:,:,8) = g12g34;
gens(:,:,9) = g13g24;
gens(:,:,10) = g14g23;

% Representations of 3-cycles
g123 = g12*g23;
g132 = g13*g23;
g124 = g12*g24;
g142 = g14*g24;
g134 = g13*g34;
g143 = g14*g34;
g234 = g23*g34;
g243 = g24*g34;

gens(:,:,11) = g123;
gens(:,:,12) = g132;
gens(:,:,13) = g124;
gens(:,:,14) = g142;
gens(:,:,15) = g134;
gens(:,:,16) = g143;
gens(:,:,17) = g234;
gens(:,:,18) = g243;

% Representations of 4-cycles
g1234 = g12*g23*g34;
g1243 = g12*g24*g34;
g1324 = g13*g23*g24;
g1342 = g13*g34*g24;
g1423 = g14*g24*g23;
g1432 = g14*g34*g23;

gens(:,:,19) = g1234;
gens(:,:,20) = g1243;
gens(:,:,21) = g1324;
gens(:,:,22) = g1342;
gens(:,:,23) = g1423;
gens(:,:,24) = g1432;

