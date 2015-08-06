% This function generates a basis of C^3, either using input angles, or randomly.
% If rand_switch is 1, the angles input is disregarded and the six angles are selected
% uniformly over their range. If not, angles must be a vector of six real numbers. 
% The first three should be on [0,pi] and the second three on [0,2pi), although this is
% not necessary. 

function[flag] = flag3d_gen(angles, rand_switch)

flag = eye(3);

% Generate random angles

if (rand_switch)
	phs = rand(3,1)*pi;
	ths = rand(3,1)*2*pi;

% Extract angles from input

else
	phs = angles(1:3);
	ths = angles(4:6);
end

% Calculate first basis vector and two auxiliary vectors

flag(1:3,1) = [cos(phs(1)/2); sin(phs(1)/2)*cos(phs(2)/2)*exp(ths(1)*I);sin(phs(1)/2)*sin(phs(2)/2)*exp(ths(1)*I+ths(2)*I)];
aux1 = [-sin(phs(1)/2); cos(phs(1)/2)*cos(phs(2)/2)*exp(ths(1)*I);cos(phs(1)/2)*sin(phs(2)/2)*exp(ths(1)*I+ths(2)*I)];
aux2 = [0;-sin(phs(2)/2);cos(phs(2)/2)*exp(ths(2)*I)];

% Calculate second and third vectors from the auxiliaries.

flag(1:3,2) = cos(phs(3)/2)*aux1 + sin(phs(3)/2)*exp(I*ths(3))*aux2;
flag(1:3,3) = -sin(phs(3)/2)*aux1 + cos(phs(3)/2)*exp(I*ths(3))*aux2;


end
