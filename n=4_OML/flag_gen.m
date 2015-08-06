% This function generates a basis of C^4, either using input angles, or randomly.
% If rand_switch is 1, the angles input is disregarded and the 12 angles are selected
% uniformly over their range. If not, angles must be a vector of 12 real numbers. 
% The first six should be on [0,pi] and the second six on [0,2pi), although this is
% not necessary. 

function[flag] = flag_gen(angles, rand_switch)

flag = eye(4);

% Generate random angles

if (rand_switch)
	phs = rand(6,1)*pi;
	ths = rand(6,1)*2*pi;

% Extract angles from input

else
	phs = angles(1:6);
	ths = angles(7:12);
end

% Calculate first basis vector and three auxiliary vectors

flag(:,1) = [cos(phs(1)/2); sin(phs(1)/2)*cos(phs(2)/2)*exp(ths(1)*I);sin(phs(1)/2)*sin(phs(2)/2)*cos(phs(3)/2)*exp(ths(1)*I+ths(2)*I);sin(phs(1)/2)*sin(phs(2)/2)*sin(phs(3)/2)*exp(ths(1)*I+ths(2)*I+ths(3)*I)];

aux1 = [-sin(phs(1)/2); cos(phs(1)/2)*cos(phs(2)/2)*exp(ths(1)*I);cos(phs(1)/2)*sin(phs(2)/2)*cos(phs(3)/2)*exp(ths(1)*I+ths(2)*I);cos(phs(1)/2)*sin(phs(2)/2)*sin(phs(3)/2)*exp(ths(1)*I+ths(2)*I+ths(3)*I)];
aux2 = [0; -sin(phs(2)/2);cos(phs(2)/2)*cos(phs(3)/2)*exp(ths(2)*I);cos(phs(2)/2)*sin(phs(3)/2)*exp(ths(2)*I+ths(3)*I)];
aux3 = [0; 0; -sin(phs(3)/2);cos(phs(3)/2)*exp(ths(3)*I)];

flag(:,2)=aux1;
flag(:,3)=aux2;
flag(:,4)=aux3;

% Generate a 3d flag to transform the auxiliaries
u3 = flag_gen_help([phs(4:6),ths(4:6)],0);

% Transform the auxiliaries
flag(:,2:4) = flag(:,2:4)*u3;

end

%% Helper function generates a flag on C^3

function[flag] = flag_gen_help(angles, rand_switch)

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
