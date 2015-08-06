% This is an increment function for sweeping over the set of flags in C^3
% This set is a 6-dim manifold: three angles live on [0,pi], three on [0,2*pi):
% 	psi1 = cos(phi1/2)*e1 + sin(phi1/2)e^(i*theta1)*cos(phi2/2)*e2
%	 	+ sin(phi1/2)e^(i*theta1)*sin(phi2/2)*e^(i*theta2)*e3
% 	psi2' = -sin(phi1/2)*e1 + cos(phi1/2)e^(i*theta1)*cos(phi2/2)*e2
%	 	+ sin(phi1/2)e^(i*theta1)*sin(phi2/2)*e^(i*theta2)*e3
%	psi3' = cos(phi2/2)*e2 + sin(phi2/2)*e^(i*theta2)*e3
%	psi2 = cos(phi3/2)*psi2' + sin(phi3/2)*e^(i*theta3)*psi3'
%	psi3 = -sin(phi3/2)*psi2' + cos(phi3/2)*e^(i*theta3)*psi3'

% The flag is (psi1,psi2,psi3) and psi2' and psi3' are auxiliary vectors.

% All six parameters are discretized by the same factor (grid_size). The function
% takes the old flag as a 3x3 matrix, the auxiliary vectors as a 3x2 matrix, and also
% requires the six indices. The 1st, 3rd and 5th indices run 0 ... grid_size since the
% phi's live on a closed interval, and the 2nd, 4th and 6th run 0 ... grid_size-1. 
% Additionally, when a psi is 0 or pi, the accompanying theta only runs one value ...
% others would be redundant.

% So the total number of points are: ( n^2- n + 2)^3.


function[new_flag,new_aux,new_indices] = ...
	flag_inc(old_flag,old_aux,old_indices,grid_size)

% Initialization

new_flag = old_flag;
new_aux = old_aux;
new_indices = old_indices;

% Some useful parameters

cs = cos(pi/grid_size/2);
sn = sin(pi/grid_size/2);
del = 2*pi/grid_size;
	
% If j3, k3 = 0, go straight to j3 = 1 (j3, k3 being the indices for phi3, theta3).

if ( new_indices(6) == 0 && new_indices(5) == 0)
	new_flag(:,2) = cs*new_aux(:,1) + sn*new_aux(:,2);
	new_flag(:,3) = -sn*new_aux(:,1) + cs*new_aux(:,2);
	new_indices(5) += 1;

% If j3 isn't 0 or grid_size, and k3 isn't maxed out, increment k3.

elseif (new_indices(5) != grid_size && new_indices(6) < grid_size-1) 
	th = del*new_indices(6);
	fact = (exp(1i*(th+del))-exp(1i*th));

	new_flag(:,2) += sin(new_indices(5)*pi/grid_size/2)*new_aux(:,2)*fact;
	new_flag(:,3) += cos(new_indices(5)*pi/grid_size/2)*new_aux(:,2)*fact;
	new_indices(6) += 1;

% If j3 isn't 0 or grid_size, and k3 is maxed out, increment j3 and zero k3.

elseif (new_indices(5) != grid_size && new_indices(6) == grid_size-1) 
	ph = pi*(new_indices(5)+1)/grid_size;
	ncs = cos(ph/2);
	nsn = sin(ph/2);
	new_flag(:,2) = ncs*new_aux(:,1)+nsn*new_aux(:,2);
	new_flag(:,3) = -nsn*new_aux(:,1)+ncs*new_aux(:,2);
	new_indices(6) = 0;
	new_indices(5) += 1;

% if j3 = grid_size, zero j3 and look at j2, k2

else
	new_indices(5:6) *= 0;

% If j2, k2 = 0, go straight to j2 = 1.

if ( new_indices(4) == 0 && new_indices(3) == 0)

	new_flag(3,1) = sn*new_flag(2,1);
	new_flag(2,1) *= cs;
	new_aux(3,1) = sn*new_aux(2,1);
	new_aux(2,1) *= cs;
	new_aux(:,2) = [0;-sn;cs];

	new_flag(:,2) = new_aux(:,1);
	new_flag(:,3) = new_aux(:,2);
	new_indices(3) += 1;

% If j2 isn't 0 or grid_size, and k2 isn't maxed out, increment k2.

elseif (new_indices(3) != grid_size && new_indices(4) < grid_size-1) 
	th = del*new_indices(6);
	fact = (exp(1i*(th+del))-exp(1i*th));

	new_flag(3,1) *= exp(1i*del);
	new_aux(3,1) *= exp(1i*del);
	new_aux(3,2) *= exp(1i*del);

	new_flag(:,2) = new_aux(:,1);
	new_flag(:,3) = new_aux(:,2);
	new_indices(4) += 1;

% If j2 isn't 0 or grid_size, and k2 is maxed out, increment j2 and zero k2.
	
elseif (new_indices(3) != grid_size && new_indices(4) == grid_size-1) 
	newph = (new_indices(3) + 1)*pi/grid_size;
	newcs = cos(newph/2);
	newsn = sin(newph/2);
	topph = new_indices(1)*pi/grid_size;
	topcs = cos(topph/2);
	topsn = sin(topph/2);
	topexp = exp(1i*del*new_indices(2));

	new_flag(2,1) = topsn*topexp*newcs;
	new_flag(3,1) = topsn*topexp*newsn;
	new_aux(2,1) = topcs*topexp*newcs;
	new_aux(3,1) = topcs*topexp*newsn;
	new_aux(2,2) = -newsn;
	new_aux(3,2) = newcs;
	
	new_flag(:,2) = new_aux(:,1);
	new_flag(:,3) = new_aux(:,2);
	new_indices(4) = 0;
	new_indices(3) += 1;

% If j2 = grid_size, zero j2 and look at j1,k1.

else
	new_indices(3:4) *= 0;

	new_aux(:,2) = [0;0;1];
	new_flag(:,3) = new_aux(:,2);

% If j1, k1 = 0, increment j1.

if (new_indices(2) == 0 && new_indices(1) == 0)
	new_flag(:,1) = [cs;sn;0];
	new_aux(:,1) = [-sn;cs;0];
	new_flag(:,2) = new_aux(:,1);
	new_indices(1) += 1;

% If j1 isn't 0 or grid_size, and k1 isn't maxed out, increment k1.

elseif (new_indices(1) != grid_size && new_indices(2) < grid_size-1) 
	
	topph = new_indices(1)*pi/grid_size;
	topcs = cos(topph/2);
	topsn = sin(topph/2);
	topexp = exp(1i*del*new_indices(2));
	
	new_flag(2,1) = topsn*topexp*exp(1i*del);
	new_flag(3,1) = 0;
	new_aux(2,1) = topcs*topexp*exp(1i*del);
	new_aux(3,1) = 0;
	new_flag(:,2) = new_aux(:,1);
	new_indices(2) += 1;

% If j1 isn't 0 or grid_size, and k1 is maxed out, increment j1.

elseif (new_indices(1) != grid_size && new_indices(2) == grid_size-1) 

	nnph = (new_indices(1) + 1)/grid_size*pi;
	nncs = cos(nnph/2);
	nnsn = sin(nnph/2);
	new_flag(:,1) = [nncs;nnsn;0];
	new_aux(:,1) = [-nnsn;nncs;0];

	new_flag(:,2) = new_aux(:,1);
	new_indices(1) += 1;
	new_indices(2) = 0;

% If j1 = grid_size, reset the flag to starting position.

else
	new_indices = zeros(6,1);
	new_aux = zeros(3,2); new_aux(2,1) = 1; new_aux(3,2) = 1;
	new_flag = eye(3);

end
end
end

end

