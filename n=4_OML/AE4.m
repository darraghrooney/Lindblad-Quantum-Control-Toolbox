% This function plots the likely boundaries for the STLC set for
% an n=4 Lindblad system. "Likely" meaning only the 92 surfaces
% corresponding to neighboring sectors. The total number of surfaces
% is 2024, but this number of meshes can overwhelm the CPU, and the 
% other 1932 surfaces are usually interior. 

function [] = AE4(w,grid_size)

% Check that the co-efficients are non-negative

if ( max(max(w < -1E-15)) == 1)
	error('Co-efficients must be non-negative');
end

% Calculate input rates and re-order if necessary

[w,Jv] = order_rates(w);

% Calculate A(w_ij)
 
[A,b,rinf] = A_b_rinf(w);

% Calculate the images of A(w_ij) under action of the symmetric group

AE = zeros(3,3,24);
bE = zeros(3,24);
AE(:,:,1) = A;
bE(:,1) = b;

S4g = S4gen();

for j=2:24
	AE(:,:,j) = S4g(:,:,j)*A*S4g(:,:,j)'; 
	bE(:,j) = S4g(:,:,j)*b;
end

% Initialize mesh

tomesh=zeros(grid_size+1,grid_size+1,3,92);

% Calculate the 12 surfaces corresponding to the (1 2)(3 4) subgroup and rotations

tomesh(:,:,:,1) = AE_helper(AE(:,:,[1 2 7 8]),bE(:,[1 2 7 8]),grid_size);
tomesh(:,:,:,2) = AE_helper(AE(:,:,[2 1 8 7]),bE(:,[2 1 8 7]),grid_size);

tomesh(:,:,:,3) = rot_surf(tomesh(:,:,:,1),11);
tomesh(:,:,:,4) = rot_surf(tomesh(:,:,:,2),11);
tomesh(:,:,:,5) = rot_surf(tomesh(:,:,:,1),12);
tomesh(:,:,:,6) = rot_surf(tomesh(:,:,:,2),12);
tomesh(:,:,:,7) = rot_surf(tomesh(:,:,:,1),13);
tomesh(:,:,:,8) = rot_surf(tomesh(:,:,:,2),13);
tomesh(:,:,:,9) = rot_surf(tomesh(:,:,:,1),14);
tomesh(:,:,:,10) = rot_surf(tomesh(:,:,:,2),14);
tomesh(:,:,:,11) = rot_surf(tomesh(:,:,:,3),17);
tomesh(:,:,:,12) = rot_surf(tomesh(:,:,:,4),17);


% Calculate the 40 surfaces corresponding to the (1 2 3) subgroup and rotations

tomesh(:,:,:,13) = AE_helper(AE(:,:,[1 2 5 3]),bE(:,[1 2 5 3]),grid_size);
tomesh(:,:,:,14) = AE_helper(AE(:,:,[1 11 12 3]),bE(:,[1 11 12 3]),grid_size);
tomesh(:,:,:,15) = AE_helper(AE(:,:,[2 1 11 3]),bE(:,[2 1 11 3]),grid_size);
tomesh(:,:,:,16) = AE_helper(AE(:,:,[2 1 12 5]),bE(:,[2 1 12 5]),grid_size);
tomesh(:,:,:,17) = AE_helper(AE(:,:,[11 2 3 12]),bE(:,[11 2 3 12]),grid_size);
tomesh(:,:,:,18) = AE_helper(AE(:,:,[2 5 11 12]),bE(:,[2 5 11 12]),grid_size);
tomesh(:,:,:,19) = AE_helper(AE(:,:,[2 5 11 12]),bE(:,[2 5 11 12]),grid_size);
tomesh(:,:,:,20) = AE_helper(AE(:,:,[2 1 3 12]),bE(:,[2 1 3 12]),grid_size);
tomesh(:,:,:,21) = AE_helper(AE(:,:,[5 2 12 11]),bE(:,[5 2 12 11]),grid_size);
tomesh(:,:,:,22) = AE_helper(AE(:,:,[1 5 11 3]),bE(:,[1 5 11 3]),grid_size);

for k = 1:10
	tomesh(:,:,:,22+k) = rot_surf(tomesh(:,:,:,12+k),17);
	tomesh(:,:,:,32+k) = rot_surf(tomesh(:,:,:,12+k),18);
	tomesh(:,:,:,42+k) = rot_surf(tomesh(:,:,:,12+k),13);
end

% Calculate the 40 surfaces corresponding to the (2 3 4) subgroup and rotations

tomesh(:,:,:,53) = AE_helper(AE(:,:,[1 5 7 6]),bE(:,[1 5 7 6]),grid_size);
tomesh(:,:,:,54) = AE_helper(AE(:,:,[1 17 18 6]),bE(:,[1 17 18 6]),grid_size);
tomesh(:,:,:,55) = AE_helper(AE(:,:,[5 1 18 6]),bE(:,[5 1 18 6]),grid_size);
tomesh(:,:,:,56) = AE_helper(AE(:,:,[5 1 17 7]),bE(:,[5 1 17 7]),grid_size);
tomesh(:,:,:,57) = AE_helper(AE(:,:,[18 5 6 17]),bE(:,[18 5 6 17]),grid_size);
tomesh(:,:,:,58) = AE_helper(AE(:,:,[17 7 6 1]),bE(:,[17 7 6 1]),grid_size);
tomesh(:,:,:,59) = AE_helper(AE(:,:,[5 18 7 17]),bE(:,[5 18 7 17]),grid_size);
tomesh(:,:,:,60) = AE_helper(AE(:,:,[17 1 6 5]),bE(:,[17 1 6 5]),grid_size);
tomesh(:,:,:,61) = AE_helper(AE(:,:,[18 5 17 7]),bE(:,[18 5 17 7]),grid_size);
tomesh(:,:,:,62) = AE_helper(AE(:,:,[1 5 18 6]),bE(:,[1 5 18 6]),grid_size);

for k = 1:10
	tomesh(:,:,:,62+k) = rot_surf(tomesh(:,:,:,52+k),11);
	tomesh(:,:,:,72+k) = rot_surf(tomesh(:,:,:,52+k),12);
	tomesh(:,:,:,82+k) = rot_surf(tomesh(:,:,:,52+k),14);
end

% Plot the meshes

hold on

for k = 1:92
	mesh(tomesh(:,:,1,k),tomesh(:,:,2,k),tomesh(:,:,3,k));	
end

% Adapt plot window
extent = norm(A\b);
if (extent != 0)
	axis([-extent,extent,-extent,extent,-1.1*extent,extent]);
end
view(3)
hold off
end

% This is a helper function that calculates surfaces
% Requires four A's and four b's as well as grid_size specification.
% Each surface requires only three, but since each surface is "triangular"
% and meshes are "square", the surfaces are produced in pairs that share
% an edge (= "diagonal" of the mesh). The 2nd and 3rd sections in As and bs
% correspond to the diagonal/shared points.

function[tomesh] = AE_helper(As,bs,grid_size)

% Initialize
tomesh = zeros(grid_size+1,grid_size+1,3);
[s,t]=meshgrid((0:grid_size)/grid_size);

% Sweep over grid
for j=1:(grid_size+1)
	for k=1:(grid_size+2-j)

		% Calculate surface point
		invA = s(j,k)*As(:,:,2) + t(j,k)*As(:,:,3) + ...
			(1-s(j,k)-t(j,k))*As(:,:,1);
	        bb = s(j,k)*bs(:,2) + t(j,k)*bs(:,3) + (1-s(j,k)-t(j,k))*bs(:,1);
		tomesh(j,k,:) = invA\bb;

		% Calculate surface point for paired surface
		if (j != grid_size + 2 - k)
			invA = s(j,k)*As(:,:,2) + t(j,k)*As(:,:,3) + ...
				(1-s(j,k)-t(j,k))*As(:,:,4);
		        bb = s(j,k)*bs(:,2) + t(j,k)*bs(:,3) + ...
				(1-s(j,k)-t(j,k))*bs(:,4);
			tomesh(grid_size+2-k,grid_size+2-j,:) = invA\bb;
		end	
	end
end

end
