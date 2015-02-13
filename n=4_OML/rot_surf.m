%% This take a meshed surface in R^3 and rotates it using a generator of 
%% the symmetric group S4. "r" specifies which generator to use.

function[rsurf] = rot_surf(srf,r)

% Initialization
S4g = S4gen();
rsurf = 0*srf;

% Sweep over points in the mesh
for j=1:size(srf(:,:,1),1)
	for k=1:size(srf(:,:,1),2)

		% Use the specified generator to rotate the given vector
		for l = 1:3
			inp = S4g(:,:,r)*vec(srf(j,k,:));
			rsurf(j,k,l)=inp(l);
		end
	end
end

