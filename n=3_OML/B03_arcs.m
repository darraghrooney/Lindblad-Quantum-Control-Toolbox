% This function takes the boundary of a SLC set for an n=3 Lindblad control system
% with discretized control set, and calculates the first-level trajectories that go
% into finding the GC set, aka B^0_3. These trajectories start at either the end-states
% or points where a vector field is parallel to an AE3 arc. Two grid_sizes are required:
% one for the number of points that discretize each AE3 arc, and one for the the number
% of points in each trajectory.

function[all_traj] = B03_arcs(w, grid_size1, grid_size2, plot_switch)

% Find the six b-vectors and A-matrices

[A1,b1,rinf] = A_b_rinf(w);
S3g = S3gen;
As = zeros(2,2,6);
bs = zeros(2,6);
rinfs = zeros(2,6);
for j = 1:6
	As(:,:,j) = S3g(:,:,j)*A1*S3g(:,:,j)';
	bs(:,j) = S3g(:,:,j)*b1;
	rinfs(:,j) = S3g(:,:,j)*rinf;
end

% Find trajectories from the end-states (there are 6*5 = 30).

es_traj = zeros(2, grid_size2+3, 30);
count = 0;

for j = 1:6
for k = 1:6
if j ~= k

	count += 1;
	es_traj(:,1:grid_size2+2,count) = trajec(rinfs(1,j), rinfs(2,j), As(:,:,k),...
			 bs(:,k), grid_size2, 0);
	es_traj(:,grid_size2+3,count) = k*[1;1];

end
end
end

% Find launch points

Jm = [[0,1];[-1,0]];
launches = [];

% Sweep over AE3 arcs

for j = 1:5
for k = j+1:6

	% Sweep over points on arc	

	for l = 1:grid_size1-1

		% Find launch points

		s1 = (l-1)/grid_size1;
		s2 = l/grid_size1;
		x1 = (s1*As(:,:,j)+(1-s1)*As(:,:,k))\(s1*bs(:,j)+(1-s1)*bs(:,k));
		x2 = (s2*As(:,:,j)+(1-s2)*As(:,:,k))\(s2*bs(:,j)+(1-s2)*bs(:,k));
		v1 = (s1*As(:,:,j)+(1-s1)*As(:,:,k))\(bs(:,j)-bs(:,k)-(As(:,:,j)-As(:,:,k))*x1);
		v2 = (s2*As(:,:,j)+(1-s2)*As(:,:,k))\(bs(:,j)-bs(:,k)-(As(:,:,j)-As(:,:,k))*x2);
		f1 = (bs(:,j)-As(:,:,j)*x1)'*Jm*v1;
		f2 = (bs(:,j)-As(:,:,j)*x2)'*Jm*v2;

		if f1*f2 <= 0 && f1 != 0
			sc = s1 - f1*(s2-s1)/(f2-f1);
			xc = (sc*As(:,:,j)+(1-sc)*As(:,:,k))\(sc*bs(:,j)+(1-sc)*bs(:,k));
			launches = [launches, [xc; j]];			
			launches = [launches, [xc; k]];			
		end
	end
end
end

% Find trajectories corresponding to launch points

sym_traj = zeros(2,grid_size2+3,size(launches,2));

for j = 1:size(launches,2);
	x0 = launches(1:2,j);
	kl = launches(3,j);
	sym_traj(:,1:grid_size2+2,j) = trajec(x0(1),x0(2),As(:,:,kl),bs(:,kl),grid_size2,0);
	sym_traj(:,grid_size2+3,j) = kl*[1;1];
end


% Plot if so desired

if (plot_switch)
	hold on
	for j = 1:30
		plot(es_traj(1,1:grid_size2+2,j),es_traj(2,1:grid_size2+2,j),'g')
	end
	for j = 1:size(launches,2)
		plot(sym_traj(1,1:grid_size2+2,j),sym_traj(2,1:grid_size2+2,j),'g')
	end
	hold off
end

% Combine all trajectories into one array and return

all_traj = zeros(2,grid_size2+3,30 + size(launches,2));
all_traj(:,:,1:30) = es_traj;
all_traj(:,:,31:30+size(launches,2)) = sym_traj;

end
