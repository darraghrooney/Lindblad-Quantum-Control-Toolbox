% This function takes a b-vector and A-matrix, and an initial orbit, and
% calculates a trajectory in orbit space s.t. xdot = b - A*x.

function[traj] = trajec(x1i, x2i, A, b, tgrid, plot_switch)

% Fetch eigeninformation
[vecs,vals]=eig(A);

% Initialize trajectory
traj=zeros(2,tgrid+2);
rinf = A\b;

% Sweep over points in trajectory
for j=1:tgrid+1
	
	% Set points in time-space	
    	t=(j-1)/tgrid*5/real(vals(2,2));

	% Find trajectories
    	traj(:,j) = rinf + expm(-A*t)*([x1i;x2i]-rinf);
end

% Set final point in trajectory
traj(:,tgrid+2) = rinf;

% Plot trajectories
if plot_switch
	hold on
	plot(traj(1,:),traj(2,:));
	hold off
end

end
