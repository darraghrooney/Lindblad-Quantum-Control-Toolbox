% This function calculates the outer boundary of the B03 set for an n=3 Lindblad system 
% with a discrete control set. The trajectories are loaded, and this function determines 
% the boundary of the interior of the region formed from them. Three grid sizes are required: 
% one for the AE3 arcs, one for the B03 trajectories and one for the boundary. 
% If plot_switch is 1, the boundary is plotted.

function [boundaryE] = B03_boundary(w, grid_size1, grid_size2, grid_size3, plot_switch)

% Fetch B03 trajectories
traj = B03_arcs(w,grid_size1,grid_size2,0);
t_no = size(traj,3);


% Initialize boundary. First row is distance from origin, second row is the cw angle from
% the +y-axis, third row represents the appropriate vector field (1 thru 6).
boundary = zeros(3,grid_size3+1);
boundary(2,:) = (0:grid_size3)/grid_size3*pi/3;

% Sweep through trajectories
for j = 1:t_no

	% Sweep across pairs of neighboring points on the trajectory
	for m=1:grid_size2+1

		% Get info for the pair
	        x1 = traj(1,m,j);
		y1 = traj(2,m,j);
		x2 = traj(1,m+1,j);
		y2 = traj(2,m+1,j);
		s1 = norm([x1,y1]);
		s2 = norm([x2,y2]);
		
		% If both points are in the first sector and away from the origin only
		if (x2 >= 0 && x1 >= 0 && sqrt(3)*y2 - x2 >= 0 &&...
			sqrt(3)*y1 - x1 >= 0 && s1 >= 1e-10 && s2 >= 1e-10 ) 
	
			% Get angle information for the points
			th1 = atan(x1/y1);
			th2 = atan(x2/y2);
			thno1 = round(th1*grid_size3*3/pi)+1;
			thno2 = round(th2*grid_size3*3/pi)+1;
	
			% If both points have the same discretized angle
			if (thno1 == thno2)
				if (max(s1,s2) > boundary(1,thno1))
					boundary(1,thno1) = max(s1,s2);
					boundary(3,thno1) = traj(1,grid_size2+3,j);
				end

			% If points' discretized angles are off by more than one, we need
			% to interpolate
			else

				% Sweep between angles and place on boundary if applicable
				for k = 1:abs(thno2-thno1)+1
					thnok = thno1 + (k-1)*sign(thno2-thno1);
					thk = (thnok-1)/grid_size3*pi/3;
					sk = (y1*(x2-x1)-x1*(y2-y1))/...
						(cos(thk)*(x2-x1)-sin(thk)*(y2-y1));
					if (sk > boundary(1,thnok))
						boundary(1,thnok) = sk;
						boundary(3,thnok) = traj(1,grid_size2+3,j);
					end
				end
			end

		% If the first point is in sector I, but the other is in sector II 
		elseif (x2 < 0 && x1 >= 0 && sqrt(3)*y2 + x2 >= 0 &&...
			sqrt(3)*y1 - x1 >= 0 && s1 >= 1e-10 && s2 >= 1e-10 ) 
			
			% Get angle information for the first point
			th1 = atan(x1/y1);
			thno1 = round(th1*grid_size3*3/pi)+1;
			
			% Calculate intercept of interpolation
			intercept = y1 - x1*(y1-y2)/(x1-x2);
	
			% Sweep over angles and interpolate
			for k = 1:thno1
				thk = (k-1)/grid_size3*pi/3;
				sk = (y1*(-x1)-x1*(intercept-y1))/...
						(cos(thk)*(-x1)-sin(thk)*(intercept-y1));

				% Add point to boundary if applicable
				if (sk > boundary(1,k))
						boundary(1,k) = sk;
						boundary(3,k) = traj(1,grid_size2+3,j);
				end

			end

		% If the second point is in sector I, but the other is in sector II 
		elseif (x2 >= 0 && x1 < 0 && sqrt(3)*y2 - x2 >= 0 &&...
			sqrt(3)*y1 + x1 >= 0 && s1 >= 1e-10 && s2 >= 1e-10 ) 

			% Get angle information for the second point
			th2 = atan(x2/y2);
			thno2 = round(th2*grid_size3*3/pi)+1;
			
			% Calculate intercept
			intercept = y1 - x1*(y1-y2)/(x1-x2);
	
			% Sweep over angles and interpolate
			for k = 1:thno2
				thk = (k-1)/grid_size3*pi/3;
				sk = (y2*(-x2)-x2*(intercept-y2))/...
						(cos(thk)*(-x2)-sin(thk)*(intercept-y2));

				% Add point to boundary if applicable
				if (sk > boundary(1,k))
						boundary(1,k) = sk;
						boundary(3,k) = traj(1,grid_size2+3,j);
				end
			end

		% If the first point is in sector I, but the other is in sector VI 
		elseif (x2 >= 0 && x1 >= 0 && sqrt(3)*y2 - x2 < 0 &&...
			sqrt(3)*y1 - x1 >= 0 && s1 >= 1e-10 && s2 >= 1e-10 ) 

			% Get angle information for the first point
			th1 = atan(x1/y1);
			thno1 = round(th1*grid_size3*3/pi)+1;

			% Calculate intercept
			m = (y1-y2)/(x1-x2);
			intercept = 2*(y1 - m*x1)/(1 - m*sqrt(3)); 
	
			% Sweep over angles and interpolate
			for k = grid_size3 + 1:-1:thno1

				thk = (k-1)/grid_size3*pi/3;
				sk = (y1*(intercept*sqrt(3)/2-x1)-x1*(intercept/2-y1))/...
					(cos(thk)*(intercept*sqrt(3)/2-x1)...
						-sin(thk)*(intercept/2-y1));

				% Add point to boundary if applicable
				if (sk > boundary(1,k))
						boundary(1,k) = sk;
						boundary(3,k) = traj(1,grid_size2+3,j);
				end
			end

		% If the second point is in sector I, but the other is in sector VI 
		elseif (x2 >= 0 && x1 >= 0 && sqrt(3)*y2 - x2 >= 0 &&...
			sqrt(3)*y1 - x1 < 0 && s1 >= 1e-10 && s2 >= 1e-10 ) 

			% Get angle information for the second point
			th2 = atan(x2/y2);
			thno2 = round(th2*grid_size3*3/pi)+1;

			% Calculate intercept			
			m = (y1-y2)/(x1-x2);
			intercept = 2*(y1 - m*x1)/(1 - m*sqrt(3)); 
	
			% Sweep over angles and interpolate
			for k = grid_size3 + 1:-1:thno2

				thk = (k-1)/grid_size3*pi/3;
				sk = (y2*(intercept*sqrt(3)/2-x2)-x2*(intercept/2-y2))/...
					(cos(thk)*(intercept*sqrt(3)/2-x2)...
						-sin(thk)*(intercept/2-y2));

				% Add point to boundary if applicable
				if (sk > boundary(1,k))
						boundary(1,k) = sk;
						boundary(3,k) = traj(1,grid_size2+3,j);
				end

			end

		end
	end
end

% This is the multiplication table for S3
j_trans = [[1,2,3,4,5,6];  [2,1,5,6,3,4]; [3,6,1,5,4,2]; [4,5,6,1,2,3,];...
	[5,4,2,3,6,1]; [6,3,4,2,1,5]];

% Create an extended boundary for all six sectors
boundaryE = [boundary(1,:).*sin(boundary(2,:)); ...
	boundary(1,:).*cos(boundary(2,:));boundary(3,:)];

% Reflect across lamb2 = lamb3 line
toadd = [boundary(1,grid_size3+1:-1:1).*sin(2*pi/3-boundary(2,grid_size3+1:-1:1)); ...
	boundary(1,grid_size3+1:-1:1).*cos(2*pi/3-boundary(2,grid_size3+1:-1:1));...
		zeros(1,grid_size3+1)];

for j = 1:grid_size3+1
	toadd(3,grid_size3+2-j) = j_trans(boundary(3,j),4);
end

boundaryE = [boundaryE, toadd];

% Rotate original by 120 degrees

toadd = [boundary(1,:).*sin(boundary(2,:) + 2*pi/3); ...
	boundary(1,:).*cos(boundary(2,:) + 2*pi/3); zeros(1,grid_size3+1)];

for j = 1:grid_size3+1
	toadd(3,j) = j_trans(boundary(3,j),6);
end

boundaryE = [boundaryE, toadd];

% Reflect original across lamb1 = lamb3 line

toadd = [boundary(1,grid_size3+1:-1:1).*sin(-2*pi/3-boundary(2,grid_size3+1:-1:1)); ...
	boundary(1,grid_size3+1:-1:1).*cos(-2*pi/3-boundary(2,grid_size3+1:-1:1)); ...
		zeros(1,grid_size3+1)];

for j = 1:grid_size3+1
	toadd(3,grid_size3+2-j) = j_trans(boundary(3,j),3);
end

boundaryE = [boundaryE, toadd];

% Rotate original by 240 degrees

toadd = [boundary(1,:).*sin(boundary(2,:) + 4*pi/3); ...
	boundary(1,:).*cos(boundary(2,:) + 4*pi/3); zeros(1,grid_size3+1)];

for j = 1:grid_size3+1
	toadd(3,j) = j_trans(boundary(3,j),5);
end

boundaryE = [boundaryE, toadd];

% Reflect original across lamb1 = lamb2 line

toadd = [boundary(1,grid_size3+1:-1:1).*sin(-boundary(2,grid_size3+1:-1:1)); ...
	boundary(1,grid_size3+1:-1:1).*cos(-boundary(2,grid_size3+1:-1:1)); ...
		zeros(1,grid_size3+1)];

for j = 1:grid_size3+1
	toadd(3,grid_size3+2-j) = j_trans(boundary(3,j),2);
end

boundaryE = [boundaryE, toadd];

% Plot extended boundary if requested

hold on
if (plot_switch)
	plot(boundaryE(1,:),boundaryE(2,:),'m')
end
hold off

end
