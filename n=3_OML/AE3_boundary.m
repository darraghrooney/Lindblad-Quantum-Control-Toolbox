% This function calculates the outer boundary of the STLC set for an n=3 Lindblad system 
% with a discrete control set. The function AE3_arcs finds the 15 candidate arcs, and
% this function determines the boundary of the interior of the region formed from those 
% arcs.	Two grid sizes are required: one for the candidate arcs, and one for the final
% boundary. If plot_switch is 1, the boundary is plotted.

function [boundary, boundaryE] = AE3_boundary(w, grid_size1, grid_size2, plot_switch)

% Fetch candidate arcs 
curves = AE3_arcs(w,grid_size1,0);

% Initialize boundary. First row is distance from origin, second row is the cw angle from
% the +y-axis, third row represents the candidate arc (1 through 15).
boundary = zeros(3,grid_size2+1);
boundary(2,:) = (0:grid_size2)/grid_size2*pi/3;

% Sweep through candidate arcs
for j = 1:15

	% Sweep across pairs of neighboring points on the arc
	for m=1:grid_size1

		% Get info for the pair
	        x1 = curves(m,1,j);
		y1 = curves(m,2,j);
		x2 = curves(m+1,1,j);
		y2 = curves(m+1,2,j);
		s1 = norm([x1,y1]);
		s2 = norm([x2,y2]);

		% If both points are in the first sector and away from the origin only
		if (x2 >= 0 && x1 >= 0 && sqrt(3)*y2 - x2 >= 0 &&...
			sqrt(3)*y1 - x1 >= 0 && s1 >= 1e-10 && s2 >= 1e-10 ) 

			% Get angle information for the points
			th1 = atan(x1/y1);
			th2 = atan(x2/y2);
			thno1 = round(th1*grid_size2*3/pi)+1;
			thno2 = round(th2*grid_size2*3/pi)+1;
	
			% If both points have the same discretized angle
			if (thno1 == thno2)
				if (max(s1,s2) > boundary(1,thno1))
					boundary(1,thno1) = max(s1,s2);
					boundary(3,thno1) = j;
				end

			% If points' discretized angles are off by more than one, we need
			% to extrapolate
			else

				% Calculate extrapolation
				extrap = s1 + (s2 - s1)*(0:abs(thno1-thno2))/abs(thno1 - thno2);

				% Sweep across extrapolation and place on boundary if applicable
				for k = 1:abs(thno2-thno1)+1
					thnok = thno1 + (k-1)*sign(thno2-thno1);
					if (extrap(k) > boundary(1,thnok))
						boundary(1,thnok) = extrap(k);
						boundary(3,thnok) = j;
					end
				end
			end

		% If the first point is in sector I, but the other is in sector II 
		elseif (x2 < 0 && x1 >= 0 && sqrt(3)*y2 + x2 >= 0 &&...
			sqrt(3)*y1 - x1 >= 0 && s1 >= 1e-10 && s2 >= 1e-10 ) 
			
			% Get angle information for the first point
			th1 = atan(x1/y1);
			thno1 = round(th1*grid_size2*3/pi)+1;
			
			% Extrapolate between point and the y-axis
			intercept = y1 - x1*(y1-y2)/(x1-x2);
			extrap = intercept + (s1-intercept)*(0:thno1-1)/thno1;			
	
			% Sweep over angles
			for k = 1:thno1

				% Add point to boundary if applicable
				if (extrap(k) > boundary(1,k))
					boundary(1,k) = extrap(k);
					boundary(3,k) = j;
				end
			end

		% If the second point is in sector I, but the other is in sector II 
		elseif (x2 >= 0 && x1 < 0 && sqrt(3)*y2 - x2 >= 0 &&...
			sqrt(3)*y1 + x1 >= 0 && s1 >= 1e-10 && s2 >= 1e-10 ) 

			% Get angle information for the second point
			th2 = atan(x2/y2);
			thno2 = round(th2*grid_size2*3/pi)+1;
			
			% Extrapolate between point and the y-axis
			intercept = y1 - x1*(y1-y2)/(x1-x2);
			extrap = intercept + (s2-intercept)*(0:thno2-1)/thno2;			
	
			% Sweep over angles
			for k = 1:thno2

				% Add point to boundary if applicable
				if (extrap(k) > boundary(1,k))
					boundary(1,k) = extrap(k);
					boundary(3,k) = j;
				end
			end

		% If the first point is in sector I, but the other is in sector VI 
		elseif (x2 >= 0 && x1 >= 0 && sqrt(3)*y2 - x2 < 0 &&...
			sqrt(3)*y1 - x1 >= 0 && s1 >= 1e-10 && s2 >= 1e-10 ) 

			% Get angle information for the first point
			th1 = atan(x1/y1);
			thno1 = round(th1*grid_size2*3/pi)+1;

			% Extrapolate between point and the y-axis			
			m = (y1-y2)/(x1-x2);
			intercept = 2*(y1 - m*x1)/(1 - m*sqrt(3)); 
			extrap = intercept + (s1 - intercept)*...
				(0:(grid_size2 + 1 - thno1))/(grid_size2 + 2 - thno1);			
	
			% Sweep over angles
			for k = grid_size2 + 1:-1:thno1

				% Add point to boundary if applicable
				if (extrap(grid_size2 + 2 - k) > boundary(1,k))
					boundary(1,k) = extrap(grid_size2 + 2 - k);
					boundary(3,k) = j;
				end
			end

		% If the second point is in sector I, but the other is in sector VI 
		elseif (x2 >= 0 && x1 >= 0 && sqrt(3)*y2 - x2 >= 0 &&...
			sqrt(3)*y1 - x1 < 0 && s1 >= 1e-10 && s2 >= 1e-10 ) 

			% Get angle information for the second point
			th2 = atan(x2/y2);
			thno2 = round(th2*grid_size2*3/pi)+1;

			% Extrapolate between point and the y-axis			
			m = (y1-y2)/(x1-x2);
			intercept = 2*(y1 - m*x1)/(1 - m*sqrt(3)); 
			extrap = intercept + (s2 - intercept)*...
				(0:(grid_size2 + 1 - thno2))/(grid_size2 + 2 - thno2);			
	
			% Sweep over angles
			for k = grid_size2 + 1:-1:thno2

				% Add point to boundary if applicable
				if (extrap(grid_size2 + 2 - k) > boundary(1,k))
					boundary(1,k) = extrap(grid_size2 + 2 - k);
					boundary(3,k) = j;
				end
			end

		end
	end
end

% Create an extended boundary for all six sectors
boundaryE = [boundary(1,:).*sin(boundary(2,:)); boundary(1,:).*cos(boundary(2,:))];

% Reflect across lamb2 = lamb3 line
boundaryE = [boundaryE, [boundary(1,grid_size2+1:-1:1).*sin(2*pi/3-boundary(2,grid_size2+1:-1:1)); ...
		boundary(1,grid_size2+1:-1:1).*cos(2*pi/3-boundary(2,grid_size2+1:-1:1))]];

% Rotate original by 120 degrees
boundaryE = [boundaryE, [boundary(1,:).*sin(boundary(2,:) + 2*pi/3); ...
		boundary(1,:).*cos(boundary(2,:) + 2*pi/3)]];

% Reflect original across lamb1 = lamb3 line
boundaryE = [boundaryE, [boundary(1,grid_size2+1:-1:1).*sin(-2*pi/3-boundary(2,grid_size2+1:-1:1)); ...
		boundary(1,grid_size2+1:-1:1).*cos(-2*pi/3-boundary(2,grid_size2+1:-1:1))]];

% Rotate original by 240 degrees
boundaryE = [boundaryE, [boundary(1,:).*sin(boundary(2,:) + 4*pi/3); ...
		boundary(1,:).*cos(boundary(2,:) + 4*pi/3)]];

% Reflect original across lamb1 = lamb2 line
boundaryE = [boundaryE, [boundary(1,grid_size2+1:-1:1).*sin(-boundary(2,grid_size2+1:-1:1)); ...
		boundary(1,grid_size2+1:-1:1).*cos(boundary(2,grid_size2+1:-1:1))]];

% Plot extended boundary if requested
hold on
if (plot_switch)
	plot(boundaryE(1,:),boundaryE(2,:),'r')
end
hold off

end
