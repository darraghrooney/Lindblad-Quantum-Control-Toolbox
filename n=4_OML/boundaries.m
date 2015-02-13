% This function plots the edges of the arena for n=4. i.e. the edges of a 
% regular tetrahedron with side length = two.

function[] = boundaries()

hold on

% Set axes
axis([-1.1,1.1,-1.3,0.7,-1.3,0.6])

% Plot top edges
plot3([1,-1,0,1],[1/sqrt(3),1/sqrt(3),-2/sqrt(3),1/sqrt(3)], ...
	[1/sqrt(6),1/sqrt(6),1/sqrt(6),1/sqrt(6)],"k");

% Plot descending edges
plot3([0,1],[0,1/sqrt(3)],[-3/sqrt(6),1/sqrt(6)],"k");
plot3([0,-1],[0,1/sqrt(3)],[-3/sqrt(6),1/sqrt(6)],"k");
plot3([0,0],[0,-2/sqrt(3)],[-3/sqrt(6),1/sqrt(6)],"k");

hold off	

end
