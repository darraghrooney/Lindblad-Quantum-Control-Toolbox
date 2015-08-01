% This function plots the edges of the arena for n=4. i.e. the edges of a 
% regular tetrahedron with side length = sqrt(2).

function[] = boundaries()

hold on

% Set axes
axis([-1.1,1.1,-1.3,0.7,-1.3,0.6]/sqrt(2))

% Plot top edges
plot3([1,-1,0,1]/sqrt(2),[1/sqrt(3),1/sqrt(3),-2/sqrt(3),1/sqrt(3)]/sqrt(2), ...
	[1/sqrt(6),1/sqrt(6),1/sqrt(6),1/sqrt(6)]/sqrt(2),"k");

% Plot descending edges
plot3([0,1]/sqrt(2),[0,1/sqrt(3)]/sqrt(2),[-3/sqrt(6),1/sqrt(6)]/sqrt(2),"k");
plot3([0,-1]/sqrt(2),[0,1/sqrt(3)/sqrt(2)],[-3/sqrt(6),1/sqrt(6)]/sqrt(2),"k");
plot3([0,0],[0,-2/sqrt(3)/sqrt(2)],[-3/sqrt(6),1/sqrt(6)]/sqrt(2),"k");
box off
hold off	

end
