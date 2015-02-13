% This function plots the 2-simplex for n=3, plus the inner boundaries.

function[] = arenas()

hold on

% Set axes
axis([-1.1,1.1,-1.3,0.7])

% Plot outer edges
plot([1,-1,0,1],[1/sqrt(3),1/sqrt(3),-2/sqrt(3),1/sqrt(3)],"k");

% Plot inner edges
plot([0,0], [-2/sqrt(3),1/sqrt(3)],"k");
plot([-1/2,1], [-1/2/sqrt(3),1/sqrt(3)],"k");
plot([1/2,-1], [-1/2/sqrt(3),1/sqrt(3)],"k");

hold off	

end
