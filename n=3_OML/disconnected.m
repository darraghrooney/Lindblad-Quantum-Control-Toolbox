% This function plots the STLC set and the GC set for a degenerate n=2 system, 
% where there are no rates into or out of one of the three nodes. WLOG, this is taken
% to be node 3, so only w(1,2) and w(2,1) are non-zero. The GC set is plotted in green,
% and the STLC set, which is smaller, is super-imposed in blue.

function [ ] = disconnected(w12,w21)

% Plot the simplex boundaries
arenas;

hold on

% If w12 and w21 are equal, the STLC and GC sets are just the center point, 
% so don't bother plotting
if (w12==w21)
	error('GC set is only a point for w12 = w21');
else

	% Useful constant
	R = (w12+w21)/(w12-w21);

	
	% Plot the GC set
	fill(4*[(R+1)/2,R,(R-1)/2,-(R-1)/2,-R,-(R+1)/2,(R+1)/2]/(1+3*R^2), ...
		4*[(3*R-1)/2,1,(-3*R-1)/2,(-3*R-1)/2,1,(3*R-1)/2,(3*R-1)/2]...
		/(1+3*R^2)/sqrt(3),'g');

	% Plot the STLC set
	fill([2*(R+1)/(1+3*R^2),-2/(1-3*R),4*R/(1+3*R^2),2/(1+3*R),2*(R-1)/...
		(1+3*R^2),0,-2*(R-1)/(1+3*R^2),-2/(1+3*R),-4*R/(1+3*R^2),2/...
		(1-3*R),-2*(R+1)/(1+3*R^2),0,2*(R+1)/(1+3*R^2)],[2*(3*R-1)/...
		(1+3*R^2)/sqrt(3),-2/(1-3*R)/sqrt(3),4/(1+3*R^2)/sqrt(3),-2/...
		(1+3*R)/sqrt(3),2*(-3*R-1)/(1+3*R^2)/sqrt(3),4/(1-3*R)/sqrt(3),...
		2*(-3*R-1)/(1+3*R^2)/sqrt(3),-2/(1+3*R)/sqrt(3),4/(1+3*R^2)/...
		sqrt(3),-2/(1-3*R)/sqrt(3),2*(3*R-1)/(1+3*R^2)/sqrt(3),4/(1+3*R)...
		/sqrt(3),2*(3*R-1)/(1+3*R^2)/sqrt(3)],'b');


	% Plot the six lines that make up the set boundaries
	plot([1,0,-1]/R,[1,-2,1]/sqrt(3),'b');
	plot([1/2-1/2/R,-1,1/2+1/2/R],[-1/2-3/2/R,1,-1/2+3/2/R]/sqrt(3),'b');
	plot(-[1/2-1/2/R,-1,1/2+1/2/R],[-1/2-3/2/R,1,-1/2+3/2/R]/sqrt(3),'b');

end

hold off

end

