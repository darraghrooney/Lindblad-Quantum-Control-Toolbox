% For a set of Lindblad rates, this function takes points on the boundary
% of the primary arena, and plots the corresponding trajectories.

function[] = trajectories(w,gno)

% Get A and b
[A,b,trash] = A_b_rinf(w);

% Sweep over initial points
for j = 0:gno
	
	% Plot points from the left boundary
	trajec(j/gno,1/sqrt(3),A,b,100,1);

	% Plot points from the top boundary
	trajec(0,j/gno/sqrt(3),A,b,100,1);

	% Plot points from the diagonal boundary
	trajec(j/gno,j/gno/sqrt(3),A,b,100,1);
end

% Plot the boundary
hold on
plot([0,0,1,0],[0,1/sqrt(3),1/sqrt(3),0],"k")
hold off

end
