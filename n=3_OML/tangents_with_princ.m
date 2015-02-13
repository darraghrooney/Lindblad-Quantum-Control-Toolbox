% This function plots the discretized range of tangent vectors, given
% a collection of Lindblad operators and an orbit point, and compares 
% them to the six principal tangent vectors. 

function[tngs] = tangents_with_princ(Ls, x1, x2, gno)

% Plot the discretized range of tangent vectors
tangents(Ls, x1, x2, gno);

% Fetch relevant info about the principal controls
[A,b,rinf] = A_b_rinf(principal_w(Ls));

% Get symmetry generators
S3g = S3gen;

hold on

% Sweep through vector fields
for j = 1:6

	% Get tangent vector and plot
	v = S3g(:,:,j)*(b - A*S3g(:,:,j)'*[x1,x2]');	
	plot(v(1),v(2), "or", "markersize", 10);
end
hold off

end

