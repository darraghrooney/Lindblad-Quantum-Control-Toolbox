% This function calculates/plots the boundaries for the SLC set for
% an n=3 Lindblad system, for a discrete set of flags. The identity flag
% and its 5 permutations are automatically plotted. A specified number of 
% random flags, as well as their permutations are selected, for a total of
% 6*un_no. The total number of curves is 6*(un_no+1)*(6*(un_no+1)-1)/2. 
% If plot_switch is zero, the plot is suppressed and only the curves are 
% returned.

function [curves] = AE3_arcs(L, grid_size, un_no, plot_switch)

% Initialization

bs = zeros(2, 6*(un_no+1));
As = zeros(2, 2, 6*(un_no+1));
flag = eye(3);
count = 0;
gens = S3gen();

% Loop over number of random flags.

while (count <= un_no)

	% Calculate w, b, A for current flag

	w = Ltow(L, flag);
	[A,b] = A_b_rinf(w);

	% Calculate b, A for permutations	
	
	for j = 1:6
		bs(:,6*count+j) = gens(:,:,j)*b;
		As(:,:,6*count+j) = gens(:,:,j)*A*gens(:,:,j)';
	end

	% Calculate next random flag	

	flag = flag_gen([],1);

count++;
end

% Initialize curves

pair_count = 6*(un_no+1)*(6*(un_no+1)-1)/2;
curves = zeros(grid_size+1,2,pair_count);
count = 1;

% Sweep over pairs of vertices

for j = 1:(6*un_no+6)
for k = (j+1):(6*un_no+6)

	% Sweep over grid
	for l = 0:grid_size
		s = l/grid_size;
		Ain = (1-s)*As(:,:,j) + s*As(:,:,k);
		bin = (1-s)*bs(:,j) + s*bs(:,k);
		curves(l+1,:,count) = Ain\bin; 
	end
	count += 1;
end
end

% Plot the curves

if (plot_switch)
	hold on
	for j = 1:pair_count
		plot(curves(:,1,j),curves(:,2,j),'b');
	end


	% Adapt plot window

	extent0 = max(max(curves(:,1,:)));
	extent1 = max(max(curves(:,2,:)));
	extent2 = min(min(curves(:,2,:)));

	axis([-1.1*extent0,1.1*extent0,1.1*extent2,1.1*extent1]);

	hold off
end

end
