% This function finds the steady state of a Lindblad system, provided it is 
% connected and has no more than one sink. The un-normalized component in the 
% mth direction is equal to the sum of products w(j_1,k_1)...w(j_?,k_?) where 
% the edges (j_?,k_?) of each product form a tree rooted at j. 

function[lamb_inf] = stat_state(w)

	% Find dimension of system
	n = size(w,1);

	% Check that w is non-negative and non-zero
	if (!(w >= 0) && norm(diag(w)) > 1e-10 &&
		n != size(w,2))
		error("Invalid w");
	end

	% Check that there are no more than one sink
	wsort = sort(diag(w));
	if (!(sort(2)>1e-10))
		error("More than one sink");
	end

	% Compute the un-normalized components
	input_rates = zeros(n,1);
	for j = 1:n
		input_rates(j) = rt_w(w,j,[1:(j-1),(j+1):n]);
	end

	% Normalize the steady-state
	total_rate = sum(input_rates);
	lamb_inf = input_rates/total_rate;
end
