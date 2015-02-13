% For a Lindblad system with m1>0 sinks, there are m1 constraints on the
% eigenvalues (including the constraint they add to one). This function
% finds these constraints. The constraint matrix is such that the inner product
% of any column, with any column of B, is zero. Where d(lambda)/dt = B*lambda

% Such a matrix is not unique, but this one is chosen so that the upper m1 rows
% form a multiple of the identity matrix. This multiple happens to be the sum of
% products w(j_1,k_1)...w(j_n,k_n) where the pairs (j,k) form a multi-tree on the 
% nodes 1...n. Specifically, the multi-tree is m1 disjoint rooted trees, where
% the sinks are the roots, and the descendents are the non-sinks.

% The other elements of the constraint matrix are the sum of such products, but 
% obeying the following selection rule: the (j,k) element is formed only from
% multi-trees where the node j is in the tree with root k.

function[constraints,B] = multisink(m1,m2,w)

	% Check that m1 and m2 add up to dim(w)
	if (size(w,1) != m1+m2 && size(w,2) != m1+m2)
		error("Dimension mismatch");
	end

	% Check that rates from sinks are zero
	if (norm(w(:,1:m1)) != 0)
		disp("w(j,k) should be zero for k <= m1");	
		disp("Now setting these rates to zero");
		w(:,1:m1) *= 0;
	end


	B = B_mat(w);

	% Initializing constraint coefficients
	constraints = zeros(m1+m2,m1);

	% Use a separate variable for rates from all trees
	alltrees = 0;

	% 'assignment' will keep track of the leaf partition
	assignment = zeros(m2,1);

	% Sum over the possible partitions
	for j = 1:m1^m2

		% Initialize the contribution from this partition
		toadd = 1;
	
		% Sweep over the sinks
		for k = 1:m1

			% Figure out what leaves belong to this sink
			leaves = [];
			for l = 1:m2	
				if (assignment(l) == k-1)
					leaves = [leaves,l];
				end 
			end

			% Apply the product for this sink 
			new_prod = rt_w(w,k,leaves+m1);		
			toadd *= new_prod;
		end

		% Add this contribution to the total rate
		alltrees += toadd;

		% Add this contribution to the appropriate particular rates
		for l = 1:m2		
			constraints(m1+l,assignment(l)+1) += toadd;	
		end

		% Increment partition
		assignment = inc_nd(m1,m2,assignment);
	end

	% Add the total rate to the constraint matrix
	constraints(1:m1,:) = alltrees*eye(m1);	

end
