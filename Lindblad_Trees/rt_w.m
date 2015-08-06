% This function takes a matrix of Lindblad rates and returns the sum of 
% products w(j_1,k_1)w(j_2,k_2)...w(j_m,k_m). Each product is constructed by a 
% rooted tree composed of edges (j_i,k_i). Each unique complete rooted tree, with
% root and descendent vertices specified by input, is represented precisely once
% in the sum.

function[K_T] = rt_w(w,root,desc)

	% Find the number of descendents
	lno = length(desc);	

	% If there are no descendents, the product is one
	if (lno==0) K_T = 1;
	
	else
		% Initialize		
		K_T = 0;

		% Sweep over the possible choices of depth-one daughters
		% 'assignment' is a binary number representing possibilities	
		assignment = zeros(lno,1);		
		for j = 1:(2^lno-1)
			
			% These vectors collect daughter and deeper nodes		
			daughters = [];
			deeper = [];		
			
			% Increment before as (0,...,0) is not a possibility
			% but (1,...,1) is
			assignment = inc_nd(2,lno,assignment);

			% Split descendent collection into daughters and others
			for k = 0:(lno-1)
				if (assignment(k+1) == 1)
					daughters = [daughters,desc(k+1)];
				else
					deeper = [deeper,desc(k+1)];
				end
			end

			% Get helper function to calculate the individual product
			% Add this product to sum of products
			K_T += K_helper(w,root,daughters,deeper);
		end
	end
end

% This function for a given root, daughters, and deeper nodes, computes
% the corresponding products w(j_1,k_1)...w(j_?,k_?) and sums. Products are summed.

function[K_T_dau] = K_helper(w,root,daughters,deeper)

	% Initialize part corresponding to daughters
	K_T_dau = 1;

	% Find number of daughters and deeper nodes
	m1 = length(daughters);	
	m2 = length(deeper);	

	% Multiply the rates corresponding to daughters
	for j = 1:m1
		K_T_dau *= w(root,daughters(j));
	end

	% Initialize the part corresponding to deeper nodes
	K_T_deep = 0;

	% Sweep over possible assignments of deeper nodes to daughters
	assignment = zeros(m2,1);
	for j = 1:(m1^m2)

		% Initialize the product for each partition
		K_T_part = 1;	
	
		% Sweep over sinks
		for k = 1:m1

			% Find the deeper nodes in sink's component 
			k_nodes = [];
			for l = 0:(m2-1)
				if (assignment(l+1) == k-1)
					k_nodes = [k_nodes,deeper(l+1)];	
				end		
			end
		
			% Multiply contribution from sink's component
			K_T_part *= rt_w(w,daughters(k),k_nodes);
		end

		% Add the contribution from this partition
		K_T_deep += K_T_part;

		% Increment to find new partition
		assignment = inc_nd(m1,m2,assignment);
	end

	% Multiply the part for the daughters by the part for the deeper nodes
	K_T_dau *= K_T_deep;
end

