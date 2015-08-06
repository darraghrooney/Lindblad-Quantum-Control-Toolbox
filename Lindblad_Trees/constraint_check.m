% This function checks the constraints on a multi-sink process.
% The inner products of the columns of Omega and the constraint matrix should
% be zero. This function finds the matrix of errors and returns the norm.

function[err_norm] = constraint_check(constraints, Omega)

	% Find the number of sinks and sources
	m1 = size(constraints,2);
	m2 = size(constraints,1)-m1;

	% Check to see B and the constraints match dimension
	if (m1+m2 !=  size(Omega,1))
		error("Dimension mismatch");
	end

	% Check to see the first m1 columns represent sinks
	if (norm(Omega(:,1:m1)) > 1e-10)
		error("Sink rates should be zero");
	end

	% Calculate error matrix and norm
	errors = constraints'*Omega(:,(m1+1):(m1+m2));
	err_norm = norm(errors);

end
