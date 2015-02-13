% This function produces a randomized matrix of Lindblad rates. Only the
% off-diagonal elements are non-zero as rates in the form w_jj do not
% contribute to dynamics. For convenience, the rates have mean one and are
% uniform on [0,2]. If a second integer is provided, it is an indication to set
% the first m columns to zero as they correspond to sink states.

function [w] = rand_w(n,varargin)

	% Produce random rates
	w = rand(n)*2;
	w = w - diag(diag(w));

	% Set sink rates to zero
	if (length(varargin) > 0)
		w(:,1:varargin{1}) *= 0;
	end
end
