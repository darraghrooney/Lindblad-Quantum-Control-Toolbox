%% This function generates a random set of Lindblad rates in n=4
%% Mean is set to 1 ... distribution is uniform on [0,2]

function [w,Jvec] = rand_w()

% Generate random variables
w = 2*rand(4,4);

% Set diagonal entries to zero
w -= diag(diag(w));

% Make sure the rates are ordered properly
[w,Jvec] = order_rates(w);

end
