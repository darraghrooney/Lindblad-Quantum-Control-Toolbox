%% This function generates a random GKS matrix for n=3.
%% A randomized Hermitian matrix is produced (all elements uniform on
%% [0,2]), and then any negative eigenvalues are flipped.

function [GKS] = rand_GKS()

% Generate random Hermitian matrix

GKS = rand(8)+1i*rand(8);
GKS = GKS + GKS';

% Get eigeninfo

[vecs,vals] = eig(GKS);

% Flip negative eigenvalues

vals = abs(vals);

% Reconstruct GKS matrix

GKS = zeros(8);
for k = 1:8
	GKS += vals(k,k)*vecs(:,k)*vecs(:,k)';
end	

end
