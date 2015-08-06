%% This function generates a random collection of Lindblad matrices for n=3.
%% Each matrix element is complex whose real and imaginary parts are uniform on
%% [0,N_k] with exp(N_k) uniform on [-5,5]. The traceless part is removed.

function [Ls] = rand_Ls(K)

% Initialize

Ls = zeros(3,3,K);

% Generate random operators

for k = 1:K
	Ls(:,:,k) = exp(10*rand-5)*(rand(3)+1i*rand(3));
	Ls(:,:,k) -= trace(Ls(:,:,k))*eye(3)/3;
	Ls(:,:,k) *= exp(-1i*arg(Ls(1,1,k)));
end

end
