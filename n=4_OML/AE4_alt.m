% This is an unabridged version of AE4.m. That function is abridged in 
% that it only plots 92 of the total 2024 surfaces. This function plots
% all 2024. This tends to slow down the computer far too much, but it is put
% here for completeness.

function [] = AE4_alt(w,grid_size)

% Check that the co-efficients are non-negative

if ( max(max(w < -1E-15)) == 1)
	error('Co-efficients must be non-negative');
end

% Calculate input rates and re-order if necessary

[w,Jv] = order_rates(w);

% Calculate A(w_ij)
 
[A,b,rinf] = A_b_rinf(w);

% Calculate the images of A(w_ij) under action of the symmetric group

AE = zeros(3,3,24);
bE = zeros(3,24);
AE(:,:,1) = A;
bE(:,1) = b;

S4g = S4gen();

for j=2:24
	AE(:,:,j) = S4g(:,:,j)*A*S4g(:,:,j)'; 
	bE(:,j) = S4g(:,:,j)*b;
end

% Initialization

[s,t]=meshgrid((0:grid_size)/grid_size);

figure
hold on

% Sweep over the choices of triples in Z_24
for j = 1:22
for k = (j+1):23
for l = (k+1):24

% Set mesh to zero
tomesh=zeros(grid_size+1,grid_size+1,3);

% Sweep over grid
for jm = 1:(grid_size+1)
for km = 1:(grid_size+2-jm)

	% Calculate surface point
	invA = s(jm,km)*AE(:,:,k) + t(jm,km)*AE(:,:,l) + (1-s(jm,km)-t(jm,km))*AE(:,:,j);
	bb = s(jm,km)*bE(:,k) + t(jm,km)*bE(:,l) + (1-s(jm,km)-t(jm,km))*bE(:,j);
	tomesh(jm,km,:) = invA\bb;
	
end
end

% "Fold" the triangular mesh over so that Octave can plot it
tomesh = tomesh(:,(grid_size+1):-1:1,:);
for n = 1:3	
	tomesh(:,:,n) = tomesh(:,:,n) + (tomesh(:,:,n) - diag(diag(tomesh(:,:,n))))';
end

% Plot the surface
mesh(tomesh(:,:,1),tomesh(:,:,2),tomesh(:,:,3));	

end
end
end

% Set plot window
extent = norm(A\b);
axis([-extent,extent,-extent,extent,-extent,extent])
view(3)
hold off

end

