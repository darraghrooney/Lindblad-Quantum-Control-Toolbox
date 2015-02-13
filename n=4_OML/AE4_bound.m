% This function plots the edges of the STLC set for n=4. Edges meaning the edges
% of the bounding surfaces. Since there are 2024 such surfaces, and this can overwhelm
% the CPU when plotting, a good alternative is to plot the edges. There are 6072 such
% edges, but a curve is much easier to handle than a mesh.

function [] = AE4_bound(w,grid_size)

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


figure
hold on

% Sweep over choices of vertex pairs
for j = 1:23
for k = (j+1):23

% Initialization
toplot=zeros(grid_size+1,3);
s = (0:grid_size)/grid_size;

% Sweep over the grid
for jm = 1:(grid_size+1)

	% Calculate the location of the edge-point
	invA = s(jm)*AE(:,:,k) + (1-s(jm))*AE(:,:,j);
	bb = s(jm)*bE(:,k) + (1-s(jm))*bE(:,j);
	toplot(jm,:) = invA\bb;
	
end

% Plot the edge corresponding to the vertex pair
plot3(toplot(:,1),toplot(:,2),toplot(:,3));	

end
end

% Set plot size
extent = norm(A\b);
axis([-extent,extent,-extent,extent,-extent,extent])
view(3)
hold off

end

