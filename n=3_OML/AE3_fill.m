% This function plots the filled A3E region for the STLC set for
% an n=3 Lindblad system, given a discrete set of six Lindblad 
% rates. 

function [] = AE3_fill(w,grid_size)

% Check that the co-efficients are non-negative

if ( max(max(w < -1E-15)) == 1)
	error('Co-efficients must be non-negative');
end

% Calculate input rates and re-order if necessary

[w,Jv] = order_rates(w);

% Calculate A(w_ij)
 
[A,b,rinf] = A_b_rinf(w);

% Calculate the images of A(w_ij) under action of the symmetric group

AE = zeros(2,2,6);
bE = zeros(2,6);
AE(:,:,1) = A;
bE(:,1) = b;

S3g = S3gen();

for j=2:6
	AE(:,:,j) = S3g(:,:,j)*A*S3g(:,:,j)'; 
	bE(:,j) = S3g(:,:,j)*b;
end

% Initialize curves

curves = zeros(grid_size+1,2,15);

% Used for tracking pair indices
pairs = zeros(6,6);
count = 1;

% Sweep over pairs of vertices

for j = 1:6
for k = (j+1):6

	% Sweep over grid
	for l = 0:grid_size
		s = l/grid_size;
		Ain = (1-s)*AE(:,:,j) + s*AE(:,:,k);
		bin = (1-s)*bE(:,j) + s*bE(:,k);
		curves(l+1,:,count) = Ain\bin; 
	end
	pairs(j,k) = count;
	count += 1;
end
end

% Plot the patches

hold on
for j = 1:4
	for k = j+1:5
		for l = k+1:6
			bnd = [curves(:,:,pairs(j,k))',curves(2:grid_size+1,:,pairs(k,l))',...
				curves(grid_size:-1:1,:,pairs(j,l))'];	
			patch(bnd(1,:),bnd(2,:),'b','edgecolor','b');
		end
	end	
end

% Adapt plot window

extent0 = max(max(curves(:,1,:)));
extent1 = max(max(curves(:,2,:)));
extent2 = min(min(curves(:,2,:)));

axis([-1.1*extent0,1.1*extent0,1.1*extent2,1.1*extent1]);

hold off
end


