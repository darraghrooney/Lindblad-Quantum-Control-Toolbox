%% This take a curve in R^2 and transforms it using a generator of 
%% the symmetric group S3. "r" specifies which generator to use.

function[rcurve] = rot_curve(curve,r)

% Initialization
S3g = S3gen();
rcurve = 0*curve;

% Sweep over points in the curve
for j=1:size(curve(:,:,1),1)

		% Use the specified generator to rotate the given vector
		for l = 1:2
			inp = S3g(:,:,r)*vec(curve(j,:));
			rcurve(j,l)=inp(l);
		end
	end
end

