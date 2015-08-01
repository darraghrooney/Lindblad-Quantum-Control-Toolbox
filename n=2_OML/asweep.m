% This function sweeps through values of aj to measure various parameters
% of the distribution of rcrit over a surface in b-space. Namely the level surface
% b1/2/sqrt(a2 a3) + b2/2/sqrt(a3 a1) + b3/2/sqrt(a1 a2) = B. 

% Examination of these distributions allows us to select 4 relevant features:

% 1. There are equal local minima at (phi,theta)=(pi/4,pi/2) and (pi/4,3pi/2). 
%    Parameter p1 will measure the rcrit of this minima.
% 2. There are equal saddle points at (phi,theta)=(pi/4,0) and (pi/4,pi).
%    Parameter p2 will measure the rcrit at these saddle points.
% 3. There are saddle points at the poles phi=0,pi/2.
%    Parameter p3 will measure the rcrit there.
% 4. There are pairs of local maxima on the longitudes theta=0,pi. The rcrit of 
%    these maxima is always equal to B (for unknown reasons). The pair are symmetric
%    about the equator so we will set parameter p4 to measure the cos(phi) value of the 
%    northern element (on the interval [0,1]).

function[p1,p2,p3,p4] = asweep(B,samp1,samp2,samp3)

% Define values of aj's
a1 = 100;
a2 = (1:samp2)*100/samp2;
a3 = (0:samp2-1)*100/samp2;

% Initialize parameters to -1 so that we can cut off unassigned slots.
p1=-ones(samp2,samp2);
p2=-ones(samp2,samp2);
p3=-ones(samp2,samp2);
p4=-ones(samp2,samp2);

% Sweep through a2 and a3 on the simplex 100 >= a2 >= a3 >=0 with two exceptions:
% the vertices a2=a3=100, and a2=a3=0. The critical points are not defined there,
% and the rcrit's are trivial to determine anyhow (always 1 and 0 respectively).

for j=1:samp2
	for k=1:min(j+1,samp2)
    
		% Fill out values of Bj		
		B1=2*sqrt(a2(j)*a3(k))*B;
		B2=2*sqrt(a1*a3(k))*B;
		B3=2*sqrt(a1*a2(j))*B;

		% Calculate p1 thru p3
		p1(j,k)=rdot(a1,a2(j),a3(k),0,B2,0,samp1,0);
		p2(j,k)=rdot(a1,a2(j),a3(k),B1,0,0,samp1,0);
		p3(j,k)=rdot(a1,a2(j),a3(k),0,0,B3,samp1,0);

		% Get first estimate for p4
		maxcands = zeros(samp3+1,1);		
		cph = (0:samp3)/samp3;
		for l=0:samp3
			cph=l/samp3;
			maxcands(l+1)=rdot(a1,a2(j),a3(k),B1*sqrt(1-cph^2),0,B3*cph,samp1,0);
		end
		[max1,posmax1]=max(maxcands);

		% Refine the estimate for p4
		if posmax1==1
			maxcands = zeros(samp3+1,1);		
			cph = (0:samp3)/samp3/samp3;	
			for l=0:samp3
				maxcands(l+1)=rdot(a1,a2(j),a3(k),...
					B1*sqrt(1-cph(l+1)^2),0,B3*cph(l+1),samp1,0);
			end
			[max2,posmax2]=max(maxcands);
			p4(j,k) = cph(posmax2);
		else if posmax1==samp3+1
			maxcands = zeros(samp3+1,1);		
			cph = 1-1/samp3+(0:samp3)/samp3/samp3;	
			for l=0:samp3
				maxcands(l+1)=rdot(a1,a2(j),a3(k),...
					B1*sqrt(1-cph(l+1)^2),0,B3*cph(l+1),samp1,0);
			end
			[max2,posmax2]=max(maxcands);
			p4(j,k) = cph(posmax2);
		else
			maxcands = zeros(2*samp3+1,1);		
			cph = (posmax1-1+(-samp3:samp3)/samp3)/samp3;	
			for l=0:2*samp3
				maxcands(l+1)=rdot(a1,a2(j),a3(k),...
					B1*sqrt(1-cph(l+1)^2),0,B3*cph(l+1),samp1,0);
			end
			[max2,posmax2]=max(maxcands);
			p4(j,k) = cph(posmax2);
		end		
	end
end			
end
