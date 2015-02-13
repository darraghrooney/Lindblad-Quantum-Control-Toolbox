% This function, for a given a1, a2 and a3, finds the rcrit for b-vectors
% by sweeping through the angles. For all angles, the norm of b in the relevant
% units is assumed to be constant. rcrit is plotted as a function of the angles.

% samp1 is the mesh size for r on [0,1] and samp3 is the mesh size for b on S2

% If shouldplot is 1, the rcrits are plotted. If false, the plot is suppressed

function[rcrit,cosphis,thetas] = bsweep(a1,a2,a3,B,samp1,samp3,shouldplot)

rcrit=zeros(samp3+1,samp3+1);

% Define grid and the angles on the grid
me = [0:samp3]/samp3;
cosphis=1-2*me;
sinphis=sqrt(1-cosphis.*cosphis);
thetas=2*pi*me;

% Determine b on the grid
b1 = 2*sqrt(a2*a3)*sinphis'*cos(thetas)*B; 
b2 = 2*sqrt(a1*a3)*sinphis'*sin(thetas)*B; 
b3 = 2*sqrt(a1*a2)*cosphis'*ones(1,samp3+1)*B; 

% Sweep over grid and calculate rcrit
for j=0:samp3
	for k=0:samp3
		rcrit(j+1,k+1)=rdot(a1,a2,a3,b1(j+1,k+1),b2(j+1,k+1),b3(j+1,k+1),samp1,0);
	end
end

% Plot if requested
if (shouldplot==1)
	mesh(thetas,cosphis,rcrit);
	xlabel("theta"); ylabel("cosphi"); zlabel("rcrit");
end

end
