% This function generates a random GKS matrix for n=2 then sweeps over the
% complete flags to find values of w12 and w21. The range of w12 and w12
% is then plotted. 

% Note that w_ij := the norm-squared of sum(k) wf_i L_k wf_k, where L_k are the
% Lindblad operators and wf_i are the basis vectors of the flag.

function [w12, w21, A] = wsweep( mesh )	% A mesh size around 40 I think is optimal

% Initialize matrices
w12 = zeros(mesh+1,mesh+1);
w21 = zeros(mesh+1,mesh+1);
paulis = zeros(2,2,3);
paulis(1,2,1) = 1; paulis(2,1,1) = 1;
paulis(1,2,2) = -1i; paulis(2,1,2) = 1i;
paulis(1,1,3) = 1; paulis(2,2,3) = -1;
wf1 = zeros(2,1); wf2 = zeros(2,1);

% Generate GKS matrix.
[a,b,A] = randA;

% Sweep over set of flags
for j = 0:mesh
    for k = 0:mesh
    	phi = pi*j/mesh;
        theta = 2*pi*k/mesh;
        wf1(1) = cos(phi/2); 
	wf1(2) = sin(phi/2)*exp(1i*theta);
        wf2(1) = -sin(phi/2); 
	wf2(2) = cos(phi/2)*exp(1i*theta); 
   
% Compute w12 and w21    	
	for l = 1:3
            for m = 1:3
                w12(j+1,k+1) = w12(j+1,k+1) + ...
			A(l,m)*wf1'*paulis(:,:,l)*(wf2*wf2')*paulis(:,:,m)*wf1;
                w21(j+1,k+1) = w21(j+1,k+1) + ...
			A(l,m)*wf2'*paulis(:,:,l)*(wf1*wf1')*paulis(:,:,m)*wf2;
            end
        end
    end
end

% Plot
w12 = real(w12);
w21 = real(w21);
plot(w12,w21,"1");
hold on;
plot(w12',w21',"1");	% Plot the transpose so that "gridlines" are shown in both "dimensions"
hold off;
xlabel("w12"); ylabel("w21");

end

