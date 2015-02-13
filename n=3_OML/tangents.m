% For a given collection of Lindblad operators, and a point on orbit space,
% this function plots the discretized range of available tangent vectors.
% The number of tangent vectors goes as gno^6, where gno is the grid number,
% so only small gno's should be used. gno = 6 generally gives a good plot, and
% takes 2-8 minutes to run on a cheap laptop.

function[tngs] = tangents(Ls, x1, x2, gno)

% Initialize flag
flag = eye(3);
aux = zeros(3,2); aux(2,1) = 1; aux(3,2) = 1;
indices = zeros(6,1);

% Find number of tangent vectors
pno = (gno^2 - gno + 2)^3;

% Initialize
tngs = zeros(pno,2);

% If number of L ops is exorbitant, reduce to 8.
if (size(Ls,3) > 8)
	Ls = GKStoL(LtoGKS(Ls));
end

% Sweep over flags
for k = 1:pno

	% Get the Lindblad rates for this flag
	w = Ltow(Ls,flag);

	% Get A and b for these Lindblad rates
	[A,b,rinf] = A_b_rinf(w);

	% Calculate tangent vector
	tngs(k,:) = b - A*[x1;x2];

	% Increment flag
	[flag,aux,indices] = flag_inc(flag,aux,indices,gno);
end

% Plot tangent vectors
figure
hold on
plot(tngs(:,1),tngs(:,2),".+")

% Plot zero point for reference
plot(0,0,".@0","markersize",15)
hold off

end
