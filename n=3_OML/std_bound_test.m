%% This function is intended for Lindblad systems with jump and de-phasing ops
%% only. It plots the stat states for one-parameter rotations, in the unfolded
%% or folded context. Additionally it plots the stat states for a specified 
%% number of randomized flags (as well as their permutations in the unfolded
%% context).

function[] = std_bound_test(L, uno, fold_switch)

% Curves are discretized
grid1 = 100;
bound = zeros(2,grid1+1,9);

% Unfolded context
if !fold_switch

% Sweep over possible rotations
for j1 = 1:3
for j2 = 1:3

% Sweep over discretization
for gr = 0:grid1

  % Calculate flag
	phi = gr/grid1*pi;
	flag = zeros(3);
	flag(j1,j2) = 1;
	j3a = mod(j1,3)+1;  
	j3b = mod(j1+1,3)+1;
	j4a = mod(j2,3)+1;
	j4b = mod(j2+1,3)+1;
	flag(j3a,j4a) = cos(phi/2);
	flag(j3a,j4b) = -sin(phi/2);
	flag(j3b,j4a) = sin(phi/2);
	flag(j3b,j4b) = cos(phi/2);
	
  % Calculate stat state
	w = Ltow(L,flag);
	[A,b,rinf] = A_b_rinf(w);
	bound(:,gr+1,3*(j2-1)+j1) = rinf;

end
end
end

% Sweep over flag randomizations
for j = 1:uno

  % Get random flag
	flag = flag_gen(0,1);
	
  % Sweep over permutations
  for k = 1:6
		flagper = flag(:,perms([1,2,3])(k,:));

    % Calculate stat state and plot
		[A,b,rinf] = A_b_rinf(Ltow(L,flagper));
		plot(rinf(1),rinf(2),'g');
    hold on;
	end
end

% Plot the rotation curves
for j = 1:9
	plot(bound(1,:,j),bound(2,:,j),'b','linewidth',2);
  hold on;
end

% Plot the simplex boundary
arenas;
hold off;


% Folded context
else

% Ensure the stat state will be in the primary sector
[w,J,I] = order_rates(Ltow(L,eye(3)));
perm = eye(3)(:,I);
for j = 1:size(L,3)
  L(:,:,j) = perm'*L(:,:,j)*perm;
end

count = 1;

% Sweep over rotation combinations
for j1 = 1:2
for j2 = (j1+1):3
for j3 = 1:2
for j4 = (j3+1):3

% Sweep over curve discretization
for gr = 0:grid1

  % Fetch flag
	phi = gr/grid1*pi;
	flag = zeros(3);
	flag(j1,j3) = cos(phi/2);
	flag(j1,j4) = -sin(phi/2);
	flag(j2,j3) = sin(phi/2);
	flag(j2,j4) = cos(phi/2);
	lookup = [[0,3,2];[3,0,1];[2,1,0]];
  flag(lookup(j1,j2),lookup(j3,j4)) = 1;
  
  % Calculate stat state
	w = Ltow(L,flag);
	[A,b,rinf] = A_b_rinf(w);
	bound(:,gr+1,count) = rinf;

end 

count++;
end
end
end
end

% Sweep over randomizations
for j = 1:uno

  % Generate random flag
	flag = flag_gen(0,1);

  % Calculate and plot stat state
  w = order_rates(Ltow(L,flag));
	[A,b,rinf] = A_b_rinf(w);
	plot(rinf(1),rinf(2),'g');
  hold on;
end

% Plot rotation curves
for j = 1:9

  % Remove parts of curve not in primary sector
  pcand = restrict(bound(:,:,j));
  
  % Plot if points remaining
  if size(pcand,2)>0
  	plot(pcand(1,:),pcand(2,:),'b','linewidth',2);
  end
end
		
% Plot simplex boundary
plot([0,1,0,0]/sqrt(2),[0,1,1,0]/sqrt(6),'k');
hold off;

end

end

% Helper function that removes piece of curve outside of primary sector
function[oc] = restrict(inc)

  oc = [];
  for j = 1:size(inc,2)
    if (inc(1,j) >= -1e-3 && inc(2,j) >= inc(1,j)/sqrt(3) - 1e-3)
       oc = [oc,inc(:,j)];
    end 
  end

end