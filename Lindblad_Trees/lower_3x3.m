% This function takes a n=6 system with three sinks, and calculates the curves
% of stat-states as one makes one-parameter rotations between the sinks. 
% Additionally, it calculates the stat-states for a specified number of 
% random unitary flags (so basically Monte Carlo). The stat-states depend on the
% initial condition: we assume the initial state is a mixture of the three non-
% sinks, and the exact mixture is given as the input vector lamb.

function[stat] = lower_3x3(L, uno, lamb)

% Initialize
grid1 = 10;
count = 1;
stat = zeros(2,grid1+1,9);
Pi = [[1,-1,0]/sqrt(2);[1,1,-2]/sqrt(6)];

% Sweep over rotation combinations
for j1 = 1:2
for j2 = (j1+1):3
for j3 = 1:2
for j4 = (j3+1):3

% Sweep over curve discretization
for gr = 0:grid1

  % Fetch flag
	phi = gr/grid1*pi;
	flag = zeros(6);
  flag(4:6,4:6) = eye(3);
	flag(j1,j3) = cos(phi/2);
	flag(j1,j4) = -sin(phi/2);
	flag(j2,j3) = sin(phi/2);
	flag(j2,j4) = cos(phi/2);
	lookup = [[0,3,2];[0,0,1]];
  flag(lookup(j1,j2),lookup(j3,j4)) = 1;
  
  % Calculate rotational stat-states
	w = w_get(L,flag);
  constraints = multisink(3,3,w);
  stat(:,gr+1,count) = Pi*(lamb*constraints(4:6,:)/constraints(1,1))';

end 

count++;
end
end
end
end

% Calculate the Monte Carlo stat-states

for j = 1:uno
  flag(4:6,4:6) = eye(3);
  flag(1:3,1:3) = flag3d_gen(0,1);
	w = w_get(L,flag);
  constraints = multisink(3,3,w);
  ss = (lamb*constraints(4:6,:))'/constraints(1,1);

  % For each random flag, also use the other five permutations
  for j = 1:6
    pm = eye(3)(:,perms([1,2,3])(j,:));
    ssp = Pi*pm*sort(real(ss),'descend');
    plot(ssp(1),ssp(2),'g','markersize', 4);
    hold on
  end
end

% Plot the rotational stat-states

for j = 1:9
  pcand = stat(:,:,j);
  if (size(pcand,2)>0)
    plot(pcand(1,:),pcand(2,:),'b','linewidth', 2);
    hold on;
  end
end

% Plot arena

plot([0,0,1,-1,0]/sqrt(2),[1,-2,1,1,-2]/sqrt(6), 'k');
plot([1,-1/2]/sqrt(2),[1,-1/2]/sqrt(6), 'k');
plot([-1,1/2]/sqrt(2),[1,-1/2]/sqrt(6), 'k');
axis([-.75,.75,-1,.5]);
hold off

end
