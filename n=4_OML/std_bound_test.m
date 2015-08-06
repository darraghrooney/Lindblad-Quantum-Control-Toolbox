% This function takes a collection of "standard" Lindblad operators, and does 
% two things. First it plots the asymptotic states corresponding to one-parameter
% rotations from the main critical state in the primary sector. There are six of
% these rotations. Next it calculates random flags (uno is the input parameter)
% and plots the asymptotic states corresponding to them (permuted so that they
% appear in the primary sector). These random asymptotic states should be "inside"
% the rotational states.

function[] = std_bound_test(L, uno)

% Discretization of the rotation curves.
grid1 = 30;

% Permute labels so that the main critical state is in the primary sector.
bound = zeros(3,grid1+1,6);
[w,J,I] = order_rates(Ltow(L,eye(4)));
perm = eye(4)(:,I);
for j = 1:size(L,3)
  L(:,:,j) = perm'*L(:,:,j)*perm;
end

count = 1;
% Sweep over possible rotations
for j1 = 1:4
for j2 = (j1+1):4
for j3 = 1:4
for j4 = (j3+1):4
for j5 = 0:1

% Sweep over curve discretization
for gr = 0:grid1
  
  % Calculate flag
  phi = gr/grid1*pi/2;
	flag = zeros(4);
	flag(j1,j3) = cos(phi/2);
	flag(j1,j4) = -sin(phi/2);
	flag(j2,j3) = sin(phi/2);
	flag(j2,j4) = cos(phi/2);
	lookup = [[0,3,2,2];[3,0,1,1];[2,1,0,1];[2,1,1,0]];
  lookup2 = [[0,4,4,3];[4,0,4,3];[4,4,0,2];[3,3,2,0]];
  flag(lookup(j1,j2),lookup(j3,j4)) = 1-j5;
  flag(lookup(j1,j2),lookup2(j3,j4)) = j5;
  flag(lookup2(j1,j2),lookup(j3,j4)) = j5;
  flag(lookup2(j1,j2),lookup2(j3,j4)) = 1-j5;

  % Calculate asymptotic state and insert into rotation curve
	w = Ltow(L,flag);
	[A,b,rinf] = A_b_rinf(w);
	lookup = [[0,1,2,3];[0,0,4,5];[0,0,0,6]];		
  bound(:,gr+1,count)= rinf;

	
end

count++;

end
end
end
end
end

% Plot the arena
plot3([0,1,0,0,0,1,0,0]/sqrt(2),[0,1,1,0,0,1,1,0]/sqrt(6),[1,1,1,1,0,1,1,0]/sqrt(12),'k');
hold on

% Plot the rotation curves
for j = 1:72

  % Remove parts of curves not in the primary sector
  pcand = restrict(bound(:,:,j));
  if (size(pcand,2)>0)
    plot3(pcand(1,:),pcand(2,:),pcand(3,:),'b',"linewidth",2);
  end
  hold on;
end

% Sweep over random flags
for j = 1:uno

  % Find flag
	flag = flag_gen(0,1);
	w = Ltow(L,flag); 
	w = order_rates(w);
	
  % Calculate asymptotic state and plot
  [A,b,rinf] = A_b_rinf(w);
	plot3(rinf(1),rinf(2),rinf(3),'g','markersize',4,'marker','.');
end

box off;

end

% Helper function to remove curve parts not in the primary sector
function[outc] = restrict(inc)
  outc = [];
  for j = 1:size(inc,2)
    if (inc(1,j) >= -1e-3 && sqrt(3)*inc(2,j) >= inc(1,j) -1e-3 && ...
                sqrt(2)*inc(3,j) >= inc(2,j) - 1e-3)
      outc = [outc, inc(:,j)];
    end
  end
end