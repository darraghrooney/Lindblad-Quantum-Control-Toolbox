% This function creates randomized standard Lindblad operators for 
% a n=6 system, where the first three nodes are sinks. If the dephase_switch
% switch is one, a single de-phasing operator is included. The sinks
% are then sorted so that an initial state completely mixed in the non-sink 
% subspace decays to the primary sector.

function[L] = create_3x3(dephase_switch)
  
  % Create jump operators
  L = zeros(6,6,15);
  L(1,4,1) = rand()*20;
  L(2,4,2) = rand()*20;
  L(3,4,3) = rand()*20;
  L(1,5,4) = rand()*20;
  L(2,5,5) = rand()*20;
  L(3,5,6) = rand()*20;
  L(1,6,7) = rand()*20;
  L(2,6,8) = rand()*20;
  L(3,6,9) = rand()*20;
  L(4,5,10) = rand()*20;
  L(5,4,11) = rand()*20;
  L(4,6,12) = rand()*20;  
  L(6,4,13) = rand()*20;  
  L(5,6,14) = rand()*20;  
  L(6,5,15) = rand()*20;
  
  % Create de-phasing op
  if dephase_switch
    L(:,:,16) = diag(rand(6,1)*20);  
  end
  
  % Calculate w-rates and determine where the asymptotic state is
  w = w_get(L,eye(6));
  constraints = multisink(3,3,w);
  stat = sum(constraints(4:6,:))/3/constraints(1,1);

  % Re-order sinks
  [s,i] = sort(stat, "descend");
  perm_mat = eye(6);
  perm_mat(1:3,1:3) = eye(3)(:,i);
  for l = 1:15
    L(:,:,l) = perm_mat'*L(:,:,l)*perm_mat;
  end

end

