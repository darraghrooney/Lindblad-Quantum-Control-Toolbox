% This function returns a collection of "standard" Lindblad ops, that is a complete
% set of jump operators (with amplitudue uniform on [0,10], as well as one
% dephasing operator, with diagonal entries uniform on [-10,10].

function[Ls] = rand_stdLs()

Ls = zeros(3,3,13);

Ls(1,2,1) = rand()*10;
Ls(2,1,2) = rand()*10;
Ls(2,3,3) = rand()*10;
Ls(3,2,4) = rand()*10;
Ls(3,1,5) = rand()*10;
Ls(1,3,6) = rand()*10;
Ls(1,4,7) = rand()*10;
Ls(4,1,8) = rand()*10;
Ls(2,4,9) = rand()*10;
Ls(4,2,10) = rand()*10;
Ls(3,4,11) = rand()*10;
Ls(4,3,12) = rand()*10;
Ls(1,1,13) = rand()*10;
Ls(2,2,13) = rand()*20-10;
Ls(3,3,13) = rand()*20-10;
Ls(4,4,13) = rand()*20-10;

end
