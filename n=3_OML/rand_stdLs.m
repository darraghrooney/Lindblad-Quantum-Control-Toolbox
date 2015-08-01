% Generates random Lindblad operators that are "standard". ie there are six
% jump operators and one de-phasing operator.

function[Ls] = rand_stdLs()

Ls = zeros(3,3,7);

Ls(1,2,1) = rand()*10;
Ls(2,1,2) = rand()*10;
Ls(2,3,3) = rand()*10;
Ls(3,2,4) = rand()*10;
Ls(3,1,5) = rand()*10;
Ls(1,3,6) = rand()*10;
Ls(1,1,7) = rand()*10;
Ls(2,2,7) = rand()*10;
Ls(3,3,7) = rand()*10;

end
