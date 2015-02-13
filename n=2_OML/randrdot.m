% This functions picks random initial data then runs rdot.

function [rcrit, rdot, nnM, nnm,a,b] = randrdot(samples, shouldplot)

[a,b,A] = randA();
[rcrit, rdot, nnM, nnm] = rdot(a(1),a(2),a(3),b(1),b(2),b(3),samples,shouldplot);

end
