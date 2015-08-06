% This function calculate the w-rates for a collection of Lindblad ops

function[w] = w_get(L,flag)

Lno = size(L,3);
dime = size(L,1);
w=zeros(dime,dime);

for l = 1:Lno
  for j = 1:dime
  for k = setdiff(1:6, j)
    w(j,k) += flag(:,j)'*L(:,:,l)*flag(:,k)*flag(:,k)'*L(:,:,l)'*flag(:,j);
  end
  end
end

end
