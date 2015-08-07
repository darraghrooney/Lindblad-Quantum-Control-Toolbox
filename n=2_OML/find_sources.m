# This function searches for alternate threads that pop up. Either zero, one
# or two points in the Bloch ball will serve as a new thread beginning. The
# function returns either an empty matrix or a 1x3 or 2x3 matrix depending.

function [sources] = find_sources(a1,a2,a3,b1,b2,b3)

# Calculate polynomial coefficients

c(7) = b1^2*a2^3*a3^3 + b2^2*a3^3*a1^3 + b3^2*a1^3*a2^3;
c(6) = - 3*b1^2*a2^2*a3^2*(a2+a3)- 3*b2^2*a1^2*a3^2*(a1+a3)- 3*b3^2*a2^2*a1^2*(a2+a1);
c(5) = 3*b1^2*(a2*a3^3+3*a2^2*a3^2+a2^3*a3)+3*b2^2*(a1*a3^3 ...
      +3*a1^2*a3^2+a1^3*a3)+3*b3^2*(a1*a2^3+3*a1^2*a2^2+a1^3*a2);
c(4) = -b1^2*(a3^3+a2^3+9*a2*a3*(a2+a3))-b2^2*(a3^3+a1^3+9*a1*a3*(a1+a3)) ...
      -b3^2*(a1^3+a2^3+9*a1*a2*(a1+a2));
c(3) = b1^2*(3*a3^2+9*a2*a3+3*a2^2)+b2^2*(3*a3^2+9*a1*a3+3*a1^2) ...
      +b3^2*(3*a1^2+9*a2*a1+3*a2^2);
c(2) = -3*b1^2*(a2+a3)-3*b2^2*(a1+a3)-3*b3^2*(a2+a1);
c(1) = b1^2 + b2^2 + b3^2;

# Calculate roots of polynomials
mu = roots(c);

# Use roots to find sources
sources = [];
for j = 1:size(mu)

  # Ignore complex roots
  if (abs(imag(mu(j))) < 1e-8)
    muj = real(mu(j));
    n1 = -b1/2/(a1-muj);
    n2 = -b2/2/(a2-muj);
    n3 = -b3/2/(a3-muj);

    # Ignore sources outside the Bloch ball
    if (norm([n1,n2,n3]) <= 1)
      sources = [sources; [n1,n2,n3]];
    endif

   endif
endfor

endfunction
