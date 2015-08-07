# This function conducts a survey on how frequently a system has alternate threads.
# There can be zero, one or two for any system. The first return vector is 3x1
# and contains the relative frequencies, summing to N, the number of randomized
# systems. The second return value is an ?x3 matrix containing a list of the
# alternate thread sources.

function [sourceno, sources] = source_survey (N)

sourceno = zeros(4,1);
sources = [];

# Sweep over number of experiments
for j = 1:N

  # Get randomized data
  [a,b] = randA;

  # Find sources
  source_j = find_sources(a(1),a(2),a(3),b(1),b(2),b(3));

  # Add to return values
  sn = size(source_j,1);
  sources = [sources; source_j];
  sourceno(sn+1)++;
  
endfor

endfunction
