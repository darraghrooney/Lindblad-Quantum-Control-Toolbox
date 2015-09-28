# This function conducts a survey on how frequently a system has alternate threads.
# There can be zero, one or two for any system. The first return vector is 3x1
# and contains the relative frequencies, summing to N, the number of randomized
# systems. The second return value is an ?x3 matrix containing a list of the
# alternate thread sources.

function [sourceno, typeO, typeI, typeII] = source_survey (Na, Nb)

sourceno = zeros(3,1);
typeO = [];
typeI = [];
typeII = [];

# Sweep over number of experiments
for k1 = 1:Na
for k2 = 1:k1

a=[100, 100*(2*k1-1)/2/Na, 100*(2*k2-1)/2/Na]';

for j = 1:Nb

  # Get randomized data
  b = randb(a);

  # Find sources
  source_j = find_sources(a(1),a(2),a(3),b(1),b(2),b(3));

  # Add to return values
  sn = size(source_j,1);
  if (sn == 0)
    typeO = [typeO; a', b'];
  elseif(sn == 1) 
    typeI = [typeI; a', b', source_j(1,:)];
  else
    if ( norm( source_j(1,:)) < norm( source_j(2,:)) )
      typeII = [typeII; a', b', source_j(1,:), source_j(2,:)];
    else
      typeII = [typeII; a', b', source_j(2,:), source_j(1,:)];
    endif
  endif
  sourceno(sn+1)++;
  
endfor
endfor
endfor
endfunction
