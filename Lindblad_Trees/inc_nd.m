% This function is for incrementing on the vector space Z_(d_1)^(d_2)
% The leftmost element has the highest value (the highest digit). If 
% all digits are equal to (d1-1), the function returns all zeros.
 
function[inc_data] = inc_nd(d1,d2,data)

	% Start at rightmost digit
	place = d2;
	while (place > 0)

		% If digit is less than d1-1, increment and break			
		if (data(place) < (d1-1))
			data(place) += 1;
			break;
		% Else reset to zero and go to next digit		
		else
			data(place) = 0;
			place -= 1;
		end
	end
	inc_data = data;
end
