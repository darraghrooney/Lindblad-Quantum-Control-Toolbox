% This function plots the filled SLC region for an n=3 Lindblad system, given a
% discrete set of flags. The identity flag and its five other permutations are done
% by default. A further un_no random flags are generated with their permutations.
% Both the unfolded and folded versions can be plotted.

function [] = AE3_fill(L, grid_size, un_no, fold_toggle)

% Get curves

curves = AE3_arcs(L, grid_size, un_no, 0);
pair_count = 6*(un_no+1);

% Plot the patches

hold on
for j = 1:(pair_count-2)
	for k = j+1:(pair_count-1)
		for l = k+1:pair_count
			bnd = [curves(:,:,pairs(j,k, un_no))',curves(2:grid_size+1,:,pairs(k,l,un_no))',...
				curves(grid_size:-1:1,:,pairs(j,l,un_no))'];	
			patch(bnd(1,:),bnd(2,:),'b','edgecolor','none');
		end
	end	
end

% Adapt plot window

if (fold_toggle == 0)
	extent0 = max(max(curves(:,1,:)));
	extent1 = max(max(curves(:,2,:)));
	extent2 = min(min(curves(:,2,:)));
	axis([-1.1*extent0,1.1*extent0,1.1*extent2,1.1*extent1]);
else
	patch([0,1,0,0,-1,0],[-2/sqrt(3),1/sqrt(3),0,1/sqrt(3),1/sqrt(3),-2/sqrt(3)],'w','edgecolor','w');
	axis([-.1,1.1,-.1,1.1/sqrt(3)]);
	plot([0,1,0,0],[0,1,1,0]/sqrt(3),'k')
end

hold off
end

% Helper function that retrieves count number for a vertex pair.

function [count] = pairs(j,k, un_no)
	count = -j*(j+1)/2 + 6*(un_no+1)*(j-1) + k;
end
