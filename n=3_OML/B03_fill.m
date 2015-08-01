%% This function uses the B03 arcs and fills them in to plot the B03 region.

function[] = B03_fill(w,grid_size1,grid_size2,grid_size3)

% Find boundary
boundaryE = B03_boundary(w, grid_size1, grid_size2, grid_size3, 0);

% Fill in
patch(boundaryE(1,:),boundaryE(2,:),'g','edgecolor','g');

end
