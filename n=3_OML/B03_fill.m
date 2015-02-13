function[] = B03_fill(w,grid_size1,grid_size2,grid_size3)

boundaryE = B03_boundary(w, grid_size1, grid_size2, grid_size3, 0);

patch(boundaryE(1,:),boundaryE(2,:),'g','edgecolor','g');

end
