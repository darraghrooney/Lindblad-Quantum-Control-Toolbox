% This function is for plotting the boundary and inner edges of mathcal{T} for n = 4

function[] = arenas()


% Plot boundary
plot3([1/sqrt(2),-1/sqrt(2),0,1/sqrt(2),0,-1/sqrt(2)],[1/sqrt(6),1/sqrt(6),-sqrt(2/3),1/sqrt(6),0,1/sqrt(6)],[1/sqrt(12),1/sqrt(12),1/sqrt(12),1/sqrt(12),-sqrt(3/4),1/sqrt(12)],'k');
hold on
plot3([0,0],[0,-sqrt(2/3)],[-sqrt(3/4),1/sqrt(12)],'k')

% Plot interior edges
plot3([0,0],[0,0],[-sqrt(3/4),1/sqrt(12)],'k')
plot3([1/sqrt(2),-1/3/sqrt(2)],[1/sqrt(6),-1/3/sqrt(6)],[1/sqrt(12),-1/3/sqrt(12)],'k')
plot3([-1/sqrt(2),1/3/sqrt(2)],[1/sqrt(6),-1/3/sqrt(6)],[1/sqrt(12),-1/3/sqrt(12)],'k')
plot3([0,0],[-sqrt(2/3),sqrt(2)/3/sqrt(3)],[1/sqrt(12),-1/6/sqrt(3)],'k')

% Plot inner edges on the boundary
plot3([1/sqrt(2),-1/2/sqrt(2)],[1/sqrt(6),-1/2/sqrt(6)],[1,1]/sqrt(12),'k')
plot3([-1/sqrt(2),1/2/sqrt(2)],[1/sqrt(6),-1/2/sqrt(6)],[1,1]/sqrt(12),'k')
plot3([0,0],[1,-2]/sqrt(6),[1,1]/sqrt(12),'k')
plot3([0,0],[0,1/sqrt(6)],[-sqrt(3/4),1/sqrt(12)],'k')
plot3([0,-1/2/sqrt(2)],[0,-1/2/sqrt(6)],[-sqrt(3/4),1/sqrt(12)],'k')
plot3([0,1/2/sqrt(2)],[0,-1/2/sqrt(6)],[-sqrt(3/4),1/sqrt(12)],'k')
plot3([1/sqrt(2),0,-1/sqrt(2)],[1/sqrt(6),-1/sqrt(6),1/sqrt(6)],[1/sqrt(12),-1/sqrt(12),1/sqrt(12)],'k')
plot3([1/sqrt(2),-1/2/sqrt(2),0],[1/sqrt(6),1/2/sqrt(6),-2/sqrt(6)],[1/sqrt(12),-1/sqrt(12),1/sqrt(12)],'k')
plot3([-1/sqrt(2),1/2/sqrt(2),0],[1/sqrt(6),1/2/sqrt(6),-2/sqrt(6)],[1/sqrt(12),-1/sqrt(12),1/sqrt(12)],'k')

box off
axis([-.75,.75,-.83,.5,-.9,.3])

endfunction