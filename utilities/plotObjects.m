function [] = plotObjects(objects,color)

nOb = size(objects,1);
scatter(0,0,'filled')
grid on
axis square
hold on

% plot rings on floor
angle = linspace(0,pi,40);
for ring = 1:4
    range = .5*ring;
    plot(range*cos(angle),range*sin(angle),'k:','LineWidth',1.5)
end

for iBackground = 1:nOb
    plane_point = objects{iBackground,2};
    plane_span1 = objects{iBackground,3};
    plane_span2 = objects{iBackground,4};
    
    vert_1_obj1 = plane_point;
    vert_2_obj1 = plane_point+plane_span1;
    vert_3_obj1 = plane_point+plane_span1+plane_span2;
    vert_4_obj1 = plane_point+plane_span2;
    objVert = [vert_1_obj1; vert_2_obj1; vert_3_obj1; vert_4_obj1];
    
    

    patch(objVert(:,1),objVert(:,2),objVert(:,3),color,'FaceAlpha',.15)
end
plot([-10,0],[0,0],'k','LineWidth',4) % plot wall footprint

set(gca, 'XDir','reverse')

xlim([-1.5 1.5])
ylim([-.5 3])
zlim([0 2.5])

objVert = [-.25,0,0;...
    .25,0,0;...
    .25,-.5,0;...
    -.25,-.5,0];

patch(objVert(:,1),objVert(:,2),objVert(:,3),'k','FaceAlpha',.15)
view(-15,25)





end

