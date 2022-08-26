function [objectsOut] = mergeRegions(objectsIn)


plane_point = objectsIn{2,2};
plane_span1 = objectsIn{2,3};
plane_span2 = objectsIn{2,4};

vert_1_obj1 = plane_point;
vert_2_obj1 = plane_point+plane_span1;
vert_3_obj1 = plane_point+plane_span1+plane_span2;
vert_4_obj1 = plane_point+plane_span2;
objVert1 = [vert_1_obj1; vert_2_obj1; vert_3_obj1; vert_4_obj1];
h1 = max(objVert1(:,3));

% keep only vertices on the ground
ifilt = find(objVert1(:,3)>1e-10);
objVertGround1 = objVert1(ifilt,1:2);
theta1 = atan2(objVertGround1(:,2),objVertGround1(:,1));
objVertGround1 = [theta1, objVertGround1];
objVertGround1 = sortrows(objVertGround1,'ascend');


plane_point = objectsIn{3,2};
plane_span1 = objectsIn{3,3};
plane_span2 = objectsIn{3,4};

vert_1_obj1 = plane_point;
vert_2_obj1 = plane_point+plane_span1;
vert_3_obj1 = plane_point+plane_span1+plane_span2;
vert_4_obj1 = plane_point+plane_span2;
objVert2 = [vert_1_obj1; vert_2_obj1; vert_3_obj1; vert_4_obj1];
h2 = max(objVert2(:,3));

ifilt = find(objVert2(:,3)>1e-10);
objVertGround2 = objVert2(ifilt,1:2);
theta2 = atan2(objVertGround2(:,2),objVertGround2(:,1));
objVertGround2 = [theta2, objVertGround2];
objVertGround2 = sortrows(objVertGround2,'ascend');


%%
if min(theta2)<min(theta1)
    % swap them so we can use the logic below
    h1temp = h1;
    h2temp = h2;
    objVertGround1temp = objVertGround1;
    objVertGround2temp = objVertGround2;
    objectsInTemp = objectsIn;
    theta1temp = theta1;
    theta2temp = theta2;
    
    h1 = h2temp;
    h2 = h1temp;
    objVertGround1 = objVertGround2temp;
    objVertGround2 = objVertGround1temp;
    objectsIn(2,:) = objectsInTemp(3,:);
    objectsIn(3,:) = objectsInTemp(2,:);
    theta1 = theta2temp;
    theta2 = theta1temp;
end
%%
if isempty(theta2) || isempty(theta1)
    stopAndThink = 1;
end
if (min(theta2)>=max(theta1)) || max(theta2)<=min(theta1) % no overlap
    
    objectsOut = objectsIn;
else
    objectsOut(1,:) = objectsIn(1,:);

%     if min(theta2)<max(theta1) && max(theta2)>max(theta1) % recently
%     removed - check for issues

        
%         objectsOut{2,1} = 3;
        if h1>h2
            % keep first facet the same
            objectsOut(2,:) = objectsIn(2,:);
            
            % remove overlap from second
            objectsOut{3,1} = 3;
            objectsOut{3,2} = [objVertGround1(2,2), objVertGround1(2,3), 0]; % end of first facet starts the second
            objectsOut{3,3} = [objVertGround2(2,2), objVertGround2(2,3), 0]-objectsOut{3,2}; % span is the end of second - end of first
            objectsOut{3,4} = [0,0,h2];
            objectsOut{3,5} = objectsIn{3,5};
        else
            % truncate the first facet
            objectsOut{2,1} = 3;
            objectsOut{2,2} = [objVertGround1(1,2),objVertGround1(1,3),0];
            objectsOut{2,3} = [objVertGround2(1,2),objVertGround2(1,3),0] - objectsOut{2,2}; % start of second facet - start of first
            objectsOut{2,4} = [0,0,h1];
            objectsOut{2,5} = objectsIn{2,5};
            
            % keep the second facet the same
            objectsOut(3,:) = objectsIn(3,:);

        end

end






% figure
% scatter(objVertGround2(:,2),objVertGround2(:,3),'r')
% hold on
% scatter(objVertGround1(:,2),objVertGround1(:,3),'b')
% compute distance from corner






end

