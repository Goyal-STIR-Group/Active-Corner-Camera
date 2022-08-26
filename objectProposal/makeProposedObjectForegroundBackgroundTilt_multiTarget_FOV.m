function [proposedObjects,rejectFlag] = makeProposedObjectForegroundBackgroundTilt_multiTarget_FOV(...
    phiBounds, foregroundHeights,foregroundRange,backgroundRange,backgroundTilt)
global wall_v0
global epsilon
global roomHeightMax
global l_xyz0
global fovCenter
% alpha - tilt of background facet

nTarget = size(phiBounds,1);
rejectFlag = 0;

objCount = 1;
for p = 1:nTarget
    facet_phi1 = phiBounds(p,1);
    facet_phi2 = phiBounds(p,2);
    objectHeight = foregroundHeights(p);
    proposedRange = foregroundRange(p);
    proposedRangeBackground = backgroundRange(p);
    alpha = backgroundTilt(p);
    
    
  
    % handle foreground object
    deltaPhi = max([facet_phi1, facet_phi2])-min([facet_phi1, facet_phi2]); % angular extent of foreground object
    % halfWidth = tan(deltaPhi/2)*proposedRange;
    meanPhi = (facet_phi1 + facet_phi2)/2;
    
    edgeRange = proposedRange/cos(deltaPhi/2); % distance to edges of facet
    
    plane_point = [edgeRange*cos(facet_phi1), edgeRange*sin(facet_phi1), 0];
    plane_point2 = [edgeRange*cos(facet_phi2), edgeRange*sin(facet_phi2), 0];
    plane_span1 = plane_point2-plane_point; % points from plane_point to plain_point2;
    
    % plane_point = [-halfWidth, proposedRange, 0];
    % plane_span1 = [2*halfWidth,0,0];
    plane_span2 = [0,0,objectHeight];
    plane_normal = -[cos(meanPhi), sin(meanPhi), 0];
    
    proposedObjects{objCount,1} = 2;
    proposedObjects{objCount,2} = plane_point;
    proposedObjects{objCount,3} = plane_span1;
    proposedObjects{objCount,4} = plane_span2;
    proposedObjects{objCount,5} = plane_normal;
    objCount = objCount + 1;
    
    maxR = max(norm(plane_point-l_xyz0),norm(plane_point2-l_xyz0));
    dR = .1;
    %%%%%
    % handle background object
    
    
    
    
    % and a point on the background plane
    x0 = proposedRangeBackground*cos(meanPhi)- l_xyz0(1);
    y0 = proposedRangeBackground*sin(meanPhi)- l_xyz0(2);
    z0 = 0- l_xyz0(3);
    
    % the normal to the background, include tilt angle
    normal_b = -[cos(meanPhi+alpha), sin(meanPhi+alpha), 0];
    
    % normal_b = plane_normal;
    A = normal_b(1);
    B = normal_b(2);
    C = normal_b(3);
    
    % foreground object vertices
    vert_f(1,:) = plane_point - wall_v0 - l_xyz0;
    vert_f(2,:) = vert_f(1,:)+plane_span1;
    vert_f(3,:) = vert_f(1,:)+plane_span1+plane_span2;
    vert_f(4,:) = vert_f(1,:)+plane_span2;
    
    % find foreground vertices projected onto the background
    
    numerator = A*x0 + B*y0 + C*z0;
    backgroundVertices = zeros(4,3);
    
    x = vert_f(:,1);
    y = vert_f(:,2);
    z = vert_f(:,3);
    
    theta_f = acos(z./sqrt(x.^2+y.^2+z.^2)); % inclination, measured down from z
    phi_f = atan2(y,x); % azimuth
    
    r = numerator./(A*cos(phi_f).*sin(theta_f) + B*sin(phi_f).*sin(theta_f) + C*cos(theta_f));
    
    x_b = r.*cos(phi_f).*sin(theta_f);
    y_b = r.*sin(phi_f).*sin(theta_f);
    z_b = r.*cos(theta_f);
    
    z_b(z_b>roomHeightMax) = roomHeightMax;
    
    
    backgroundVertices(:,1) = x_b + wall_v0(1) + l_xyz0(1);
    backgroundVertices(:,2) = y_b + wall_v0(2) + l_xyz0(2);
    backgroundVertices(:,3) = z_b + wall_v0(3) + l_xyz0(3);
    
    
    ifiltTop = find(abs(backgroundVertices(:,3))>epsilon); % indices of not-ground vertices
    meanOccludedHeight = mean(backgroundVertices(ifiltTop,3));
    backgroundVertices(ifiltTop,3) = meanOccludedHeight;
    
    
    
    
    
    plane_point_o  = backgroundVertices(1,:);
    plane_span1_o = backgroundVertices(2,:) - plane_point_o;
    plane_span2_o = backgroundVertices(4,:) - plane_point_o;
    
    proposedObjects{objCount,1} = 3; % 0, to indicate this is an occluded region
    proposedObjects{objCount,2} = plane_point_o;
    proposedObjects{objCount,3} = plane_span1_o;
    proposedObjects{objCount,4} = plane_span2_o;
    [ proposedObjects ] = compute_normal( proposedObjects, objCount );
    objCount = objCount + 1;
    
    if min(r)<=(maxR+dR)
        rejectFlag = 1;
        disp('Rejecting because of extreme geometry.')
    end
    
    %%%%%% repeat again to compute the region occluded from view of the floor
    
    
    % and a point on the background plane
    x0 = proposedRangeBackground*cos(meanPhi)-fovCenter(1);
    y0 = proposedRangeBackground*sin(meanPhi)-fovCenter(2);
    z0 = 0- fovCenter(3);
    
    % the normal to the background, include tilt angle
    normal_b = -[cos(meanPhi+alpha), sin(meanPhi+alpha), 0];
    
    % normal_b = plane_normal; % we can get more complicated later, for now the facet points towards the corner like the foreground facet
    A = normal_b(1);
    B = normal_b(2);
    C = normal_b(3);
    
    % foreground object
    vert_f(1,:) = plane_point - wall_v0 - fovCenter;
    vert_f(2,:) = vert_f(1,:)+plane_span1;
    vert_f(3,:) = vert_f(1,:)+plane_span1+plane_span2;
    vert_f(4,:) = vert_f(1,:)+plane_span2;
    
    % find foreground vertices projected onto the background
    numerator = A*x0 + B*y0 + C*z0;
    backgroundVertices = zeros(4,3);
    
    x = vert_f(:,1);
    y = vert_f(:,2);
    z = vert_f(:,3);
    
    theta_f = acos(z./sqrt(x.^2+y.^2+z.^2)); % inclination, measured down from z
    phi_f = atan2(y,x); % azimuth
    
    r = numerator./(A*cos(phi_f).*sin(theta_f) + B*sin(phi_f).*sin(theta_f) + C*cos(theta_f));
    
    x_b = r.*cos(phi_f).*sin(theta_f);
    y_b = r.*sin(phi_f).*sin(theta_f);
    z_b = r.*cos(theta_f);
    
    z_b(z_b>roomHeightMax) = roomHeightMax;
    
    
    backgroundVertices(:,1) = x_b + wall_v0(1) + fovCenter(1);
    backgroundVertices(:,2) = y_b + wall_v0(2) + fovCenter(2);
    backgroundVertices(:,3) = z_b + wall_v0(3) + fovCenter(3);
    
    if min(r)<=(maxR + dR)
        rejectFlag = 1;
        disp('Rejecting because of extreme geometry.')
    end
    
    if min(backgroundVertices(:,3))<-1e-6
        stopAndThink = 1;
    end
    if min(backgroundVertices(:,2))<0
        stopAndThink = 1;
    end
    % for now, we assume that the background occluded region is a rectangular
    % facet. In reality, it could be a trapezoid or something. Use the mean height
    % of the middle as the facet height:
    ifiltTop = find(abs(backgroundVertices(:,3))>epsilon); % indices of not-ground vertices
    meanOccludedHeight = mean(backgroundVertices(ifiltTop,3));
    backgroundVertices(ifiltTop,3) = meanOccludedHeight;
    
    
    % ifiltTop = find(abs(backgroundVertices(ii,3))>epsilon); % indices of not-ground vertices
    % meanOccludedHeight = mean(backgroundVertices(ifiltTop,3));
    % backgroundVertices(ifiltTop,3) = meanOccludedHeight;
    
    plane_point_o  = backgroundVertices(1,:);
    plane_span1_o = backgroundVertices(2,:) - plane_point_o;
    plane_span2_o = backgroundVertices(4,:) - plane_point_o;
    
    proposedObjects{objCount,1} = 3; % 0, to indicate this is an occluded region
    proposedObjects{objCount,2} = plane_point_o;
    proposedObjects{objCount,3} = plane_span1_o;
    proposedObjects{objCount,4} = plane_span2_o;
    [ proposedObjects ] = compute_normal( proposedObjects, objCount );
    objCount = objCount + 1;
end






end

