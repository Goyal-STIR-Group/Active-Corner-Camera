function [angularBounds,nIntegral,nArcSeg] = findAngularBounds_verticalFacet_condition_par(small_vertices_rsu,...
    A,B,plotFlag,maxDist,params)
epsilon = params.epsilon;
lineSegs = [1 2; ...
    2 3; ...
    3 4;...
    1 4]; % combos of vertices that make plane boundaries
nseg = size(lineSegs,1);
% compute eqs for different line segments
mLine = zeros(nseg,1);
cLine = zeros(nseg,1);
intersectFlag = zeros(nseg,1);
intersections = zeros(2,2);
count = 1;

% figure(20)
% patch(small_vertices_rsu(:,1),small_vertices_rsu(:,2),'r')
% grid on
% hold on
% xlabel('x')
% ylabel('y')
% % 
% nAngle = 90;
% theta = linspace(0,2*pi,nAngle); % polar angles for plotting
% rpol = A*B./sqrt((B*cos(theta)).^2 +(A*sin(theta)).^2);
% ellipse_r = rpol.*cos(theta);
% ellipse_s = rpol.*sin(theta);
% plot(ellipse_r,ellipse_s);
for ii = 1:nseg
    % vertices of line segment
    vertex1 = lineSegs(ii,1);
    vertex2 = lineSegs(ii,2);
    
    % pick out vertex coords
    x1 = small_vertices_rsu(vertex1,1);
    y1 = small_vertices_rsu(vertex1,2);
    
    x2 = small_vertices_rsu(vertex2,1);
    y2 = small_vertices_rsu(vertex2,2);
    
%     scatter(x1,y1,'filled')
%     scatter(x2,y2,'filled')
    
    
%     line([x1,y1],[x2,y2])
    if abs(x1 - x2) < epsilon % we are on vertical edge of facet
        % check for intersection
        temp1 = x1^2/A^2;
        if temp1<=1
            ycand1 = B*sqrt(1-temp1);
            ycand2 = -ycand1;
            % check the candidates to make sure they are on the SEGMent
            if (ycand1<max(y1,y2)) && (ycand1>=min(y1,y2))
                if plotFlag == 1
                    figure(20)
                    scatter(x1,ycand1,'filled')
                end
                intersections(count,:) = [x1,ycand1];
                count = count+1;
            end
            if (ycand2<max(y1,y2)) && (ycand2>=min(y1,y2))
                if plotFlag == 1
                    figure(20)
                    scatter(x1,ycand2,'filled')
                end
                intersections(count,:) = [x1,ycand2];
                count = count+1;
            end
        end
    else % we are on horizontal edge of facet, y2 = y1
        temp2 = y2^2/B^2;
        if temp2<=1
            xcand1 = A*sqrt(1-temp2);
            xcand2 = -xcand1;
            if  (xcand1<max(x1,x2)) && (xcand1>=min(x1,x2))
                if plotFlag == 1
                    figure(20)
                    scatter(xcand1,y1,'filled')
                end
                intersections(count,:) = [xcand1,y1];
                count = count+1;
            end
            if  (xcand2<max(x1,x2)) && (xcand2>=min(x1,x2))
                if plotFlag == 1
                    figure(20)
                    scatter(xcand2,y1,'filled')
                end
                intersections(count,:) = [xcand2,y1];
                count = count+1;
            end
        end
    end
end
thetaIntersections = atan2(intersections(:,2),intersections(:,1));
thetaIntersections = mod(thetaIntersections+2*pi,2*pi); % make positive

[thetaIntersections,sortI] = sort(thetaIntersections,'ascend');
intersections = intersections(sortI,:);

% deal with the case where we have two segments (i.e. a fat facet with top
% corners in two pieces
nArcSeg = length(thetaIntersections)/2;

if (round(nArcSeg)-nArcSeg)>0
    disp('nArcSeg not an integer!')
end
angularBounds = zeros(2,nArcSeg);
nIntegral = zeros(nArcSeg,1);
for ns = 1:nArcSeg
    
    distanceSpanned = sqrt((intersections(2*(ns-1)+1,1) - intersections(2*(ns-1)+2,1))^2 ...
        + (intersections(2*(ns-1)+1,2) - intersections(2*(ns-1)+2,2))^2);
    if distanceSpanned>maxDist
        nIntegral(ns) = ceil(distanceSpanned/maxDist);
    else
        nIntegral(ns) = 1;
    end
    
    % each column of angularBounds will be bounds for a different segment
    angularBounds(1,ns) = thetaIntersections(2*(ns-1)+1);
    angularBounds(2,ns) = thetaIntersections(2*(ns-1)+2);

end

    

% % find distance between two points and decide if we need to split the
% % integral!
% distanceSpanned = sqrt((intersections(1,1) - intersections(2,1))^2 ...
%         + (intersections(1,2) - intersections(2,2))^2);
% if distanceSpanned>maxDist
%     nIntegral = ceil(distanceSpanned/maxDist);
% else
%     nIntegral = 1;
% end
% 
% thetaIntersections = atan2(intersections(:,2),intersections(:,1));
% thetaIntersections = mod(thetaIntersections+2*pi,2*pi); % make positive
% if (thetaIntersections(2)-thetaIntersections(1))>pi
%     thetaIntersectionsOld = thetaIntersections;
%     thetaIntersections(1) = -(2*pi-thetaIntersectionsOld(2));
%     thetaIntersections(2) = thetaIntersectionsOld(1);
% end
end

