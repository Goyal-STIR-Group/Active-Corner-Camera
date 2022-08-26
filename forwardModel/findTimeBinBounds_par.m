function [dSmallest, dLargest] = findTimeBinBounds_par(params,vertices,m,l,w, plotFlag)
epsilon = params.epsilon;
if plotFlag == 1
   
    figure(77)
    hold on
    scatter3(l(1),l(2),0,'r','filled') % plot laser spot
    scatter3(w(1),w(2),0,'b','filled') % plot detector spot
    scatter3(0,0,0,'filled') % plot origin
    grid on
    xlabel('x')
    ylabel('y')
    zlabel('z')
    zlim([-.5 4])
    % plot the facet
    patch(vertices(:,1),vertices(:,2),vertices(:,3),'red', 'FaceAlpha', .1)
    
end


% compute the closest time
% the closest parts of the plane correspond to the shortest times
minZ = min(vertices(:,3));
ifilt = find(abs(vertices(:,3)- minZ)< epsilon); % find lowest (in z) vertices

% line segment between closest vertices
x1 = vertices(ifilt(1),1);
x2 = vertices(ifilt(2),1);
y1 = vertices(ifilt(1),2);
y2 = vertices(ifilt(2),2);

% find slope and intersept
% come back and deal with the case where x1 = x2
if abs(x1 - x2)<epsilon
    if sign(y1) ~= sign(y2)  % line segment crosses y = 0
        dSmallest = norm([x1, 0, minZ]'-w) + norm([x1, 0, minZ]-l);
    else % then, closest point is one of the two bottom vertices
        if abs(y1)<abs(y2)
            dSmallest = norm(vertices(ifilt(1),:)-w) + norm(vertices(ifilt(1),:)-l);
        else
            dSmallest = norm(vertices(ifilt(2),:)-w) + norm(vertices(ifilt(2),:)-l);
        end
        
    end
    disp('Make sure that the fix for findTimeBinBounds.m to accomodate special case where x1 = x2 works. See line 17. ')
else
    
    mLine = (y2-y1)/(x2-x1);
    cLine = -x1*(y2-y1)/(x2-x1)+y1;
    if plotFlag == 1
        plot([x1,x2],[y1,y2])
    end
    
    dSmallest = 2*sqrt((cLine^2+(mLine*m/2)^2)/(mLine^2 + 1));
    
    
    d = dSmallest; % this is for 2D ellipse
    a = 0.5*sqrt(d^2 - m^2);
    b = d/2;
    
    if plotFlag == 1
        figure(77)
        % compute in cart coords
        nAngle = 90;
        theta = linspace(0,2*pi,nAngle);
        rpol = a*b./sqrt((b*cos(theta)).^2 +(a*sin(theta)).^2);
        ellipse_x = rpol.*cos(theta);
        ellipse_y = rpol.*sin(theta);
        plot3(ellipse_x,ellipse_y,min(vertices(:,3))*ones(size(ellipse_y)),'LineWidth',4)
    end
    
    % find intersection point of ellipse and line and plot it
    x_intersection = (-(a^2)*mLine*cLine) / (b^2+a^2*mLine^2);
    y_intersection = mLine*x_intersection+cLine;
    if plotFlag == 1
        scatter3(x_intersection,y_intersection,min(vertices(:,3)),'filled')
    end
    
    % if that point isnt on the facet segment, pick one that is!
    if x_intersection<min(x1,x2)
        x_intersection = min(x1,x2);
        y_intersection = mLine*x_intersection+cLine;
    elseif x_intersection>max(x1,x2)
        x_intersection = max(x1,x2);
        y_intersection = mLine*x_intersection+cLine;
    end
    intersection_point = [x_intersection,y_intersection,min(vertices(:,3))];
    
    if plotFlag == 1
        scatter3(intersection_point(1),intersection_point(2),intersection_point(3),'filled')
    end
    
    dSmallest = norm(intersection_point-w)+norm(intersection_point-l);
end

if plotFlag == 1
    figure(77)
    d = dSmallest;
    a = 0.5*sqrt(d^2 - m^2);
    b = d/2;
    
    % plot inner ellipsoid just to confirm that nothing is broken
    [X,Y,Z] = ellipsoid(0,0,0,a,b,a);
    surf(X,Y,Z, 'FaceAlpha', .1)
end


% now find the max distance on the plane
ifilt = find(abs(vertices(:,3)-max(vertices(:,3)))<epsilon);
dMaxCand1 = norm(vertices(ifilt(1),:)-w)+norm(vertices(ifilt(1),:)-l);
dMaxCand2 = norm(vertices(ifilt(2),:)-w)+norm(vertices(ifilt(2),:)-l);
if dMaxCand1 >= dMaxCand2
    dLargest = dMaxCand1;
else
    dLargest = dMaxCand2;
end
d = dLargest;
a = 0.5*sqrt(d^2 - m^2);
b = d/2;


if plotFlag == 1
    figure(77)
    [X,Y,Z] = ellipsoid(0,0,0,a,b,a);
    % plot outer ellipsoid
    surf(X,Y,Z, 'FaceAlpha', .1)
end

end

