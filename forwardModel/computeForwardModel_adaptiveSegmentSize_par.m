function [data,data_background,data_foreground,data_occluded] = ...
    computeForwardModel_adaptiveSegmentSize_par(objects,params,plotFlag)
% global c_light
% global bin_size
% global bin_bounds_mid


num_time_bins = params.num_time_bins;
n_wall = params.n_wall;
wall_v0 = params.wall_v0;
wall_v1 = params.wall_v1;
X = params.X;
Y = params.Y;
theta_mid_vec = params.theta_mid_vec;
transparency = params.transparency;
nPix = params.nPix;
l_xyz0 = params.l_xyz0;
maxDist = params.maxDist;
epsilon = params.epsilon;
c_light = params.c_light;



% and now, for each object, find phi_min and phi_max
for oi = 1:size(objects,1) % loop through scene objects
    % find vector pointing from corner edge to bottom vertices
    
    plane_point = objects{oi,2};
    plane_span1 = objects{oi,3};
    plane_span2 = objects{oi,4};
    
    % we don't want to have to compute this inside of the loop so store here
    vert_1_obj1 = plane_point;
    vert_2_obj1 = plane_point+plane_span1;
    vert_3_obj1 = plane_point+plane_span1+plane_span2;
    vert_4_obj1 = plane_point+plane_span2;
    vertices_obj1 = [vert_1_obj1; vert_2_obj1; vert_3_obj1; vert_4_obj1];
    
    % find closest distance of facet to origin and use that to compute a
    % reasonable maxDist
    closestDist = min(sqrt(vertices_obj1(:,1).^2+vertices_obj1(:,2).^2+vertices_obj1(:,3).^2));
    if closestDist<.8
        maxDist = .1;
    else
        maxDist = .4;
    end
    
    objects{oi,6} = vertices_obj1;
%     
%     figure(88)
%     hold on
%     %plot the facet
%     if objects{oi,1} == 1 % if background
%         patch(vertices_obj1(:,1),vertices_obj1(:,2),vertices_obj1(:,3),'red', 'FaceAlpha', 1.5*transparency)
%     elseif objects{oi,1} == 2 % foreground
%         patch(vertices_obj1(:,1),vertices_obj1(:,2),vertices_obj1(:,3),'green', 'FaceAlpha', 1.5*transparency)
%     else % occluded region
%         patch(vertices_obj1(:,1),vertices_obj1(:,2),vertices_obj1(:,3),'blue', 'FaceAlpha', 1.5*transparency)
%     end
    
    objVert = objects{oi,6};
    
    % pick out just the two bottom vertices
    ifilt = find(abs( objVert(:,3) - min(objVert(:,3)))<=epsilon);
    
    % we care only about azimuthal angle, so ignore height here
    obj_v1 = [objVert(ifilt(1),1), objVert(ifilt(1),2),0];
    obj_v2 = [objVert(ifilt(2),1), objVert(ifilt(2),2),0];
    
    % find the vector pointing from corner to object vertices
    facet_corner_vec1 = -(wall_v0 - obj_v1);
    facet_corner_vec2 = -(wall_v0 - obj_v2);
    
    % find angles of facets
    facet_phi1 = acos(dot(facet_corner_vec1,n_wall)/norm(facet_corner_vec1));
    facet_phi2 = acos(dot(facet_corner_vec2,n_wall)/norm(facet_corner_vec2));
    
    phiMin = min(facet_phi1, facet_phi2);
    phiMax = max(facet_phi1, facet_phi2);
    objects{oi,7} = [phiMin, phiMax]; % save these angles
    
    % and tag/save the indices of vertices with the largest phi angle
    % doing this in advance so we don't have to repeat this computation
    if facet_phi1>facet_phi2 % angle 1 is larger, find those indices
        ifilt_largestPhi = find( ((abs(objVert(:,1) - objVert(ifilt(1),1))<epsilon) & (abs(objVert(:,2) - objVert(ifilt(1),2))< epsilon) )) ;
        v_big = obj_v1; % vertex with biggest angle into scene
        span_vec_facet = -(obj_v1 - obj_v2);  % vector pointing from v_big to the other vertex
    else
        ifilt_largestPhi = find( ((abs(objVert(:,1) - objVert(ifilt(2),1))<epsilon) & (abs(objVert(:,2) - objVert(ifilt(2),2))<epsilon))) ;
        v_big = obj_v2; % vertex with biggest angle into scene
        span_vec_facet = -(obj_v2 - obj_v1);  % vector pointing from v_big to the other vertex
    end
    
    %     indexLargestPhi = zeros(4,1);
    %     indexLargestPhi(ifilt_largestPhi) = 1;
    objects{oi,8} = ifilt_largestPhi; % save the indices of vertices with the largest phi angle
    objects{oi,9} =  v_big;
    objects{oi,10} =  span_vec_facet;
    
end


data_background = zeros(nPix,num_time_bins); 
data_foreground = zeros(nPix,num_time_bins);
data_occluded = zeros(nPix,num_time_bins);


% 
% occlusionIndicator = zeros(size(X));


parfor fp = 1:nPix% loop through camera pixels
    

    
    theta_mid = theta_mid_vec(fp);
%     theta_max = theta_max_vec(fp);
%     theta_min = theta_min_vec(fp);
    
    pixX = X(fp);
    pixY = Y(fp);
    w_xyz0 = [pixX, pixY, 0]; % in the original coordinate system (w in notes)
    
    % change coordinate system so that we are oriented with l on negative y
    % and w on positive y axis
    % the new basis vectors
    m = norm(w_xyz0-l_xyz0); % distance between l and w
    ey = (w_xyz0-l_xyz0)/m; % unit vector pointing from l_xyz0 to pix_xyz0
    ez = [0 0 1];
    ex = cross(ey,ez);
    
    R = [ex; ey; ez];
    o = l_xyz0+(m/2)*ey; % vector to the center of shifted system
    
    if plotFlag == 1
        figure(88)
        hold on
        scatter3(w_xyz0(1),w_xyz0(2),w_xyz0(3),'*') % plot the detector
        scatter3(o(1),o(2),o(3),'*') % plot the center of new coordinate system
    end

    
    

    w = xyz2rsu(R,o,w_xyz0,1);
    l = xyz2rsu(R,o,l_xyz0,1);
    
    % plot the new coordinate system
    if plotFlag == 1
        figure(99)
        hold on
        axis square
        scatter3(l(1),l(2),0,'filled') % plot laser spot
        scatter3(w(1),w(2),0,'filled') % plot detector spot
        scatter3(0,0,0,'filled') % plot origin
        grid on
        xlabel('x')
        ylabel('y')
        zlabel('z')
        zlim([-.5 4])
        xlim([-1 4])
        ylim([-1 4])
    end
    
    for oi = 1:size(objects,1) % loop through scene objects
%         if   fp == 982 && plotFlag ==982 && oi == 2
%             stophere = 1;
%         end
        
        
        phi_min = objects{oi,7}(1);
        phi_max = objects{oi,7}(2);
        
        if (phi_min + theta_mid) >= pi % completely occluded
%             occlusionIndicator(fp) = 0;
            % no nothing!
        else % determine if partly occluded, adjust the bounds of the facet
            vertices_xyz0 = objects{oi,6}; % vertices of facet
            
            if (phi_max + theta_mid) <= pi % completely visible
%                 occlusionIndicator(fp) = 1;
            else % partly visible, adjust vertices
                new_phi_max = pi - theta_mid;
                
                if new_phi_max < phi_max
                    indexLargestPhi = objects{oi,8}; % indices to change, switched to have origin at the corner
                    v_big = objects{oi,9} - wall_v0; % vertex to change
                    
                    if plotFlag == 1
                        figure(33)
                        hold off
                        scatter3(0,0,0,'*')
                        hold on
                        scatter3(l_xyz0(1),l_xyz0(2),0,'filled','r') % plot laser spot
                        plot3([wall_v0(1) wall_v1(1)],[wall_v0(2) wall_v1(2)],[0 0],'lineWidth', 3) % plot wall
                        scatter3(w_xyz0(1),w_xyz0(2),0,'filled','b') % plot pixel
                        patch(vertices_obj1(:,1),vertices_obj1(:,2),vertices_obj1(:,3),'red', 'FaceAlpha', 1.5*transparency)
                        xlim([-1 4])
                        ylim([-1 4])
                        axis square
                        scatter3(v_big(1),v_big(2),v_big(3),'g','filled') % plot initial vertex to change
                    end
                    
                    span_vec_facet = objects{oi,10}; % vector pointing from vertex to change to other vertex
                    
                    if plotFlag == 1
                        figure(33)
                        scatter3(span_vec_facet(1),span_vec_facet(2),0,'filled','o')
                        test = v_big + 0.5*span_vec_facet;
                        scatter3(test(1),test(2),test(3),'g','*') % plot a point in the middle, just to check parameterization
                    end
                    
                    tt = (sin(new_phi_max) * v_big(1) - cos(new_phi_max) * v_big(2))/...
                        (cos(new_phi_max) * span_vec_facet(2) - sin(new_phi_max) * span_vec_facet(1)); % line parameter value for when line is at angle new_phi_max
                    
                    % range at that parameter value
                    RR = (v_big(1) + tt*span_vec_facet(1))/cos(new_phi_max);
                    vertices_xyz0(indexLargestPhi(1),:) = [RR*cos(new_phi_max),RR*sin(new_phi_max),...
                        vertices_xyz0(indexLargestPhi(1),3)] + wall_v0;
                    vertices_xyz0(indexLargestPhi(2),:) = [RR*cos(new_phi_max),RR*sin(new_phi_max),...
                        vertices_xyz0(indexLargestPhi(2),3)] + wall_v0;
                    
                    % test
                    if plotFlag == 1
                        figure(33)
                        scatter3(RR*cos(new_phi_max),RR*sin(new_phi_max),0,'filled','b') % plot new vertex, make sure its on the edge
                        scatter3(vertices_xyz0(indexLargestPhi(1),1),vertices_xyz0(indexLargestPhi(1),2),...
                            vertices_xyz0(indexLargestPhi(1),3),'filled','g')
                        scatter3(vertices_xyz0(indexLargestPhi(2),1),vertices_xyz0(indexLargestPhi(2),2),...
                            vertices_xyz0(indexLargestPhi(2),3),'filled','g')
                    end
%                     occlusionIndicator(fp) = .5;
                end
            end
            vertices = xyz2rsu(R,o,vertices_xyz0,4);
            plane_v1 = vertices(1,:)-vertices(2,:); % points from 2 to 1
            plane_v2 = vertices(1,:)-vertices(3,:); % points from 3 to 1
            plane_n = cross(plane_v1, plane_v2);
            plane_n = plane_n./norm(plane_n);
            
            
            plane_v1 = vertices(1,:)-vertices(2,:); % points from 2 to 1, in new coordinates
            
            perp_to_plane_n = plane_v1./norm(plane_v1); % unit 1 vector perpendicular to plane_n
            
            % distance to plane from origin
            dd = dot(plane_n,vertices(1,:));
            
            % point in plane
            q = dd*plane_n/norm(plane_n)^2;
            
            
            if plotFlag == 1
                figure(99)
                patch(vertices(:,1),vertices(:,2),vertices(:,3),'red', 'FaceAlpha', .1)
                hold on
                scatter3(q(1),q(2),q(3),'filled')
            end
            
            
            [dSmallest, dLargest] = findTimeBinBounds_par(params,vertices,m, l,w,plotFlag);
%             if dSmallest>10
%                 % this happens sometimes when the closer parts of the facet
%                 % are occluded from this pixel
%                 stopHere = 1;
%                 figure(88)
%                 plotObjects(objects(oi,:),'red')
%                 xlabel('x')
%                 ylabel('y')
%             end
             [histogram] = makeHistogram_condition_circle_par(params,vertices,dSmallest,dLargest,q,m,l,w,...
                perp_to_plane_n, plane_n,plotFlag);

%             [histogram,snapOut] = makeHistogram_condition_fast_troubleshoot(vertices,dSmallest,dLargest,q,m,l,w,...
%                 perp_to_plane_n, plane_n,plotFlag,maxDist);
%             if snapOut == 1
%                 stopAndFix = 1;
%             end
%             histogram = makeHistogram(vertices,dSmallest,dLargest,q,m,l,w,...
%                 perp_to_plane_n, plane_n,plotFlag);
            
            if objects{oi,1} == 1% background
                data_background(fp,:) = data_background(fp,:) + histogram;
            elseif  objects{oi,1} == 2% foreground
                data_foreground(fp,:) = data_foreground(fp,:) + histogram;
            else % occluded
                data_occluded(fp,:) = data_occluded(fp,:) + histogram;
            end
        end
    end
    
%     if plotFlag == 1
%         figure(34)
%         plot(histogram);
%         close all
%     end
    
end


data = data_background + data_foreground - data_occluded;

end

