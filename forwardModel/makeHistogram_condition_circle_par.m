function [histogram] = makeHistogram_condition_circle_par(params,vertices,dSmallest,dLargest,q,m,l,w,...
    perp_to_plane_n, plane_n, plotFlag)

bin_bounds_mid = params.bin_bounds_mid;
num_time_bins = params.num_time_bins;
maxDist = params.maxDist;
bin_size = params.bin_size;
n_w = params.n_w;
n_l = params.n_l;
c_light = params.c_light;

bin_size_m = c_light*bin_size; % bin size in meter

% find the first time bin affected by this facet
start_bin = max(floor((dSmallest/c_light)/bin_size)+1,1);
stop_bin = min(floor((dLargest/c_light)/bin_size)+1,num_time_bins);

if stop_bin>num_time_bins
    stop_bin = num_time_bins;
end
if start_bin>num_time_bins
    start_bin=num_time_bins;
end
% check the last bin and make sure the center of it isnt past dLargest
% stop_bin
% start_bin
% dLargest
% dSmallest

if bin_bounds_mid(stop_bin)>=dLargest %  if the middle of the last bin is past the facet ignore
    stop_bin = stop_bin - 1;
end

if bin_bounds_mid(start_bin)<=dSmallest % if the middle of the first bin is before the facet 
    start_bin = start_bin + 1;

end

nAngle = 90;
theta = linspace(0,2*pi,nAngle); % polar angles for plotting




histogram = zeros(1,num_time_bins);
for bin_i = start_bin:stop_bin % now loop through all affected time bins

    d = bin_bounds_mid(bin_i); % the distance to bin center
    
    % ellipse parameters
    a = 0.5*sqrt(d^2 - m^2);
    b = d/2;
    c = a;
    
    % find A and B for end of time bin
    dmax = bin_bounds_mid(bin_i) + bin_size_m/2;
    amax = 0.5*sqrt(dmax^2 - m^2);
    bmax = dmax/2;
    cmax = amax;
    [Amax,Bmax,~,~] = ellipsoidPlaneIntersection(dmax,amax,bmax,cmax,q,perp_to_plane_n,plane_n,l,w);
    
    if bin_i==start_bin % first bin
      
        % find A and B for start of time bin
        dmin = bin_bounds_mid(bin_i) - bin_size_m/2;
        amin = 0.5*sqrt(dmin^2 - m^2);
        bmin = dmin/2;
        cmin = amin;
        [Amin,Bmin,~,~] = ellipsoidPlaneIntersection(dmin,amin,bmin,cmin,q,perp_to_plane_n,plane_n,l,w);
 
    end
    
    [A,B,o,R] = ellipsoidPlaneIntersection(d,a,b,c,q,perp_to_plane_n,plane_n,l,w);
    
    if plotFlag == 1
        rpol = A*B./sqrt((B*cos(theta)).^2 +(A*sin(theta)).^2);
        
        % compute ellipse in cart coords, for plotting
        ellipse_r = rpol.*cos(theta);
        ellipse_s = rpol.*sin(theta);
    end
    
    vertices_rsu = xyz2rsu(R,o,vertices,size(vertices,1));
    
    if plotFlag == 1
        % plot the rsu plane
        figure(20)
        patch(vertices_rsu(:,1),vertices_rsu(:,2),'blue', 'FaceAlpha', .1)
        
        hold on
        plot(ellipse_r,ellipse_s)
        grid on
        xlabel('r')
        ylabel('s')
    end
    
    
    
    if plotFlag == 1
        
        % convert ellipse to xyz for plotting
        ellipse_rsu = [ellipse_r',ellipse_s',zeros(nAngle,1)];
        
        % plot ellipse in xyz 3D plot
        figure(77)
        hold on
        ellipse_xyz = rsu2xyz(R,o,ellipse_rsu,nAngle);
        scatter3(ellipse_xyz(:,1), ellipse_xyz(:,2), ellipse_xyz(:,3))
    end
    
    
    % find angular bounds
    [angleBounds,nIntegralVec,nseg] = findAngularBounds_verticalFacet_condition_par(...
        vertices_rsu,A,B,plotFlag,maxDist,params); % bigger angle comes first
    
    for ns = 1:nseg
        nIntegral = nIntegralVec(ns);
        minAngle = min(angleBounds(:,ns));
        maxAngle = max(angleBounds(:,ns));
        
        deltaAngle = (maxAngle-minAngle)/nIntegral;
        midAngles = linspace(minAngle + deltaAngle/2, maxAngle - deltaAngle/2,nIntegral);
        startAngles = midAngles - deltaAngle/2;
        stopAngles = midAngles + deltaAngle/2;
        
        
%          fun = @(x) A*B*sqrt((A^4*sin(x).^2+B^4*cos(x).^2)./(A^2*sin(x).^2+B^2*cos(x).^2).^3);
        
        for angleI = 1:nIntegral% maybe speed this up?
            midAngle = midAngles(angleI);
            startAngle = startAngles(angleI);
            stopAngle = stopAngles(angleI);
            
            midRpol = A*B./sqrt((B*cos(midAngle)).^2 +(A*sin(midAngle)).^2);
            midP_rsu = [midRpol*cos(midAngle),midRpol*sin(midAngle),0];
            midP_xyz = rsu2xyz(R,o,midP_rsu,1);
            
            
            % figure out range extent
            thetaMid = midAngle;
            
            
%              estimated_delta_range = Amax*Bmax./sqrt((Bmax*cos(thetaMid)).^2 +(Amax*sin(thetaMid)).^2)...
%                 - Amin*Bmin./sqrt((Bmin*cos(thetaMid)).^2 +(Amin*sin(thetaMid)).^2);

            
%             if bin_i>start_bin
%                 estimated_delta_range = A*B./sqrt((B*cos(thetaMid)).^2 +(A*sin(thetaMid)).^2)...
%                     - Apast*Bpast./sqrt((Bpast*cos(thetaMid)).^2 +(Apast*sin(thetaMid)).^2);
%             else % just use 0 - its probably super small anyway. Fix later.
%                 estimated_delta_range = 1;
%             end
            
            midVal = (n_w*( -(w-midP_xyz) )')*(plane_n*(w-midP_xyz)')*...
                (plane_n*(l-midP_xyz)')*(n_l*( -(l-midP_xyz) )')/(norm(w-midP_xyz)^4*norm(l-midP_xyz)^4);
            
            if plotFlag == 1
                figure(77)
                hold on
                % plot the center point
                scatter3(midP_xyz(1), midP_xyz(2), midP_xyz(3),'filled')
            end
            
            
            
            
            % compute annulus area
            Rmax = Amax*Bmax./sqrt((Bmax*cos(thetaMid)).^2 +(Amax*sin(thetaMid)).^2);
            Rmin = Amin*Bmin./sqrt((Bmin*cos(thetaMid)).^2 +(Amin*sin(thetaMid)).^2);
            deltaAngle = abs(stopAngle-startAngle);
            annulusSize = (deltaAngle/2)*(Rmax^2-Rmin^2);
            histogram(bin_i) = histogram(bin_i)+midVal*annulusSize;
            
        end
    end

    Amin = Amax;
    Bmin = Bmax;

    
%     Apast = A;
%     Bpast = B;
end


end

