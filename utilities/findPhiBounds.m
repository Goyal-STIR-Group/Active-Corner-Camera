function [phi1min,phi1max,phi2min,phi2max] = findPhiBounds(phi_bounds_init,nTarget)
    
    % phi1 is used here to describe the smallest angle in a given targets
    % angular extent. phi2 describes teh largest angle in a targets angular
    % extent
    % The vectors phi1min, phi1max, phi2min, and phi2max are of length n
    % where n is the number of detected objects. These vectors contain the
    % min or max (depending on name) allowable value for phi1 or phi2 (depending on
    % name), for each detected object
    
    phi1min = zeros(nTarget,1);
    phi1max = zeros(nTarget,1);
    phi2min = zeros(nTarget,1);
    phi2max = zeros(nTarget,1);
    
    
    
    for t = 1:nTarget
        if t == 1
            phi1min(t) = 0;
        else
            phi1min(t) = phi_bounds_init(t,1) - (phi_bounds_init(t,1)-phi_bounds_init(t-1,2))/2;
        end
        
        targetWidthInit = phi_bounds_init(t,2)-phi_bounds_init(t,1);
        
        phi1max(t) = phi_bounds_init(t,1)+targetWidthInit/2;
        phi2min(t) = phi_bounds_init(t,2)-targetWidthInit/2;
        
        if t ==  nTarget
            phi2max(t) = pi;
        else
            phi2max(t) = phi_bounds_init(t,2) + (phi_bounds_init(t+1,1)-phi_bounds_init(t,2))/2;
        end
    end
end

