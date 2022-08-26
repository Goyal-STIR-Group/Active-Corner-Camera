function [A,B,o,R] = ellipsoidPlaneIntersection(d,a,b,c,q,perp_to_plane_n,plane_n,l,w)
global epsilon
% check that q is an interior point of the ellipse
%     if norm(l-q)+norm(q-w)>d
%         disp('q is not an interior point of the ellipsoid.')
%     end
    
    
    r = perp_to_plane_n; % pick r to be a vector perpendicular to n - here I just pick a random vector in the plane
    s = cross(r,plane_n);

    D = diag([1/a, 1/b, 1/c]);
%     dot(D*r',D*s')


    % check condition in eq (7)
    if dot(D*r',D*s')>epsilon
        %     disp('Eq 7 not satisfied - make sure everything works.')
        if (dot(D*r',D*r')-dot(D*s',D*s')) > epsilon % i.e. not z
            omega = .5*atan( 2*dot(D*r',D*s') / ( dot(D*r',D*r')-dot(D*s',D*s') ) );
        else
            omega = pi/4;
        end
        r_old = r;
        s_old = s;
        r = cos(omega)*r_old + sin(omega)*s_old;
        s = -sin(omega)*r_old + cos(omega)*s_old;
    end

    
    % find ellipse information
    ellipsoid_params = [a, b, c];
    k = dot(q,plane_n);
    kt = sqrt(sum((plane_n.^2).*(ellipsoid_params.^2)));
    d = dot(D*q',D*q') - ( dot(D*q',D*r')^2 )/dot(D*r',D*r') - ( dot(D*q',D*s')^2 )/dot(D*s',D*s');
    o = (k/kt^2)*(ellipsoid_params.^2).*plane_n; % ellipse center
%     % plot the ellipse center, as computed by Klein
%     scatter3(o(1),o(2),o(3),30,'r','filled')
    
    % semi-axes
    A = sqrt((1-d)/dot(D*r',D*r'));
    B = sqrt((1-d)/dot(D*s',D*s'));
        % new basis vectors
    er = r./norm(r);
    es = s./norm(s);
    eu = cross(er,es); eu = eu./norm(eu);
    
    % rotation matrix
    R = [er; es; eu]; % project
end

