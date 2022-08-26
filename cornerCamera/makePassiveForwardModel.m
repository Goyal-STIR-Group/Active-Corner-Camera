function [A] = makePassiveForwardModel(cam_pixel_dim,cam_span_x,...
    cam_span_y,cam_pixel_corner,nDS,hidden_scene_height,hidden_scene_range,...
    scene_recon_dim)
% pixel centers
delta_floor = cam_span_x/cam_pixel_dim;
delta_ds = delta_floor/nDS;
n_meas_pix = cam_pixel_dim^2;
pixel_size_x = cam_span_x./cam_pixel_dim; 
pixel_size_y = cam_span_y./cam_pixel_dim;
pixel_x = linspace(cam_pixel_corner(1) + pixel_size_x/2,...
    cam_pixel_corner(1) + cam_span_x - pixel_size_x/2,cam_pixel_dim);
pixel_y = linspace(cam_pixel_corner(2) + pixel_size_y/2,...
    cam_pixel_corner(2) + cam_span_y - pixel_size_y/2,cam_pixel_dim);

[X, Y] = meshgrid(pixel_x, pixel_y);
X = X(:);
Y = Y(:);
floor_coords = [X,Y,zeros(size(X))];
phi = atan2(Y,X);

% now the hidden side
n_recon_pix = scene_recon_dim;
% hidden_scene_height = 1;
% hidden_scene_range = 2;

% find info on hidden scene pixels
delta_theta = pi/n_recon_pix;
theta = delta_theta/2:delta_theta:(pi-delta_theta/2);
x_scene = hidden_scene_range*cos(theta);
y_scene = hidden_scene_range*sin(theta);
z_scene = hidden_scene_height *ones(1,n_recon_pix);
scene_coords = [x_scene',y_scene',z_scene'];

A = zeros(n_meas_pix, n_recon_pix);

for floorI = 1:n_meas_pix
    floor_pixel_center = floor_coords(floorI,:);
    floor_coord_x = linspace((floor_pixel_center(1)-delta_floor/2+delta_ds/2),...
        (floor_pixel_center(1)+delta_floor/2+delta_ds/2),nDS);
    floor_coord_y = linspace((floor_pixel_center(2)-delta_floor/2+delta_ds/2),...
        (floor_pixel_center(2)+delta_floor/2+delta_ds/2),nDS);
    [X_ds, Y_ds] = meshgrid(floor_coord_x, floor_coord_y);
    X_ds = X_ds(:);
    Y_ds = Y_ds(:);
    this_phi = atan2(Y_ds,X_ds);
%     this_phi = phi(floorI);
    for sceneI = 1:n_recon_pix
        scene_coord = scene_coords(sceneI,:);
        this_theta = theta(sceneI);
        ifilt = find((abs(this_phi)+this_theta)<pi);
        if ~isempty(ifilt)% something is not occluded
            this_range = norm(scene_coord-floor_pixel_center);
            A(floorI,sceneI) = length(ifilt)/this_range^2;
        end

    end
end
A = A./nDS^2;

end

