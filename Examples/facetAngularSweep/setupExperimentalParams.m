roomHeightMax = 3; % the max allowable room height 
transparency = .1; % transparency of facets in plots

plotFlag = 0;
% plotFlag - turn to 1 to trouble shoot the fast forward model computation.
% Best to otherwise keep set to 0

maxDist = .5; 
% maxDist - corresponds to dmax defined in Section 1.2 of the Supplementary 
% material. Larger values result in a courser discretization of hidden 
% facets in our fast forward model computation

nshift = 6;
% nshift - number of time bins(positive numbers shift model later in time)
% required to shift measurement histograms so that the first bounce is at
% time zero


epsilon = 1e-8; % a very small reference number
c_light = 299792458; %Speed of light


% define occluding edge
wall_v0 = [0,0,0]; % location of the bottom of the occluding edge
n_wall = [1,0,0];  % unit vector pointing away from corner, parallel to wall
wall_v1 = [-5,0,0]; % a vector pointing along the extent of the wall; used for plotting the wall

% SPAD Camera Parameters
n_w = [0,0,1]; % floor normal; we assume the camera is pointed at the floor. Do not change.
bin_size = params.bin_size; % time resolution of the sensor
fovWidth = params.camera_FOV; % width of the camera FOV on the ground

% populate camera structure
camera.cam_span_x = fovWidth;
camera.cam_span_y = fovWidth;
camera.num_time_bins = params.num_time_bins;
camera.cam_pixel_dim = params.cam_pixel_dim;
camera.cam_pixel_corner = [params.camera_FOV_center(1) - camera.cam_span_x/2,...
    params.camera_FOV_center(2) - camera.cam_span_y/2,0];
camera.bin_size = params.bin_size;

cam_span_x = camera.cam_span_x;
cam_span_y = camera.cam_span_y;
num_time_bins = camera.num_time_bins;
cam_pixel_dim = camera.cam_pixel_dim;
cam_pixel_corner = camera.cam_pixel_corner;
fovCenter = cam_pixel_corner + [cam_span_x/2,cam_span_y/2,0];


bin_size_m = c_light*bin_size; % camera resolution in meter
pix_area = cam_span_x*cam_span_y;
bin_bounds_start = linspace(0,bin_size_m*(num_time_bins-1),num_time_bins);


% Laser
l_xyz0 = params.laser_pos; % position of laser spot
n_l = [0,0,1]; % laser spot normal
laser_intensity = .2*4e4;


% Compute the location of each camera pixel
nPix = cam_pixel_dim^2;
pixel_size_x = cam_span_x./cam_pixel_dim;
pixel_size_y = cam_span_y./cam_pixel_dim;
pixel_x = linspace(cam_pixel_corner(1) + pixel_size_x/2,...
    cam_pixel_corner(1) + cam_span_x - pixel_size_x/2,cam_pixel_dim);
pixel_y = linspace(cam_pixel_corner(2) + pixel_size_y/2,...
    cam_pixel_corner(2) + cam_span_y - pixel_size_y/2,cam_pixel_dim);

[X, Y] = meshgrid(pixel_x, pixel_y);
X = X(:);
Y = Y(:);

pixel_mid_vec = -(wall_v0 - [X,Y,zeros(nPix,1)]); % the center of each pixel
norm_pixel_mid_vec = sqrt(sum(pixel_mid_vec(:,1).^2 + pixel_mid_vec(:,2).^2 + pixel_mid_vec(:,3).^2 , 2));
theta_mid_vec = acos((pixel_mid_vec*n_wall')./norm_pixel_mid_vec); % azimuthal angle of the center of each pixel, measured clockwise around the corner

% set up time bin bounds
bin_size_m = c_light*bin_size; % bin size in meter
bin_bounds_mid = linspace(0+bin_size_m/2,bin_size_m*(num_time_bins-1)+bin_size_m/2,num_time_bins);


globalParams.num_time_bins = num_time_bins;
globalParams.n_wall = n_wall;
globalParams.wall_v0 = wall_v0;
globalParams.wall_v1 = wall_v1;
globalParams.X = X;
globalParams.Y = Y;
globalParams.theta_mid_vec = theta_mid_vec;
globalParams.transparency = transparency;
globalParams.nPix = nPix;
globalParams.l_xyz0 = l_xyz0;
globalParams.maxDist = maxDist;
globalParams.epsilon = epsilon;
globalParams.c_light = c_light;
globalParams.bin_size = bin_size;
globalParams.bin_bounds_mid = bin_bounds_mid;
globalParams.bin_size = bin_size;
globalParams.n_w = n_w;
globalParams.n_l = n_l;
globalParams.dirpath = dirpath;
globalParams.posNum = posNum;
