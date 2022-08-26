threshMult = 6;

angle_window = 3*pi/180; % [rad]
thresh_mult_time = .5;%20;


scene_recon_dim = 64;   %number of angles
% parameters for determining rough TIME extent of foreground object
window_grow = 5; % the number of time bins wider we go to make sure we didnt miss any content with threshold

% parameters for forward model A matrix
hidden_scene_height = .5; % meter
nds = 3; % how much down sampling we need to do on floor pixel to smooth out penumbra forward model

% corner camera parameters
nIter = 10000;
stepsize = 10e-7;

% object counting parameters
minSeparation = 4; % threshold crossings that are closer than minSeparation angular bins are merged into a single object

% populate passiveParams
passiveParams.minSeparation = minSeparation;
passiveParams.thresh_mult_time = thresh_mult_time;
passiveParams.window_grow = window_grow;
passiveParams.hidden_scene_height = hidden_scene_height;
passiveParams.nds = nds;
passiveParams.lam2  = lam2;
passiveParams.nIter = nIter;
passiveParams.stepsize = stepsize;
passiveParams.angle_window = angle_window;
passiveParams.threshMult = threshMult;
passiveParams.anglePad = 15*pi/180; % how close to the boundaries (i.e. 0 and pi) to allow threshold crossingns