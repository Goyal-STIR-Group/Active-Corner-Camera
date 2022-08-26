

varMult =.3;
std_prop_r = varMult*.2*2;
std_prop_a = varMult*1;
std_prop_h = varMult*2*1.5;
std_prop_s = varMult*.005;%.1
std_prop_phi = varMult*.05;% varMult*.002;

% set up broad reasonable bounds on parameter estimates
% range bounds [m]
rmax = 5;
rmin = .2;

% albedo bounds
amin = .3;
amax = 3;

% height bounds [m]
hmin = .5;
hmax = 3;


% initial values
rpast = range_init_foreground*ones(nTarget,1)-.5;
apast = 1*ones(nTarget,1);
hpast = .5*ones(nTarget,1);
spast = scaleFactorInit;


% parameters
init.rpast = rpast;
init.apast = apast;
init.hpast = hpast;
init.spast = spast;
init.phiBounds = phi_bounds_init;


dist1.rmax = rmax;
dist1.rmin = rmin;
dist1.amax = amax;
dist1.amin = amin;
dist1.hmax = hmax;
dist1.hmin = hmin;
% dist1.smin = smin;
% dist1.smax = smax;

[phi1min,phi1max,phi2min,phi2max] = findPhiBounds(phi_bounds_init,nTarget);


dist1.phi2max = phi2max;
dist1.phi2min = phi2min;
dist1.phi1max = phi1max;
dist1.phi1min = phi1min;

dist1.std_prop_r = std_prop_r;
dist1.std_prop_a = std_prop_a;
dist1.std_prop_h = std_prop_h;
dist1.std_prop_s = std_prop_s;
dist1.std_prop_phi = std_prop_phi;


% trueBackgroundRange = 2.3;%objects{iBackground,2}(2);
% trueBackgroundAlbedo = 1;
% truth.trueRange = trueRange;
% truth.trueHeight = trueHeight;
% truth.trueAlbedo = 1;
% truth.trueBackgroundRange = trueBackgroundRange;
% truth.trueBackgroundAlbedo = trueBackgroundAlbedo ;
% truth.trueB = 0; % tilt from edge pointing normal

setup.laser_intensity = laser_intensity;
setup.backgroundRate = backgroundRate;

std_mult = .3;
adapt_after = 99;
max_iter = 800;