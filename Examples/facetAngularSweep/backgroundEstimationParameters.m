
%     % initialize
init.rpast = rhat;
init.apast = ahat;
init.hpast = hhat;
init.rbpast = max(range_init_background-.5,max(rhat)+.5)*ones(nTarget,1);
init.abpast =  1*ahat; % background albedo
init.bpast = 0*ones(nTarget,1); % background tilt, 0 because we assume its facing the edge
init.spast = scaleFactorInit;
init.phiBounds  = phiBoundsHat;

% proposal density
varMult = .1;
std_prop_rb = varMult*2;
std_prop_ab = varMult*8;



% prior dist
abmin = .1;
abmax = 2;

rbmin = max(rhat+.1,1);
rbmax = 7;


dist2.std_prop_rb = std_prop_rb;
dist2.std_prop_ab = std_prop_ab;
dist2.abmin = abmin;
dist2.abmax = abmax;
dist2.rbmin = rbmin;
dist2.rbmax = rbmax;


max_iter =600;

std_mult = .3;
adapt_after = 200;

globalParams.twoBackgroundRegionFlag = 1;