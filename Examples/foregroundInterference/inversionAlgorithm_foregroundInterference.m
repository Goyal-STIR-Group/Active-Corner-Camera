clc; clear all; close all;



% set up paths to code
codepath = '/Users/sheilawerth/Dropbox/Lab/ActiveCornerCamera/SheilaHooverShare/Experimental Results/paperCode/';
dirpath = '/Users/sheilawerth/Dropbox/Lab/ActiveCornerCamera/SheilaHooverShare/Experimental Results/paperCode/Examples/foregroundInterference/';
addpath(genpath(codepath))

declareGlobalVariables


load('data_020322_params.mat') % the parameters file for this collection
load('hotPixels.mat') % a saved list of disfunctional pixels in the array
load('roomGroundTruth.mat') % the ground truth stationary hidden scene


%% plot the ground truth stationary scene

figure(4)
if ~isempty(objects)
    plotObjects(objects,'red')
end
hold off
%% set up some parameters
% call script to setup experimental parameters
setupExperimentalParams


numBackgroundInt = 15; 
% numBackgroundInt - The number of background measurements to integrate. 
% Each background measureement has a 2 min acquisition time (acquisition 
% time is defined in Supp. material)

nsec = 30; 
% The length, in acquisition time, of each new measured frame

nFrameUse = nsec/60; % the fraction of a minute [acquisition time] used to acquire each new measurement frame
laser_intensity = laser_intensity*nFrameUse;

cond = ['_nsec' num2str(nsec) '_bg' num2str(numBackgroundInt) '/']; % the name of this test condition
dirpathres = [dirpath, 'Results', cond]; % save results here
dirpath = [dirpath, 'Figures', cond]; % save interim figures here

%% load data
    
load('data_031122_StationaryWideFacet_EmptyRoom.mat')
load('data_foregroundInterference.mat')
    
% posNum = 1: object at 100cm
% posNum = 2: object at 125cm 
% posNum = 3: object at 150cm
    
for posNum =1:3 % loop through positions
    
    
    
    close all
    clc
      
    % load data and adjust parameters that vary based upon the length of
    % the frame

    y_meas = data{posNum}.t30;
    lam2 = .000001;

    
   
    % reshape the measurement
    meas = reshapeData(y_meas,cam_pixel_dim,num_time_bins);
    
    % estimate the rates due to the stationary hidden scene, using the
    % initial background measurement
    backgroundRate = nFrameUse*reshapeData(sum(y_empty(:,(end-numBackgroundInt+1):end),2),...
        cam_pixel_dim,num_time_bins)./(2*numBackgroundInt);


    % compute the scale factor that describes the change in laser intensity
    % between the stationary scene measurement and the 
    scaleFactorInit = sum(sum(meas(indices2D,1:15)))/sum(sum(backgroundRate(indices2D,1:15)));
    meas_dif_scale = meas - scaleFactorInit*backgroundRate;
    

    % produce some plots of the measurements
    plotMeasurements
    
    %% parameter initialization using the passive corner camera
    
    passiveCornerCameraParameters
       
    tic
    rng(8)
    [phi_bounds_init,range_init_foreground,range_init_background] = ...
        parameterInitialization(meas_dif_scale,camera,...
        scene_recon_dim,passiveParams,indices2D, indices2Dthrow);
    toc
    
    nTarget = size(phi_bounds_init,1);
  

    
    %% foreground parameter estimation step

    foregroundEstimationParameters
     
    
    tic
    rng(8) % make results repeatable
    [rhat,ahat,hhat,phiBoundsHat,hist_foreground,proposedForegroundRates] = ...
        estimateForegroundParameters(globalParams,meas,init,dist1,setup,max_iter,...
        std_mult,adapt_after,objects, indices2Dthrow,nshift,[]);
    toc
    
    
    plotForegroundEstimation
    
    
    
    %% Now find background

    
    backgroundEstimationParameters
    
    tic
    rng(9) % make results repeatable
    
    [abpast,rbpast,hist,modeEst] =...
        estimateBackgroundParameters( meas,globalParams,dist2,init,setup,max_iter,objects,...
        std_mult, adapt_after,nshift,indices2D);
    
    toc



    
    
    %%
    
    save([dirpathres '/Facet_pos' num2str(posNum) '.mat'],'objects','camera',...
        'phiBoundsHat','hist','modeEst','hist_foreground','dist1','dist2')
end

%