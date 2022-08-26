function [phiBounds, range_init_foreground,range_init_background] = parameterInitialization(frame1_dif,camera,scene_recon_dim,passiveParams,...
    indices2D,indices2Dthrow)
global c_light
global dirpath
global posNum

% load parameters
anglePad = passiveParams.anglePad;
cam_pixel_dim = camera.cam_pixel_dim;
num_time_bins = camera.num_time_bins;

% parameters for determining rough TIME extent of foreground object
thresh_mult_time = passiveParams.thresh_mult_time;
window_grow = passiveParams.window_grow;% the number of time bins wider we go to make sure we didnt miss any content with threshold

% parameters for forward model A matrix
hidden_scene_height = passiveParams.hidden_scene_height;
nds = passiveParams.nds;


% corner camera parameters
lam2 = passiveParams.lam2;
nIter = passiveParams.nIter;
stepsize = passiveParams.stepsize;

angle_window = passiveParams.angle_window;
threshMult = passiveParams.threshMult;

% parameters for counting targets
minSeparation = passiveParams.minSeparation; % min separation between targets in units of bins


bin_size = camera.bin_size;
bin_size_m = c_light*bin_size; % bin size in meter

% spatially averaged measurement
frame1_spatial_sum = sum(frame1_dif(indices2D,:),1)/cam_pixel_dim^2;
sigPow = max(frame1_spatial_sum)^2;
lam2 = sigPow*lam2;



%% Working with SPATIALLY averaged data, pull out time bins of interest

thresh_time = thresh_mult_time*max((frame1_spatial_sum)); % threshold for foreground
thresh_time_neg = thresh_mult_time*min((frame1_spatial_sum)); % threshold for background


% find start time bin of background
frame1_spatial_sum_smooth = filter([1/3,1/3,1/3], 1, frame1_spatial_sum); % added for extremely noisy data
i_background = find(frame1_spatial_sum_smooth==min(frame1_spatial_sum_smooth));
i_background_start = min(find(frame1_spatial_sum<thresh_time_neg)); % index where background starts

% find time bin range of foreground
thresh_passes = find(frame1_spatial_sum>thresh_time);
thresh_passes(thresh_passes>i_background) = [];
start_time_bin = max(min(thresh_passes)-window_grow,1);
stop_time_bin = min(min(max(thresh_passes)+window_grow,num_time_bins),i_background_start);
stop_time_bin = max(stop_time_bin,max(thresh_passes));
% the strategy here is to look after the smallest threshold crossing, and
% stop looking when we start seeing negative difference

% find approximate ranges to initialize facets at
range_init_foreground = (mean([start_time_bin, stop_time_bin])*bin_size_m)/2;
range_init_background = i_background*bin_size_m/2;

% for plotting purposes only:
thresh_profile = zeros(1,num_time_bins); % all elements greater than the threshold
thresh_profile(thresh_passes) = 2*thresh_time; % 
thresh_pad_profile = zeros(1,num_time_bins); %
thresh_pad_profile(start_time_bin:stop_time_bin)= 2.5*thresh_time;

% plot spatially averaged data and show threshold crossings
figure(31)
hold on
grid on
axis square
% xlim([0 bin_size_m*100])
xlabel('Range [m]','Interpreter','latex','FontSize',12 )
plot(bin_size_m*(1:num_time_bins),frame1_spatial_sum)
plot(bin_size_m*(1:num_time_bins),thresh_time*ones(1,num_time_bins))
plot(bin_size_m*(1:num_time_bins),thresh_profile,'g','LineWidth',1.5)
plot(bin_size_m*(1:num_time_bins),thresh_pad_profile,'g','LineWidth',1.5)

title('Spatial Average of Difference Measurement','Interpreter','latex','FontSize',14)
legend('Spatial average','Threshold','Crosses','All included bins')


set(gcf,'PaperUnits','centimeters'); 
set(gcf,'PaperSize',3*[5 4]); 
fig = gcf; 
fig.PaperUnits = 'centimeters';  
fig.PaperPosition = 3*[0 0 5 4]; 
fig.Units = 'centimeters'; 
fig.PaperSize=3*[5 4]; 
fig.Units = 'centimeters'; 
print(fig,[dirpath 'spatialMeasSum_' num2str(posNum)],'-dpdf','-r200'); 
% 

%% Now, form passive measurement by TEMPORALLY integrating over foreground time bins

passiveMeasure = sum(frame1_dif(:,start_time_bin:stop_time_bin),2);
passiveMeasure = passiveMeasure(:)+150;
satLim = max(passiveMeasure(indices2D));
satLimMin = min(passiveMeasure(indices2D));
passiveMeasure(passiveMeasure<0)=0;


% plot the passive measurement
figure(303)
imagesc(reshape(passiveMeasure-150,cam_pixel_dim,cam_pixel_dim),[satLimMin-150 satLim-150])
title('Passive Measurement','Interpreter','latex','FontSize',14)
colormap('gray')
axis off
axis square
colorbar

fig = gcf; 
fig.PaperUnits = 'centimeters';  
fig.PaperPosition = 2*[-.4 -.3 5 4]; 
fig.PaperSize=2*[4.2 3.6]; 
fig.Units = 'centimeters'; 
print(fig,[dirpath 'PassiveMeas_' num2str(posNum)],'-dpdf','-r200'); 


%% Run the passive corner camera


hidden_scene_range = range_init_foreground; % approximate range of hidden scene, used to form passive foward model
angular_res = pi/scene_recon_dim; % angular resolution of hidden scene

% form the forward model
A_penumbra = makePassiveForwardModel(camera.cam_pixel_dim,camera.cam_span_x,...
    camera.cam_span_y,camera.cam_pixel_corner,nds,hidden_scene_height,...
    hidden_scene_range,scene_recon_dim);

A_penumbra(indices2Dthrow,:) = 0; % don't use data from damaged pixels

S =  A_penumbra; % A_fan
ST = S';



% We will promote sparsity of the scene estimate in this wavelet basis:
N2 = size(S,2); % number of scene elements
a = ones(N2,1);
level = fix(log2(N2));
wavename = 'db4';
[Coef,L] = wavedec(a,level,wavename);
N1 = length(Coef);
wavebasis = zeros(N1,N2);
for ii = 1:N2
    x = zeros(N2,1);
    x(ii) = 1;
    wavebasis(:,ii) = wavedec(x,level,wavename);
end

B =@(x) wavebasis*x;
invB =@(theta) wavebasis'*theta;

% run corner camera
[~,v_total] = passive1DCC(S, ST, B, invB, N1, passiveMeasure(:),...
    lam2,nIter,stepsize);

%% Use 1D scene estimate to count number of objects in motion

% compute a threshold
v_total_sum = mean(v_total);
thresh = threshMult*v_total_sum; % threshold used to find objects

% find threshold crossings
maxAngleIndex = floor((pi-anglePad)/angular_res); % don't look at angles past this angle
minAngleIndex = ceil((anglePad)/angular_res); % don't look at angles less than this angle
v_total_pad = v_total;
v_total_pad(1:minAngleIndex) = 0;
v_total_pad(maxAngleIndex:end) = 0;


ifilt_cross = find(v_total_pad>thresh);

if length(ifilt_cross)<1 % make sure we have some crossing, assuming there is some motion in the scene
    disp('no threshold crossing! Motion not detected.')    
end

% merge regions that are too close to eachother in angle
iGapFill = find((diff(ifilt_cross)<minSeparation)& diff(ifilt_cross)>1); %indices of gaps to fill

% fill the gaps
for t = 1:length(iGapFill)
    startGap = ifilt_cross(iGapFill(t));
    stopGap = ifilt_cross(iGapFill(t)+1);
    ifilt_cross = [ifilt_cross; [startGap+1:stopGap-1]'];
end

ifilt_cross = sort(ifilt_cross,'ascend'); % a list of threshold crossings

% remove entries that are one bin wide - don't count those as threshold
% crossings
% COME BACK AND COMMENT THIS OUT IF IT ISNT A PROBLEM
% if length(ifilt_cross) == 1
% %     ifilt_cross = [ifilt_cross; ifilt_cross+1];
% % REMOVE THIS CASE
%     disp('See parameterInitialization, line 220 for a check.')
% else % extend single bin threshold crosses to be two bins wide so that later logic works


% the logic used later requires different that a target not be a single
% angular bin wide. If only one bin crosses the threshold, add a neighbor
i_add = [];
if length(ifilt_cross) == 1
    i_add = [i_add,ifilt_cross+1];
else
    for p = 1:length(ifilt_cross)
        if p == 1 % first entry, only compare to next entry
            if (ifilt_cross(p+1)-ifilt_cross(p))>1 % if they are further apart than one bin
                i_add = [i_add, ifilt_cross(p)+1];
                %                 i_remove = [i_remove,p];
            end
        elseif p == length(ifilt_cross) % last entry, compare only to prior entry
            if (ifilt_cross(p)-ifilt_cross(p-1))>1
                i_add = [i_add,ifilt_cross(p)-1];
                %                 i_remove = [i_remove,p];
            end
        else
            if ((ifilt_cross(p)-ifilt_cross(p-1))>1)&&((ifilt_cross(p+1)-ifilt_cross(p))>1)
                i_add = [i_add,ifilt_cross(p)+1];
                %                 i_remove = [i_remove,p];
            end
        end
    end
    if ~isempty(i_add)
        disp('See parameterInitialization, line 220 for a check.')
    end
end
ifilt_cross = [i_add; ifilt_cross];
ifilt_cross = sort(ifilt_cross,'ascend');
%     ifilt_cross(i_remove) = [];
% end

gapIndices = find(diff(ifilt_cross)>1); % indices where there is a gap between targets
nCross = length(gapIndices)+1;
midIndices = find(((ifilt_cross-circshift(ifilt_cross,1)-1)==0)&...
    ((ifilt_cross-circshift(ifilt_cross,-1)+1)==0)); % indices that are not at the boundaries
ifilt_cross(midIndices) = [];
    

ifilt_cross = reshape(ifilt_cross,2,nCross); ifilt_cross = ifilt_cross';

phiBounds = angular_res*ifilt_cross;
phiBounds(:,1) = phiBounds(:,1) - angle_window;
phiBounds(:,2) = phiBounds(:,2) + angle_window;

% for plotting only
scene_cross_profile = zeros(size(v_total));
scene_cross_profile(ifilt_cross) = thresh;


figure(1111)
plot(linspace(0,180,length(v_total)),v_total)
title('1D Scene Estimate','Interpreter','latex','FontSize',14)
xlabel('Angle [deg]','Interpreter','latex','FontSize',14)
ylabel('1D Scene Estimate','Interpreter','latex','FontSize',14)
set(gca, 'YTickLabel', [])
ax = gca;
ax.XAxis.FontSize = 14;
grid on
xlim([0 180])

fig = gcf;
fig.PaperUnits = 'centimeters';
fig.PaperPosition = 2*[-.3 0 5 4];
fig.Units = 'centimeters';
fig.PaperSize=2*[4.4 3.9];
print(fig,[dirpath '1DsceneEst_' num2str(posNum)],'-dpdf','-r200');



figure(33)
subplot(232)
plot(linspace(0,180,length(v_total)),v_total,'b--o')
hold on
plot(linspace(0,180,length(v_total)),thresh*ones(size(v_total)))
plot(linspace(0,180,length(v_total)),scene_cross_profile,'g','LineWidth',1.5)
% scatter((180/pi)*phi_1,thresh,'r','filled')
% scatter((180/pi)*phi_2,thresh,'r', 'filled')

xlabel('Angle [deg]','Interpreter','latex','FontSize',12)
title('1D Scene Estimate w/ Threshold','Interpreter','latex','FontSize',14)
axis square
xlim([0 180])
grid on
legend('Scene estimate','Threshold','Crosses','Location','southwest')


figure(302)
plot(linspace(0,180,length(v_total)),v_total)
hold on
plot(linspace(0,180,length(v_total)),thresh*ones(size(v_total)))
plot(linspace(0,180,length(v_total)),scene_cross_profile,'g','LineWidth',1.5)
xlabel('Angle [deg]','Interpreter','latex','FontSize',12)
title('1D Scene Estimate','Interpreter','latex','FontSize',14)
axis square
xlim([0 180])
grid on
legend('Scene estimate','Threshold','Crosses')

set(gcf,'PaperUnits','centimeters'); 
set(gcf,'PaperSize',3*[5 4]); 
fig = gcf; 
fig.PaperUnits = 'centimeters';  
fig.PaperPosition = 3*[0 0 5 4]; 
fig.Units = 'centimeters'; 
fig.PaperSize=3*[5 4]; 
fig.Units = 'centimeters'; 
print(fig,[dirpath 'CornerCameraRun_' num2str(posNum)],'-dpdf','-r200'); 



end

