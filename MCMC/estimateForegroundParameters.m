function [rpast,apast,hpast,phiPast,hist,proposedRates] = estimateForegroundParameters(params,meas,...
    init,dist,setup,max_iter,stdMult,adapt_after,objects, indices2Dthrow,nshift,truth)
% this function assumes the only thing that changes is in the foreground!
% There will be some model mismatch because this is approximately true
% backgroundRate - everything that was there before introduction of
% foreground
global dirpath

l_xyz0 = params.l_xyz0;
posNum = params.posNum;
nPix = params.nPix;
num_time_bins = params.num_time_bins;


plotFlag = 0;


% init - initial values
rpast = init.rpast;
apast = init.apast;
hpast = init.hpast;
spast = init.spast;
phiPast = init.phiBounds;

% dist - proposal and prior distribution parameters
rmax = dist.rmax;
rmin = dist.rmin;
amin = dist.amin;
amax = dist.amax;
hmin = dist.hmin;
hmax = dist.hmax;
phi1max = dist.phi1max;
phi1min = dist.phi1min;
phi2max = dist.phi2max;
phi2min = dist.phi2min;

nTarget = size(phiPast,1);


std_prop_r = dist.std_prop_r;
std_prop_a = dist.std_prop_a;
std_prop_h = dist.std_prop_h;
std_prop_phi = dist.std_prop_phi;

laser_intensity = setup.laser_intensity;
backgroundRate = setup.backgroundRate;
backgroundRate(indices2Dthrow,:) = [];
meas(indices2Dthrow,:) = [];

% truth, if we have it
if ~isempty(truth)
    trueRange = truth.trueRange;
    trueAlbedo = truth.trueAlbedo;
    trueHeight = truth.trueHeight;
end

% use initial values to propose first object
proposedObject = makeProposedObject_multiTarget(phiPast, hpast, ...
    rpast);

pastRates = zeros(nPix,num_time_bins);
for t = 1:nTarget
    thisTarget = proposedObject(t,:);
    [thisTargetRates,~,~,~] = computeForwardModel_adaptiveSegmentSize_par(thisTarget,params,plotFlag);
    pastRates = pastRates + apast(t)*thisTargetRates;
end

pastRates = circshift(pastRates,nshift,2);
pastRates(indices2Dthrow,:) = [];
pastRates = laser_intensity*pastRates + spast*backgroundRate;

pastRatesNoZeros = pastRates;
pastRatesNoZeros(pastRatesNoZeros<1)=1;
log_f_x_given_theta_past =   sum(sum(meas.*log(pastRatesNoZeros) - pastRates - gammaln(meas + 1))) ;
 


n_accept = 0;
n_throw = 0;
n_total = 0;

rhist = [];
ahist = [];
hhist = [];
phi1hist = [];
phi2hist = [];

% keep track of current acceptance rate
last100 = zeros(1,100);

% figure handles
fgResult = figure(3);
fgProgress = figure(4);
for it = 1:max_iter
    
     % add an adaptive component
    if (it>adapt_after) && (mod(it,100)==1)
        if sum(last100)<20
            disp('Reducing proposal variance to increase acceptance rate.')
            std_prop_r = stdMult*std_prop_r;
            std_prop_a = stdMult*std_prop_a;
            std_prop_h = stdMult*std_prop_h;
            std_prop_phi = stdMult*std_prop_phi;
        elseif sum(last100)>40
            disp('Increasing proposal variance to decrease acceptance rate.')
            std_prop_r = 1/stdMult*std_prop_r;
            std_prop_a = 1/stdMult*std_prop_a;
            std_prop_h = 1/stdMult*std_prop_h;
            std_prop_phi = 1/stdMult*std_prop_phi;
        end
    end
    
    
    % generate candidate
    rprop = std_prop_r.*randn(nTarget,1) + rpast;
    aprop = std_prop_a.*randn(nTarget,1) + apast;
    hprop = std_prop_h.*randn(nTarget,1) + hpast;
    phiProp = std_prop_phi.*randn(nTarget,2) + phiPast;

    % to save time do this first
    % evaluate prior
    g_theta = (prod(rprop<rmax)&&prod(rprop>rmin)&&prod(aprop<amax)&&prod(aprop>amin)...
        &&prod(hprop<hmax)&&prod(hprop>hmin)...
        &&prod(phi1max>phiProp(:,1))&&prod(phi1min<phiProp(:,1))&&prod(phi2max>phiProp(:,2))...
        &&prod(phi2min<phiProp(:,2))); % indicator function
  
    if g_theta~= 0
        % form the corresponding object structure
        [proposedObject] = makeProposedObject_multiTarget(phiProp,hprop, ...
            rprop);
        
        % compute photon rates corresponding to this object
        proposedRates = zeros(nPix,num_time_bins);
        for t = 1:nTarget
            thisTarget = proposedObject(t,:);
            [thisTargetRates,~,~,~] = computeForwardModel_adaptiveSegmentSize_par(thisTarget,params,plotFlag);
            proposedRates = proposedRates + aprop(t)*thisTargetRates;
        end

        proposedRates = circshift(proposedRates,nshift,2);
        proposedRates(indices2Dthrow,:) = [];
        proposedRates = laser_intensity*proposedRates+spast*backgroundRate;

        proposedRatesNoZeros = proposedRates;
        proposedRatesNoZeros(proposedRatesNoZeros<1) = 1;
        
        % evaluate likelihood
        log_f_x_given_theta =   sum(sum(meas.*log(proposedRatesNoZeros) - proposedRates - gammaln(meas + 1))) ;      
        
        % compute acceptance ratio
        alpha = g_theta*exp( log_f_x_given_theta -  log_f_x_given_theta_past);
        
        u = rand; % uniform random number
        
        if u<=alpha % accept
            % update more recent count
            last100(1:99) = last100(2:100); % get rid of oldest entry, i.e. the first one
            last100(end) = 1; % update last entry with most recent result
            
            n_accept = n_accept + 1;
            n_total = n_total + 1;
            
            % update values
            rpast = rprop;
            rhist = [rhist, rprop];
            apast = aprop;
            ahist = [ahist, aprop];
            hpast = hprop;
            hhist = [hhist, hprop];

            
            phi1hist = [phi1hist,phiProp(:,1)];
            phi2hist = [phi2hist,phiProp(:,2)];
            phiPast = phiProp;

            
            log_f_x_given_theta_past = log_f_x_given_theta;
            
            
%             figure(3)
            set(0,'CurrentFigure',fgResult)
            subplot(321)
            hold off
            plot(rhist','b--o')
            hold on
            if ~isempty(truth)
                plot(trueRange*ones(size(rhist)))
            end
            title({['Range draw'];['${n_{\rm accept}} = $' num2str(sum(last100))]},'Interpreter','latex','FontSize',16)
            xlabel('index of accepted iterations','Interpreter','latex','FontSize',14)
            ylabel('r','Interpreter','latex','Interpreter','latex','FontSize',14)
            grid on
            
            
            subplot(322)
            hold off
            plot(ahist','b--o')
            hold on
            if ~isempty(truth)
                plot(trueAlbedo*ones(size(ahist)))
            end
            title({['Albedo draw'];['$n_{\rm total} = $' num2str(n_total)]},'Interpreter','latex','FontSize',16)
            xlabel('index of accepted iterations','Interpreter','latex','FontSize',14)
            ylabel('a','Interpreter','latex','FontSize',14)
            grid on
            
            subplot(323)
            hold off
            plot(hhist','b--o')
            hold on
            if ~isempty(truth)
                plot(trueHeight*ones(size(hhist)))
            end
            title({['Height draw']},'Interpreter','latex','FontSize',16)
            xlabel('index of accepted iterations','Interpreter','latex','FontSize',14)
            ylabel('h','Interpreter','latex','FontSize',14)
            grid on
            
            
            
            
            subplot(325)
            hold off
            plot((180/pi)*phi1hist(1,:),'b--o')
            hold on
            plot((180/pi)*phi2hist(1,:),'b--o')
            title({['phi target 1']},'Interpreter','latex','FontSize',16)
            xlabel('index of accepted iterations','Interpreter','latex','FontSize',14)
            ylabel('phi','Interpreter','latex','FontSize',14)
            grid on
            
            if nTarget>1
                subplot(326)
                hold off
                plot((180/pi)*phi1hist(2,:),'b--o')
                hold on
                plot((180/pi)*phi2hist(2,:),'b--o')
                title({['phi target 2']},'Interpreter','latex','FontSize',16)
                xlabel('index of accepted iterations','Interpreter','latex','FontSize',14)
                ylabel('phi','Interpreter','latex','FontSize',14)
                grid on
            end
            
            set(0,'CurrentFigure',fgProgress)
            plotObjects(proposedObject,'blue')
            if ~isempty(objects)
                plotObjects(objects,'red')
            end
            scatter3(l_xyz0(1),l_xyz0(2),l_xyz0(3),'r','filled')

            hold off
            
            drawnow
            
        else % reject
            % the last entry of last100 is the most recent
            last100(1:99) = last100(2:100); % get rid of oldest entry, i.e. the first one
            last100(end) = 0; % update last entry with most recent result
            
            n_throw = n_throw + 1;
            n_total = n_total + 1;
            acceptance_ratio = n_accept/n_total;
            set(0,'CurrentFigure',fgResult)
            subplot(321)
            title({['Range draw'];['${n_{\rm accept}} = $' num2str(sum(last100))]},'Interpreter','latex','FontSize',16)
            grid on
            
            subplot(322)
            title({['Albedo draw'];['$n_{\rm total} = $' num2str(n_total)]},'Interpreter','latex','FontSize',16)
            grid on
            
            subplot(323)
            title({['Height draw']},'Interpreter','latex','FontSize',16)
            grid on
            drawnow
        end
    else
        % the last entry of last100 is the most recent
        last100(1:99) = last100(2:100); % get rid of oldest entry, i.e. the first one
        last100(end) = 0; % update last entry with most recent result
        
        n_throw = n_throw + 1;
        n_total = n_total + 1;
        acceptance_ratio = n_accept/n_total;
    end

    
end

set(0,'CurrentFigure',fgResult)
set(gcf,'PaperUnits','centimeters'); 
set(gcf,'PaperSize',4*[7 7]); 
fig = gcf; 
fig.PaperUnits = 'centimeters';  
fig.PaperPosition = 4*[0 0 7 7]; 
fig.Units = 'centimeters'; 
fig.PaperSize=4*[7 7]; 
fig.Units = 'centimeters'; 
print(fig,[dirpath 'ProposalHistoryForeground_' num2str(posNum)],'-dpdf','-r200'); 


hist.rhist = rhist;
hist.ahist = ahist;
hist.hhist = hhist;


end

