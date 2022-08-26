function [abpast,rbpast,hist,modeEst] = ...
    estimateBackgroundParameters(meas,params,dist,init,setup,...
    max_iter,objects,stdMult, adapt_after, nshift,indices2D)
global dirpath
plotFlag = 0;

% dist - proposal and prior distribution parameters
l_xyz0 = params.l_xyz0;
posNum = params.posNum;
nPix = params.nPix;
num_time_bins = params.num_time_bins;
twoBackgroundRegionFlag = params.twoBackgroundRegionFlag;


abmin = dist.abmin;
abmax = dist.abmax;
rbmin = dist.rbmin;
rbmax = dist.rbmax;


% init - initial values
rprop = init.rpast;
apast = init.apast;
hprop = init.hpast;
rbpast = init.rbpast;
abpast = init.abpast;
bpast = init.bpast;
spast = init.spast;
phiPast = init.phiBounds;


nTarget = size(phiPast,1);



std_prop_rb = dist.std_prop_rb;
std_prop_ab = dist.std_prop_ab;

laser_intensity = setup.laser_intensity;
backgroundRate = setup.backgroundRate;
backgroundRate(indices2D,:) = [];
meas(indices2D,:) = [];

n_accept = 0;
n_throw = 0;
n_total = 0;


rbhist = [];
abhist = [];



% figure handles
fgResult = figure(3);
fgProgress = figure(4);


proposedObject = makeProposedObjectForegroundBackgroundTilt_multiTarget_FOV(phiPast, hprop, ...
    rprop, rbpast, bpast);

pastForegroundRates = zeros(nPix,num_time_bins);
proposedOccludedRates = zeros(nPix,num_time_bins);

for t = 1:nTarget
    
    if twoBackgroundRegionFlag == 1
        % pull out object
        thisProposedTarget = proposedObject((t-1)*3+1:t*3,:);
        thisProposedTargetMerged = mergeRegions(thisProposedTarget);
        [~,~,thisForegroundRates,thisOccludedRates] = ...
            computeForwardModel_adaptiveSegmentSize_par(thisProposedTargetMerged,params,plotFlag);
        pastForegroundRates = pastForegroundRates + apast(t)*thisForegroundRates;
        proposedOccludedRates = proposedOccludedRates + abpast(t)*thisOccludedRates;
    else
        thisProposedTarget = proposedObject((t-1)*3+1,:);
        thisProposedTarget(2,:) = proposedObject(t*3,:);
        
        [~,~,thisForegroundRates,thisOccludedRates] = ...
            computeForwardModel_adaptiveSegmentSize_par(thisProposedTarget,params,plotFlag);
        pastForegroundRates = pastForegroundRates + apast(t)*thisForegroundRates;
        proposedOccludedRates = proposedOccludedRates + abpast(t)*thisOccludedRates;

        
    end
end

pastForegroundRates = circshift(pastForegroundRates,nshift,2);
proposedOccludedRates = circshift(proposedOccludedRates,nshift,2);

pastForegroundRates(indices2D,:) = [];
proposedOccludedRates(indices2D,:) = [];

pastRates = laser_intensity*(pastForegroundRates-proposedOccludedRates)...
    +spast*backgroundRate;


pastRatesNoZeros = pastRates;
pastRatesNoZeros(pastRatesNoZeros<1)=1;
log_f_x_given_theta_past =   sum(sum(meas.*log(pastRatesNoZeros) - pastRates - gammaln(meas + 1))) ;


last100 = zeros(1,100);

for it = 1:max_iter

    % add an adaptive component
    if (it>adapt_after) && (mod(it,100)==1)
        if sum(last100)<20
            disp('Reducing proposal variance to increase acceptance rate.')
            std_prop_rb = stdMult*std_prop_rb;
            std_prop_ab = stdMult*std_prop_ab;

        elseif sum(last100)>40
            disp('Increasing proposal variance to decrease acceptance rate.')
            std_prop_rb = 1/stdMult*std_prop_rb;
            std_prop_ab = 1/stdMult*std_prop_ab;

        end
    end
    
    % generate candidate
    rbprop = std_prop_rb.*randn(nTarget,1) + rbpast;
    abprop = std_prop_ab.*randn(nTarget,1) + abpast;
    
    % evaluate prior first - we don't want to waste time evaluating
    % the rest of the likelihood if this is zero
    g_theta = (prod(rbprop<rbmax)&&prod(rbprop>rbmin)...
        &&prod(abprop<abmax)&&prod(abprop>abmin)); % indicator function
    

    [proposedObject,rejectFlag] = makeProposedObjectForegroundBackgroundTilt_multiTarget_FOV(phiPast, hprop, ...
    rprop, rbprop, bpast);

    if (g_theta == 0) || rejectFlag % reject
        n_throw = n_throw + 1;
        n_total = n_total + 1;
        
        % the last entry of last100 is the most recent
        last100(1:99) = last100(2:100); % get rid of oldest entry, i.e. the first one
        last100(end) = 0; % update last entry with most recent result
        
        acceptance_ratio = n_accept/n_total;
        % do some plotting
        set(0,'CurrentFigure',fgResult)

        subplot(121)
        title({['Background albedo draw'];['${n_{\rm accept}} = $' num2str(sum(last100))]},'Interpreter','latex','FontSize',16)
        grid on
        
        subplot(122)
        title({['Background range draw'];['$n_{\rm total} = $' num2str(n_total)]},'Interpreter','latex','FontSize',16)
        grid on
        

        
        
        drawnow
    else
        
        proposedForegroundRates = zeros(nPix,num_time_bins);
        proposedOccludedRates = zeros(nPix,num_time_bins);
        
        
            if twoBackgroundRegionFlag == 1
                for t = 1:nTarget
                    % pull out object
                    thisProposedTarget = proposedObject((t-1)*3+1:t*3,:);
                    thisProposedTargetMerged = mergeRegions(thisProposedTarget);
                    [~,~,thisForegroundRates,thisOccludedRates] = ...
                        computeForwardModel_adaptiveSegmentSize_par(thisProposedTargetMerged,params,plotFlag);
                    proposedForegroundRates = proposedForegroundRates + apast(t)*thisForegroundRates;
                    proposedOccludedRates = proposedOccludedRates + abprop(t)*thisOccludedRates;
                end
            else
                proposedObjectStore ={};
                for t = 1:nTarget
                    thisProposedTarget = proposedObject((t-1)*3+1,:);
                    thisProposedTarget(2,:) = proposedObject(t*3,:);
                    proposedObjectStore((t-1)*2+1:2*t,:) = thisProposedTarget;
                    [~,~,thisForegroundRates,thisOccludedRates] = ...
                        computeForwardModel_adaptiveSegmentSize_par(thisProposedTarget,params,plotFlag);
                    proposedForegroundRates = proposedForegroundRates + apast(t)*thisForegroundRates;
                    proposedOccludedRates = proposedOccludedRates + abprop(t)*thisOccludedRates;
                end
                proposedObject = proposedObjectStore;
            end
            
        

        proposedForegroundRates = circshift(proposedForegroundRates,nshift,2);
        proposedOccludedRates = circshift(proposedOccludedRates,nshift,2);
        
        proposedForegroundRates(indices2D,:) = [];
        proposedOccludedRates(indices2D,:) = [];
        
        proposedRates = laser_intensity*(proposedForegroundRates-proposedOccludedRates)...
            +spast*backgroundRate;


        proposedRatesNoZeros = proposedRates;
        proposedRatesNoZeros(proposedRatesNoZeros<1) = 1;
        
        
        
        if ~isempty(find(proposedRates<0))
            disp('Rejecting due to negative rates.')
            n_throw = n_throw + 1;
            n_total = n_total + 1;
            
            % the last entry of last100 is the most recent
            last100(1:99) = last100(2:100); % get rid of oldest entry, i.e. the first one
            last100(end) = 0; % update last entry with most recent result
            
            acceptance_ratio = n_accept/n_total;
            % do some plotting
            set(0,'CurrentFigure',fgResult)
            
            subplot(121)
            title({['Background albedo draw'];['${n_{\rm accept}} = $' num2str(sum(last100))]},'Interpreter','latex','FontSize',16)
            grid on
            
            subplot(122)
            title({['Background range draw'];['$n_{\rm total} = $' num2str(n_total)]},'Interpreter','latex','FontSize',16)
            grid on
            
            
            drawnow
        else
            % evaluate likelihood
            log_f_x_given_theta =   sum(sum(meas.*log(proposedRatesNoZeros) - proposedRates - gammaln(meas + 1))) ;
            
            
            % compute acceptance ratio
            alpha = min(g_theta*exp( log_f_x_given_theta -  log_f_x_given_theta_past),1);
            
            
            u = rand; % uniform random number
            
            if u<=alpha % accept
                n_accept = n_accept + 1;
                n_total = n_total + 1;
                acceptance_ratio = n_accept/n_total;
                
                % update more recent count
                last100(1:99) = last100(2:100); % get rid of oldest entry, i.e. the first one
                last100(end) = 1; % update last entry with most recent result
                
                % update values
                rbpast = rbprop;
                rbhist = [rbhist, rbprop];
                abpast = abprop;
                abhist = [abhist, abprop];

                
                log_f_x_given_theta_past = log_f_x_given_theta;
                
                modeEst.ab_mode = zeros(nTarget,1);
                modeEst.rb_mode = zeros(nTarget,1);
                modeEst.b_mode = zeros(nTarget,1);

                figure(333)
                count = 1;
                for t = 1:nTarget
                    subplot(nTarget,2,count)
                    histgm = histogram(abhist(t,:));
                    i_max = max(find(histgm.Values==max(histgm.Values)));
                    modeEst.ab_mode(t) = mean(histgm.BinEdges(i_max:i_max+1));
                    xlabel(['albedo, mode = ' num2str(modeEst.ab_mode(t))],'Interpreter','latex','FontSize',14)
                    count = count + 1;
                    
                    subplot(nTarget,2,count)
                    histgm = histogram(rbhist(t,:));
                    i_max = max(find(histgm.Values==max(histgm.Values)));
                    modeEst.rb_mode(t) = mean(histgm.BinEdges(i_max:i_max+1));
                    xlabel(['range, mode = ' num2str(modeEst.rb_mode(t))], 'Interpreter','latex','FontSize',14)                    
                    count = count + 1;

                end
                

                
                proposedObjectMode = makeProposedObjectForegroundBackgroundTilt_multiTarget_FOV(phiPast, hprop, ...
                    rprop, modeEst.rb_mode, modeEst.b_mode);

                set(0,'CurrentFigure',fgProgress)
                subplot(122)
                plotObjects(proposedObjectMode,'blue')
                if ~isempty(objects)
                    plotObjects(objects,'red')
                end
                scatter3(l_xyz0(1),l_xyz0(2),l_xyz0(3),'r','filled')
                hold off
                
                set(0,'CurrentFigure',fgResult)
                subplot(121)
                hold off
                plot(abhist(1,:),'b--o')
%                 plot(abhist','b--o')
                hold on
                if nTarget>1
                    plot(abhist(2,:),'r--o')
                end
                title({['Background albedo draw'];['${n_{\rm accept}} = $' num2str(sum(last100))]},'Interpreter','latex','FontSize',16)
                xlabel('index of accepted iterations','Interpreter','latex','FontSize',14)
                ylabel('a','Interpreter','latex','FontSize',14)
                grid on
                
                subplot(122)
                hold off
                plot(rbhist(1,:),'b--o')
%                 plot(rbhist','b--o')
                hold on
                if nTarget>1
                    plot(rbhist(2,:),'r--o')
                end
                title({['Background range draw'];['$n_{\rm total} = $' num2str(n_total)]},'Interpreter','latex','FontSize',16)
                %             title({['Background range draw']},'Interpreter','latex','FontSize',16)
                xlabel('index of accepted iterations','Interpreter','latex','FontSize',14)
                ylabel('r','Interpreter','latex','FontSize',14)
                grid on

                set(0,'CurrentFigure',fgProgress)
                subplot(121)
                plotObjects(proposedObject,'blue')
                if ~isempty(objects)
                    plotObjects(objects,'red')
                end
                scatter3(l_xyz0(1),l_xyz0(2),l_xyz0(3),'r','filled')            
                hold off
                
                drawnow
            else % reject
                n_throw = n_throw + 1;
                n_total = n_total + 1;
                acceptance_ratio = n_accept/n_total;
                
                last100(1:99) = last100(2:100); % get rid of oldest entry, i.e. the first one
                last100(end) = 0; % update last entry with most recent result
                
                % do some plotting
                set(0,'CurrentFigure',fgResult)
                subplot(121)
                title({['Background albedo draw'];['${n_{\rm accept}} = $' num2str(sum(last100))]},'Interpreter','latex','FontSize',16)
                grid on
                
                subplot(122)
                title({['Background range draw'];['$n_{\rm total} = $' num2str(n_total)]},'Interpreter','latex','FontSize',16)
                grid on

                drawnow
            end
        end
        
    end
    
end
set(0,'CurrentFigure',fgResult)
% set(gcf,'Position',[10 10 50 30])
set(gcf,'PaperUnits','centimeters');
set(gcf,'PaperSize',4*[7 7]);
fig = gcf;
fig.PaperUnits = 'centimeters';
fig.PaperPosition = 4*[0 0 7 7];
fig.Units = 'centimeters';
fig.PaperSize=4*[7 7];
fig.Units = 'centimeters';
print(fig,[dirpath 'ProposalHistoryBackground_' num2str(posNum)],'-dpdf','-r200');

figure(333)
% set(gcf,'Position',[10 10 50 30])
set(gcf,'PaperUnits','centimeters');
set(gcf,'PaperSize',4*[7 7]);
fig = gcf;
fig.PaperUnits = 'centimeters';
fig.PaperPosition = 4*[0 0 7 7];
fig.Units = 'centimeters';
fig.PaperSize=4*[7 7];
fig.Units = 'centimeters';
print(fig,[dirpath 'McmcHistogramBackground_' num2str(posNum)],'-dpdf','-r200');


hist.abhist = abhist;
hist.rbhist = rbhist;


end

