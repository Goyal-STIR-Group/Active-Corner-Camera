
% plot some histograms
leg = [1,512,nPix];
for ix = leg
    
    figure(200)
    hold on
    plot(bin_size_m*(1:num_time_bins),meas(ix,:))
    xlabel('Range [m]','Interpreter','latex','FontSize',16)
    ylabel('Photon counts','Interpreter','latex','FontSize',16)
    grid on
    axis square
    
    
    figure(203)
    hold on
    plot(bin_size_m*(1:num_time_bins),backgroundRate(ix,:))
    xlabel('Range [m]','Interpreter','latex','FontSize',16)
    ylabel('Photon rate','Interpreter','latex','FontSize',16)
    grid on
    axis square
    
end

legendStrings = "Pixel " + string(leg);

figure(200)
title('New Frame Measurements','FontSize',14,'Interpreter','latex')
legend(legendStrings,'FontSize',12,'Interpreter','latex')
fig = gcf;
fig.PaperUnits = 'centimeters';
fig.PaperPosition = 2*[-.3 0 5 4];
fig.Units = 'centimeters';
fig.PaperSize=2*[4 4];
print(fig,[dirpath 'measuredHistograms_' num2str(posNum)],'-dpdf','-r200');

figure(203)
title('Stationary Scene Rates','FontSize',14,'Interpreter','latex')
legend(legendStrings,'FontSize',12,'Interpreter','latex')
fig = gcf;
fig.PaperUnits = 'centimeters';
fig.PaperPosition = 2*[-.3 0 5 4];
fig.Units = 'centimeters';
fig.PaperSize=2*[4 4];
print(fig,[dirpath 'StationarySceneRate_' num2str(posNum)],'-dpdf','-r200');
legend(legendStrings,'FontSize',12,'Interpreter','latex')


figure
plot(bin_size_m*(1:num_time_bins),mean(meas(indices2D,:),1))
hold on
grid on
plot(bin_size_m*(1:num_time_bins),mean(backgroundRate(indices2D,:),1))
title('Spatially Averaged Measurements','FontSize',14,'Interpreter','latex')
xlabel('Range [m]','FontSize',16,'Interpreter','latex')
ylabel('Counts','FontSize',16,'Interpreter','latex')
xlim([0 11])
ylim([0 13])
legend('Stationary scene','Object in motion','Interpreter','latex','FontSize',13,'NumColumns',2)


set(gcf,'PaperUnits','centimeters');
set(gcf,'PaperSize',1.5*[10 4]);
fig = gcf;
fig.PaperUnits = 'centimeters';
fig.PaperPosition = 1.5*[0 0 10 4];
fig.Units = 'centimeters';
fig.PaperSize=1.5*[10 4];
fig.Units = 'centimeters';
print(fig,[dirpath 'spatialSumBackgroundCompare_' num2str(posNum)],'-dpdf','-r200');


figure
plot(bin_size_m*(1:num_time_bins),...
    mean(meas_dif_scale(indices2D,:),1))
grid on
xlabel('Range [m]','FontSize',16,'Interpreter','latex')
ylabel('Counts','FontSize',16,'Interpreter','latex')
xlim([0 11])
ylim([-3 9.5])
title('Spatially Averaged Difference Meas.','FontSize',14,'Interpreter','latex')



set(gcf,'PaperUnits','centimeters');
set(gcf,'PaperSize',1.5*[10 4]);
fig = gcf;
fig.PaperUnits = 'centimeters';
fig.PaperPosition = 1.5*[0 0 10 4];
fig.Units = 'centimeters';
fig.PaperSize=1.5*[10 4];
fig.Units = 'centimeters';
print(fig,[dirpath 'difference_' num2str(posNum)],'-dpdf','-r200');

