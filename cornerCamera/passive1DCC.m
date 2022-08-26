function [k,v] = passive1DCC(S, ST, B, invB, N1, penumbraIncident,...
    lam2,nIter,stepsize)

% This code runs the passive corner camera
% v is the 1D scene estimate as a function of azimuthal angle
% k is the ambient light estimate

cost = [];

N2 = size(S,2); % number of scene elements, i.e. number of angles
npix = size(S,1);
pix_dim = sqrt(npix);

vp = zeros(N2,1);
u = vp;

kp = 0;
m = 0;

tp = 1;

p = zeros(N1,1);
kHist = [];
v = zeros(N2,1);


for iter = 1:nIter
    temp = (m + S*u - penumbraIncident);
    v = u - stepsize*(ST*temp);

    [v,p] = constrained_ell1(v, lam2, B, invB, 100, [0,inf],p);


    k = m - stepsize*(ones(1,npix)*temp);
    v(v<0) = 0;
    k(k<0) = 0;
    
    
    t = (sqrt(1+4*tp^2)+1)/2;
    
    u =  v + (tp-1)/t*(v-vp);
    m =  k + (tp-1)/t*(k-kp);
%     
%     if norm(v-vp)<1e-14
%         break
%     end
    tp = t;
    vp = v;
    kp = k;
    
    kHist = [kHist, k];
    
    
    if mod(iter,200)== 0
        cost(iter/200) = norm(penumbraIncident - S*v - k,2)^2 +lam2*norm(B(v),1);
        figure(33)
        subplot(231)
        plot(cost)
        grid on
        axis square
        title('Cost','Interpreter','latex','FontSize',14)
        xlabel('Iteration','Interpreter','latex','FontSize',12)

        subplot(232)
        plot(linspace(0,180,length(v)),v)
        grid on
        xlabel('Angle [deg]','Interpreter','latex','FontSize',12)
        title('1D Scene Estimate','Interpreter','latex','FontSize',14)          
        axis square
        
        subplot(233)

        plot(kHist)
        grid on
        title('Estimated bias','Interpreter','latex','FontSize',14)
        axis square
        
        reconstructedData = S*v+k;
        clim = [0 max(reconstructedData)];
        subplot(234)
        imagesc(reshape(reconstructedData,pix_dim,pix_dim),clim)
        title('Reconstructed Data','Interpreter','latex','FontSize',14)
        axis square
        axis off
        colorbar
        
        subplot(235)
        imagesc(reshape(penumbraIncident,pix_dim,pix_dim),clim)
        title('Measurement','Interpreter','latex','FontSize',14)
        axis square
        axis off
        colorbar
        
        subplot(236)
        imagesc(reshape(reconstructedData -penumbraIncident,pix_dim,pix_dim))
        title('Dif','Interpreter','latex','FontSize',14)
        axis square
        axis off
        colorbar
        
        drawnow
        
        
    end
end




end

