close all; clc; clear all; format compact %#ok<*NOPTS>;

st = 1; et = 3;

fdir = 'stats'
N = 64; 
np = 4;
visc = 0.000185;

nx = [N N N];
nxc = [N/2+1 N/2+1 N/2+1];
k = 1:N/2+1;

u1 = zeros([N,N,N]); u2 = zeros([N,N,N]); u3 = zeros([N,N,N]);
simulationTime = zeros(et-st+1,1); timeStep = zeros(et-st+1,1);


%% Load Data and Record Statistics Over Time
time = (st:et)*2/3;
for i = st:et+1 
    % Load Data -----------------------------------------------------------
    disp(['step = ' num2str((i))])
    ii = i-st+1;
    
    if (i<11);       id=['00' num2str(i-1)];
    elseif (i<101);  id=['0' num2str(i-1)];
    elseif (i<1001); id=num2str(i-1);
    elseif (i<10001); id=num2str(i-1);
    end
    
    simulationTime(ii) = readNPY(['SimulationTime_' id '.npy']);
%    timeStep(ii) = readNPY(['TimeStep_' id '.npy']);

    for j = 1:np
        if (j<11);       proc=['00' num2str(j-1)];
        elseif (j<101);  proc=['0' num2str(j-1)];
        elseif (j<1001); id=num2str(j-1);
        end
        
        fid1 = ['Velocity1_' id '_' proc '.npy'];
        
        fid2 = ['Velocity2_' id '_' proc '.npy'];
        fid3 = ['Velocity3_' id '_' proc '.npy'];
      
       
        tu1(j,:,:,:) = readNPY(fid1);
        tu2(j,:,:,:) = readNPY(fid2);
%         tu3(j,:,:,:) = readNPY(fid3);


    end
    
    for j = 1:np
        var = size(squeeze(tu1(j,:,:,:)));
        u1((j-1)*N/np+1:(N/np)*j,:,:) = squeeze(tu1(j,:,:,:));
        u2((j-1)*N/np+1:(N/np)*j,:,:) = squeeze(tu2(j,:,:,:));
%         u3((j-1)*N/np+1:(N/np)*j,:,:) = squeeze(tu3(j,:,:,:));

    end
    
    [duxdx, duxdy, duxdz] = gradient(u1);
    [duydx, duydy, duydz] = gradient(u2);
    [duzdx, duzdy, duzdz] = gradient(u3);

    % Enstrophy
    w1 = 2*(duzdy-duydz); 
    w2 = 2*(duxdz-duzdx);
    w3 = 2*(duydx-duxdy);
    enstrophy = 0.5*(w1.^2 + w2.^2 + w3.^2);
    
    % Reynolds Stresses
    uxux(ii) = mean(mean(mean(u1.^2)));
    uyuy(ii) = mean(mean(mean(u2.^2)));
    uzuz(ii) = mean(mean(mean(u3.^2)));    
    uxuy(ii) = mean(mean(mean(u1.*u2)));
    uyuz(ii) = mean(mean(mean(u2.*u3)));
    uzux(ii) = mean(mean(mean(u3.*u1)));
    ke(ii) = 0.5*mean(mean(mean((u1.^2)+(u2.^2)+(u3.^2))));
    kefield = 0.5*((u1.^2)+(u2.^2)+(u3.^2));
    
    % Reynold's Stresses Variances
    uxux2(ii) = mean(mean(mean(u1.^4)));
    uyuy2(ii) = mean(mean(mean(u2.^4)));
    uzuz2(ii) = mean(mean(mean(u3.^4)));    
    uxuy2(ii) = mean(mean(mean(u1.*u2.*u1.*u2)));
    uyuz2(ii) = mean(mean(mean(u2.*u3.*u2.*u3)));
    uzux2(ii) = mean(mean(mean(u3.*u1.*u3.*u1)));
    
    % Turbulent Dissipation Rate
    s12 = 0.5*(duxdy+duydx);  s33 = 0.5*(duzdz+duzdz);
    s11 = 0.5*(duxdx+duxdx);  s22 = 0.5*(duydy+duydy);
    s23 = 0.5*(duydz+duzdy);  s13 = 0.5*(duzdx+duxdz);
    S = s11.*s11 + s22.*s22 + s33.*s33 + 2*s12.*s12 + 2*s23.*s23 + 2*s13.*s13;
    Sbar(ii) = mean(mean(mean(S)));
    epsilon(ii) = 2*visc*mean(mean(mean(S))); 
    
    % Subgrid Production
    % Positive is backscatter, negative is forward scatter
%     P(:,:,:) = tau11.*s11 + tau22.*s22 + tau33.*s33 + 2*tau12.*s12 + 2*tau23.*s23 + 2*tau13.*s13;
%     Pp = P.*(P>0);
%     Pm = P.*(P<0);
%     sigmaP = sqrt(mean(mean(mean(P.^2))));
%     sigmaPp = sqrt(mean(mean(mean(Pp.^2))));
%     sigmaPm = sqrt(mean(mean(mean(Pm.^2))));
%     
%     sssP(ii) = trapz(trapz(trapz(P)))/8/pi/pi/pi;
%     Ppbar(ii) = mean(mean(mean(P(P>0))));
%     Pmbar(ii) = mean(mean(mean(P(P<0))));
%     Pbar(ii) = mean(mean(mean(P)));
%     
%     sigP(ii) = sqrt(mean(mean(mean(P.^2))));
%     maxP(ii) = max(max(max(P)));
    
    
%     hold on
%     contourf(squeeze(u1(:,:,33)),20); contour(squeeze(u1(:,:,33)),20)
%     colorbar; colormap jet; %caxis([0 150])
%     title(['x-direction velocity field at ' num2str(i)],'fontsize',16);  
%     set(gca,'fontsize',14); set(gcf,'Position', [1000, 300, 800, 645])
%     shg
%     drawnow
    
    
%     hold on
%     contourf(squeeze(Pp(:,:,33)),20); contour(squeeze(Pp(:,:,33)),20)
%     colorbar; colormap jet; %caxis([0 150])
%     title(['Production field at ' num2str(simulationTime(i))],'fontsize',16);  
%     set(gca,'fontsize',14); set(gcf,'Position', [1000, 300, 800, 645])
%     shg
%     drawnow
    
%     clf
%     sigmaP = sqrt(mean(mean(mean(P.^2))));
%     Pt = P.*(P>(3*sigmaP));
%     x = 1:N; y = 1:N; z = 1:N;
%     [X,Y,Z] = meshgrid(x,y,z);
%     xslice = 1:1:32; yslice = 1:1:32; zslice = 1:1:32;
%     contourslice(X,Y,Z,Pt,xslice,yslice,zslice)
%     title(['P+ field at ' num2str(i)],'fontsize',16); 
%     axis([0 32 0 32 0 32]); %caxis([15 150])
%     colorbar; colormap jet; set(gcf,'Position', [900, 100, 1000, 800]); view([45 55 45]); grid on
%     shg
%     drawnow
    
    
    % RMS Velocity
%     vrms(ii) = sqrt(mean(mean(mean(u1.*u1 + u2.*u2 + u3.*u3))));
    
end

