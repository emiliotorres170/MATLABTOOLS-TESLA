function [u1t, u2t, u3t, ke_mean, simulationTime] = npy2mat_2(N, Nt, np, loc)


%#ok<*NOPTS>;
close all; 
clc; 
format compact 
warning('off');
here = pwd; 
chdir(loc);
st = 1; 
et = Nt;


u1t = zeros([N,N,N,et+1]);
u2t = zeros([N,N,N,et+1]);
u3t = zeros([N,N,N,et+1]);

u1 = zeros([N,N,N]); 
u2 = zeros([N,N,N]); 
u3 = zeros([N,N,N]);
simulationTime = zeros(et-st+1,1); 
for i = st:et+1 
    ii = i-st+1;
    
    if i<11       
        id =    ['00' num2str(i-1)];
    elseif i<101  
        id =    ['0' num2str(i-1)];
    elseif i<1001
        id =    num2str(i-1);
    elseif i<10001 
        id =    num2str(i-1);
    end
    
    simulationTime(ii) = readNPY(['SimulationTime_' id '.npy']);
    
    for j = 1:np
        
        if j<11       
            proc=['00' num2str(j-1)];
        elseif j<101  
            proc=['0' num2str(j-1)];
        elseif j<1001 
            id=num2str(j-1);
        end
        
        fid1 = ['Velocity1_' id '_' proc '.npy'];
        fid2 = ['Velocity2_' id '_' proc '.npy'];
        fid3 = ['Velocity3_' id '_' proc '.npy'];
        tu1(j,:,:,:) = readNPY(fid1);
        tu2(j,:,:,:) = readNPY(fid2);
        tu3(j,:,:,:) = readNPY(fid3);
    end
    
    for j = 1:np
        u1((j-1)*N/np+1:(N/np)*j,:,:) = squeeze(tu1(j,:,:,:));
        u2((j-1)*N/np+1:(N/np)*j,:,:) = squeeze(tu2(j,:,:,:));
        u3((j-1)*N/np+1:(N/np)*j,:,:) = squeeze(tu3(j,:,:,:));
    end

    u1t(:,:,:,i)    = u1;
    u2t(:,:,:,i)    = u2;
    u3t(:,:,:,i)    = u3;
end

dim     = size(u1t);
ke_mean = zeros(1, dim(end));
for i = 1:dim(end)
    ke_mean(i) = 0.5*mean(mean(mean(u1t(:,:,:,i).^2 + u2t(:,:,:,i).^2 + u3t(:,:,:,i).^2)));
end
chdir(here);

end 
