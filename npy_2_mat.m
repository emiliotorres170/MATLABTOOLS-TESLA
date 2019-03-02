function [ux, uy, uz, time] = npy_2_mat(N, Nt, np, Loc)
%=========================================================================%
% Purpose:                                                                %
%   The purpose of thsi scriopt to read the .npy files and store the data %
%   into .mat                                                             %
%                                                                         %
%   Inputs:                                                               %
%       1. N (integer) number of spatial steps                            %
%       2. Nt (integer) number of temporal steps                          %
%       3. np (integer) number of proccessors                             %
%       4. time                                                           %
%                                                                         %
%   Outputs:                                                              %
%       1. u_x                                                            %
%       2. u_y                                                            %
%       3. u_y                                                            %
%       4. time                                                           %
%                                                                         %
% Author:                                                                 %
%   Emilio Torres                                                         %
%=========================================================================%
here    = pwd;
ux      = zeros([N, N, N, Nt+1]);
uy      = zeros(N, N, N, Nt+1);
uz      = zeros(N, N, N, Nt+1);
time    = zeros(1, Nt+1);
cd(Loc);
for i = 0:Nt
    if i < 10
        vel_id  = ['00' num2str(i)];  
    elseif i > 99
        vel_id  = num2str(i);
    elseif i >= 10
        vel_id  = ['0' num2str(i)];
    end
    
    for k = 0:np-1
        if k < 10
            np_id  = ['00' num2str(k)];  
        elseif k >= 10
            np_id  = ['0' num2str(k)];
        elseif k > 99
            np_id  = num2str(k);
        end
        u1temp(k+1, :, :, :) = readNPY(['Velocity1_' vel_id '_' np_id '.npy']);
        u2temp(k+1, :, :, :) = readNPY(['Velocity2_' vel_id '_' np_id '.npy']);
        u3temp(k+1, :, :, :) = readNPY(['Velocity3_' vel_id '_' np_id '.npy']);
      end 
    for j = 1:np
        ux((j-1)*N/np+1:(N/np)*j,:,:,i+1) = squeeze(u1temp(j,:,:,:));
        uy((j-1)*N/np+1:(N/np)*j,:,:,i+1) = squeeze(u1temp(j,:,:,:));
        uz((j-1)*N/np+1:(N/np)*j,:,:,i+1) = squeeze(u1temp(j,:,:,:));
        time(i+1) = readNPY(['SimulationTime_' vel_id '.npy']);
    end
end
cd(here);
end 

