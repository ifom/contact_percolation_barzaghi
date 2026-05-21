% This script calculates the CM and RMS velocities for each sample at each time point and averages them over the time window
% Run after DYNAMIX_01_BuildTrajectories.m
% C. Guidolin

clear
close all
clc
tic


DataSet = 2021; % The dataset is mostly identified either by the year (YYYY) or the full date (YYYYMMDD) in case of more than one experiment per year
% Read the function GetDatasetInfo for special cases or if in doubt

HD = 'D'; % Hard disk

save_results = false;


SmoothingWindow = 21;

% It loads:
filename_Traj_C = '\MAT_Traj_C.mat';
filename_Traj_R = '\MAT_Traj_R.mat';
% It saves:
filename_SAVE = '\MAT_Vrms_Vcm_psi.mat';

% It uses the functions:
% GetDatasetInfo.m; GetVelocities.m; GetVelocities_C_R.m; getVrms.m
addpath([HD ':\GUIDOLIN_CHIARA\MIXED_POPULATIONS\LATEST_CODE_VERSION\Functions\'])



% Load dataset info
[main_folder, sample_name_cell, fov, scale, dt, img_size, dist_edge, n_frames] = GetDatasetInfo(DataSet,HD);

n_samples = numel(sample_name_cell);

% Loop over the samples
for sss = 1 : n_samples

    sample_name = sample_name_cell{sss};

    disp(['Sample: ' sample_name])

    % Loop over the fields of view
    for FOV = fov

        disp(['FOV = ' num2str(FOV)])

        % Folder containing the results
        folder = [main_folder sample_name '\FOV_' num2str(FOV)];

        switch sample_name

            case {'C100','C100_RACDN'} % Only CTRL

                load([folder filename_Traj_C],'cell_ALL_trajectories_C','tab_C_spots') % it loads the 'new' tab_C_spots with only spots corresponding to selected trajectories

%                 cell_ALL_trajectories_C = SelectTrajectories(cell_ALL_trajectories_C,n_frames/4,img_size,dist_edge);


                [CM_C, V_CM_C, V_CM_smooth_C, V_rms_CM_C, V_rms_lab_C, V_rms_lab_CORR_C, V_fake_C] = GetVelocities(n_frames,tab_C_spots,SmoothingWindow);

                dCM = [[0 0]; cumsum(V_CM_C)]; % CM displacement

                V_CM_mod_C = sqrt(V_CM_C(:,1).^2 + V_CM_C(:,2).^2);
                V_CM_mod_smooth_C = sqrt(V_CM_smooth_C(:,1).^2 + V_CM_smooth_C(:,2).^2);

                psi_C = V_CM_mod_C ./ V_rms_lab_C; % this is not used
                psi_CORR_C = V_CM_mod_smooth_C ./ V_rms_lab_CORR_C; % this is used afterwards

                % Average over time window
                V_CM_mod_AVG_C = mean(V_CM_mod_C,'omitnan');
                V_CM_mod_AVG_smooth_C = mean(V_CM_mod_smooth_C,'omitnan');
                V_rms_AVG_C = sqrt( mean(V_rms_CM_C,'omitnan')); % the sqrt is taken here after averaging over the time window
                V_rms_lab_CORR_AVG_C = mean(V_rms_lab_CORR_C,'omitnan'); % the sqrt in this case was already taken in the function GetVelocities.m
                psi_AVG_C = mean(psi_C,'omitnan');
                psi_CORR_AVG_C = mean(psi_CORR_C,'omitnan');

                V_rms_lab_CORR_AVG = V_rms_lab_CORR_AVG_C;
                V_CM_mod_AVG = V_CM_mod_AVG_C;
                V_CM_mod_AVG_smooth = V_CM_mod_AVG_smooth_C;
                V_rms_AVG = V_rms_AVG_C;
                psi_AVG = psi_AVG_C;
                psi_CORR_AVG = psi_CORR_AVG_C;


            case 'R100' % Only RAB5

                load([folder filename_Traj_R],'cell_ALL_trajectories_R','tab_R_spots')

%                 cell_ALL_trajectories_R = SelectTrajectories(cell_ALL_trajectories_R,n_frames/4,img_size,dist_edge);

                [CM_R, V_CM_R, V_CM_smooth_R, V_rms_CM_R, V_rms_lab_R, V_rms_lab_CORR_R, V_fake_R] = GetVelocities(n_frames,tab_R_spots,SmoothingWindow);

                dCM = [[0 0]; cumsum(V_CM_R)]; % CM displacement

                V_CM_mod_R = sqrt(V_CM_R(:,1).^2 + V_CM_R(:,2).^2);
                V_CM_mod_smooth_R = sqrt(V_CM_smooth_R(:,1).^2 + V_CM_smooth_R(:,2).^2);

                psi_R = V_CM_mod_R ./ V_rms_lab_R;
                psi_CORR_R = V_CM_mod_smooth_R ./ V_rms_lab_CORR_R;

                % Average over time window
                V_CM_mod_AVG_R = mean(V_CM_mod_R,'omitnan');
                V_CM_mod_AVG_smooth_R = mean(V_CM_mod_smooth_R,'omitnan');
                V_rms_AVG_R = sqrt(mean(V_rms_CM_R,'omitnan')); % the sqrt is taken here after averaging over the time window
                V_rms_lab_CORR_AVG_R = mean(V_rms_lab_CORR_R,'omitnan'); % the sqrt in this case was already taken in the function GetVelocities.m
                psi_AVG_R = mean(psi_R,'omitnan');
                psi_CORR_AVG_R = mean(psi_CORR_R,'omitnan');


                V_rms_lab_CORR_AVG = V_rms_lab_CORR_AVG_R;
                V_CM_mod_AVG = V_CM_mod_AVG_R;
                V_CM_mod_AVG_smooth = V_CM_mod_AVG_smooth_R;
                V_rms_AVG = V_rms_AVG_R;
                psi_AVG = psi_AVG_R;
                psi_CORR_AVG = psi_CORR_AVG_R;


            otherwise    % Mixed samples

                load([folder filename_Traj_C],'cell_ALL_trajectories_C','tab_C_spots')
                load([folder filename_Traj_R],'cell_ALL_trajectories_R','tab_R_spots')

%                 cell_ALL_trajectories_C = SelectTrajectories(cell_ALL_trajectories_C,n_frames/2,img_size,dist_edge);
%                 cell_ALL_trajectories_R = SelectTrajectories(cell_ALL_trajectories_R,n_frames/2,img_size,dist_edge);


                % ALL cells together

                tab_ALL_spots = [tab_C_spots; tab_R_spots];

                [CM, V_CM, V_CM_smooth, V_rms_CM, V_rms_lab, V_rms_lab_CORR, V_fake] = GetVelocities(n_frames,tab_ALL_spots,SmoothingWindow);

                [VRMS_ALL_TRAJ,vCM] = getVrms([cell_ALL_trajectories_C; cell_ALL_trajectories_R],n_frames,img_size,dist_edge);

                dCM = [[0 0]; cumsum(V_CM)]; % CM displacement

                V_CM_mod = sqrt(V_CM(:,1).^2 + V_CM(:,2).^2);
                V_CM_mod_smooth = sqrt(V_CM_smooth(:,1).^2 + V_CM_smooth(:,2).^2);

                psi = V_CM_mod ./ V_rms_lab;
                psi_CORR = V_CM_mod_smooth ./ V_rms_lab_CORR;

                % Average over time window
                V_CM_mod_AVG = mean(V_CM_mod,'omitnan');
                V_CM_mod_AVG_smooth = mean(V_CM_mod_smooth,'omitnan');
                V_rms_AVG = sqrt( mean(V_rms_CM,'omitnan')); % the sqrt is taken here after averaging over the time window
                V_rms_lab_CORR_AVG = mean(V_rms_lab_CORR,'omitnan'); % the sqrt in this case was already taken in the function GetVelocities
                psi_AVG = mean(psi,'omitnan');
                psi_CORR_AVG = mean(psi_CORR,'omitnan');


                % only CTRL

                [CM_C, V_CM_C, V_CM_smooth_C, V_rms_CM_C, V_rms_CM_all_C, V_rms_lab_C, V_rms_lab_CORR_C, V_fake_C] = GetVelocities_C_R(n_frames,tab_C_spots,SmoothingWindow,V_CM);

                [VRMS_C_TRAJ,vCM_C] = getVrms(cell_ALL_trajectories_C,n_frames,img_size,dist_edge);


                V_CM_mod_C = sqrt(V_CM_C(:,1).^2 + V_CM_C(:,2).^2);
                V_CM_mod_smooth_C = sqrt(V_CM_smooth_C(:,1).^2 + V_CM_smooth_C(:,2).^2);

                psi_C = V_CM_mod_C ./ V_rms_lab_C;
                psi_CORR_C = V_CM_mod_smooth_C ./ V_rms_lab_CORR_C;

                % Average over time window
                V_CM_mod_AVG_C = mean(V_CM_mod_C,'omitnan');
                V_CM_mod_AVG_smooth_C = mean(V_CM_mod_smooth_C,'omitnan');
                V_rms_AVG_C = sqrt( mean(V_rms_CM_all_C,'omitnan') ); % in the CM of the whole monolayer
                V_rms_lab_CORR_AVG_C = mean(V_rms_lab_CORR_C,'omitnan');
                psi_AVG_C = mean(psi_C,'omitnan');
                psi_CORR_AVG_C = mean(psi_CORR_C,'omitnan');


                % only RAB5

                [CM_R, V_CM_R, V_CM_smooth_R, V_rms_CM_R, V_rms_CM_all_R, V_rms_lab_R, V_rms_lab_CORR_R, V_fake_R] = GetVelocities_C_R(n_frames,tab_R_spots,SmoothingWindow,V_CM);

                [VRMS_R_TRAJ,vCM_R] = getVrms(cell_ALL_trajectories_R,n_frames,img_size,dist_edge);

                V_CM_mod_R = sqrt(V_CM_R(:,1).^2 + V_CM_R(:,2).^2);
                V_CM_mod_smooth_R = sqrt(V_CM_smooth_R(:,1).^2 + V_CM_smooth_R(:,2).^2);

                psi_R = V_CM_mod_R ./ V_rms_lab_R;
                psi_CORR_R = V_CM_mod_smooth_R ./ V_rms_lab_CORR_R;

                % Average over time window
                V_CM_mod_AVG_R = mean(V_CM_mod_R,'omitnan');
                V_CM_mod_AVG_smooth_R = mean(V_CM_mod_smooth_R,'omitnan');
                V_rms_AVG_R = sqrt( mean(V_rms_CM_all_R,'omitnan') ); % in the CM of the whole monolayer
                V_rms_lab_CORR_AVG_R = mean(V_rms_lab_CORR_R,'omitnan');
                psi_AVG_R = mean(psi_R,'omitnan');
                psi_CORR_AVG_R = mean(psi_CORR_R,'omitnan');



        end

        % Save results
        if save_results == true
            disp('Saving results ..')
            save([folder filename_SAVE])
        end

        toc

    end
    toc

end
toc
