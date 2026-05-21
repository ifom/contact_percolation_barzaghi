% This script calculates the Mean Square Displacement for each sample and FOV at different time delays
% MSD is calculated from both raw and corrected trajectories and also in the CM sdr
% Run after DYNAMIX_04_CorrectTrajectories.m
% C. Guidolin

clear
close all
clc
tic



DataSet = 2021; % The dataset is mostly identified either by the year (YYYY) or the full date (YYYYMMDD) in case of more than one experiment per year
% Read the function GetDatasetInfo for special cases or if in doubt

HD = 'D'; % Hard disk

save_results = false;

% It uses the functions:
% GetDatasetInfo.m, SelectTrajectories.m, GetSqDisp.m, GetSqDisp_CORR.m, GetSqDisp_CM.m
addpath([HD ':\GUIDOLIN_CHIARA\MIXED_POPULATIONS\LATEST_CODE_VERSION\Functions\'])

% Load dataset info
[main_folder, sample_name_cell, fov, scale, ~, img_size, dist_edge, n_frames] = GetDatasetInfo(DataSet,HD);



% It loads:
filename_Traj_C = '\MAT_Traj_Corr_C.mat';
filename_Traj_R = '\MAT_Traj_Corr_R.mat';
filename_Traj_ALL = '\MAT_Traj_Corr_ALL.mat';
filename_LOAD_Velocities = '\MAT_Vrms_Vcm_psi.mat'; %'\MAT_Vrms_Vcm_psi_new.mat';
% It saves:
filename_SAVE = '\MAT_MSD_AllTraj.mat';


n_samples = numel(sample_name_cell);

min_length_traj = n_frames/2;


% Loop over the samples
for sss = 1 : n_samples

    sample_name = sample_name_cell{sss};

    disp(['Sample: ' sample_name])

    % Loop over the fields of view
    for FOV = fov

        disp(['FOV=' num2str(FOV)])

        % Folder containing trajectories files
        folder = [main_folder sample_name '\FOV_' num2str(FOV)];

        switch sample_name

            case {'C100','C100_RACDN'}

                load([folder filename_Traj_C],'cell_ALL_trajectories_C_Corr')
                load([folder filename_LOAD_Velocities],'dCM');
                % cell_ALL_trajectories_C_Corr = SelectTrajectories(cell_ALL_trajectories_C_Corr,min_length_traj);

            case 'R100'

                load([folder filename_Traj_R],'cell_ALL_trajectories_R_Corr')
                load([folder filename_LOAD_Velocities],'dCM');
                % cell_ALL_trajectories_R_Corr = SelectTrajectories(cell_ALL_trajectories_R_Corr,min_length_traj);

            otherwise

                load([folder filename_Traj_C],'cell_ALL_trajectories_C_Corr')
                load([folder filename_Traj_R],'cell_ALL_trajectories_R_Corr')
                load([folder filename_Traj_ALL],'cell_ALL_trajectories_Corr')
                load([folder filename_LOAD_Velocities],'dCM');

                % cell_ALL_trajectories_C_Corr = SelectTrajectories(cell_ALL_trajectories_C_Corr,min_length_traj);
                % cell_ALL_trajectories_R_Corr = SelectTrajectories(cell_ALL_trajectories_R_Corr,min_length_traj);
                % cell_ALL_trajectories_Corr = SelectTrajectories(cell_ALL_trajectories_Corr,min_length_traj);

        end


        % Preallocation
        MSD_C = zeros(n_frames,1);
        MSD_R = zeros(n_frames,1);
        MSD_ALL = zeros(n_frames,1);

        MSD_CORR_C = zeros(n_frames,1);
        MSD_CORR_R = zeros(n_frames,1);
        MSD_CORR_ALL = zeros(n_frames,1);

        MSD_CM_C = zeros(n_frames,1);
        MSD_CM_R = zeros(n_frames,1);
        MSD_CM_ALL = zeros(n_frames,1);

        % Select delays
        sel_delays = unique(round(logspace(0,log10(round(n_frames/2)),50)));

        % Loop over the selected delays
        for dt = sel_delays %1 : n_frames-1

            disp(dt)

            switch sample_name

                case {'C100','C100_RACDN'}

                    % Only CTRL cells
                    SqDisp_C = GetSqDisp(cell_ALL_trajectories_C_Corr,n_frames,dt);
                    SqDisp_CORR_C = GetSqDisp_CORR(cell_ALL_trajectories_C_Corr,n_frames,dt);
                    [SqDisp_CM_C,~] = GetSqDisp_CM(cell_ALL_trajectories_C_Corr,n_frames,dt,dCM);

                    MSD_C(dt) = mean(SqDisp_C,'all','omitnan');
                    MSD_CORR_C(dt) = mean(SqDisp_CORR_C,'all','omitnan');
                    MSD_CM_C(dt) = mean(SqDisp_CM_C,'all','omitnan');

                    MSD_ALL(dt) = MSD_C(dt);
                    MSD_CORR_ALL(dt) = MSD_CORR_C(dt);
                    MSD_CM_ALL(dt) = MSD_CM_C(dt);



                case 'R100'

                    % Only RAB5 cells
                    SqDisp_R = GetSqDisp(cell_ALL_trajectories_R_Corr,n_frames,dt);
                    SqDisp_CORR_R = GetSqDisp_CORR(cell_ALL_trajectories_R_Corr,n_frames,dt);
                    [SqDisp_CM_R,~] = GetSqDisp_CM(cell_ALL_trajectories_R_Corr,n_frames,dt,dCM);

                    MSD_R(dt) = mean(SqDisp_R,'all','omitnan');
                    MSD_CORR_R(dt) = mean(SqDisp_CORR_R,'all','omitnan');
                    MSD_CM_R(dt) = mean(SqDisp_CM_R,'all','omitnan');

                    MSD_ALL(dt) = MSD_R(dt);
                    MSD_CORR_ALL(dt) = MSD_CORR_R(dt);
                    MSD_CM_ALL(dt) = MSD_CM_R(dt);

                otherwise

                    % CTRL cells
                    SqDisp_C = GetSqDisp(cell_ALL_trajectories_C_Corr,n_frames,dt);
                    SqDisp_CORR_C = GetSqDisp_CORR(cell_ALL_trajectories_C_Corr,n_frames,dt);
                    [SqDisp_CM_C,~] = GetSqDisp_CM(cell_ALL_trajectories_C_Corr,n_frames,dt,dCM);

                    MSD_C(dt) = mean(SqDisp_C,'all','omitnan');
                    MSD_CORR_C(dt) = mean(SqDisp_CORR_C,'all','omitnan');
                    MSD_CM_C(dt) = mean(SqDisp_CM_C,'all','omitnan');

                    % RAB5 cells
                    SqDisp_R = GetSqDisp(cell_ALL_trajectories_R_Corr,n_frames,dt);
                    SqDisp_CORR_R = GetSqDisp_CORR(cell_ALL_trajectories_R_Corr,n_frames,dt);
                    [SqDisp_CM_R,~] = GetSqDisp_CM(cell_ALL_trajectories_R_Corr,n_frames,dt,dCM);

                    MSD_R(dt) = mean(SqDisp_R,'all','omitnan');
                    MSD_CORR_R(dt) = mean(SqDisp_CORR_R,'all','omitnan');
                    MSD_CM_R(dt) = mean(SqDisp_CM_R,'all','omitnan');


                    % MSD ALL cells
                    SqDisp_ALL = [SqDisp_C; SqDisp_R];
                    SqDisp_CORR_ALL = [SqDisp_CORR_C; SqDisp_CORR_R]; % GetSqDisp_CORR(cell_ALL_trajectories_Corr,n_frames,dt);
                    [SqDisp_CM_ALL,~] = GetSqDisp_CM(cell_ALL_trajectories_Corr,n_frames,dt,dCM);

                    MSD_ALL(dt) = mean(SqDisp_ALL,'all','omitnan');
                    MSD_CORR_ALL(dt) = mean(SqDisp_CORR_ALL,'all','omitnan');
                    MSD_CM_ALL(dt) = mean(SqDisp_CM_ALL,'all','omitnan');

            end

        end

        % Save results
        if save_results == true
            disp('Saving ...')
            save([folder filename_SAVE])
        end

        toc

    end

end

toc
