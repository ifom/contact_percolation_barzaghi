% This script corrects (with a light smoothing) the trajectories for fake displacements due to oscillations of the microscope stage
% Run after DYNAMIX_02_Velocities.m
% C. Guidolin

clear
close all
clc

tic


DataSet = 2021; % The dataset is mostly identified either by the year (YYYY) or the full date (YYYYMMDD) in case of more than one experiment per year
% Read the function GetDatasetInfo for special cases or if in doubt

HD = 'D'; % Hard disk

save_results = false;


% It loads:
filename_LOAD_Traj_C = '\MAT_Traj_C.mat';
filename_LOAD_Traj_R = '\MAT_Traj_R.mat';
filename_LOAD_Velocities = '\MAT_Vrms_Vcm_psi.mat';%'\MAT_Vrms_Vcm_psi_new.mat';
% It saves:
filename_SAVE_C = '\MAT_Traj_Corr_C.mat';
filename_SAVE_R = '\MAT_Traj_Corr_R.mat';
filename_SAVE_ALL = '\MAT_Traj_Corr_ALL.mat';


% It uses the functions:
% GetDatasetInfo.m; CorrectTrajectories.m
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

            case {'C100','C100_RACDN'} % Only CTRL cells

                load([folder filename_LOAD_Traj_C],'cell_ALL_trajectories_C');
                load([folder filename_LOAD_Velocities],'V_fake_C');

                cell_ALL_trajectories_C_Corr = CorrectTrajectories(cell_ALL_trajectories_C,V_fake_C);

                cell_ALL_trajectories_Corr = cell_ALL_trajectories_C_Corr;

                if save_results == true
                    disp('Saving results ..')
                    save([folder filename_SAVE_C], "cell_ALL_trajectories_C_Corr",'V_fake_C')
                    save([folder filename_SAVE_ALL], "cell_ALL_trajectories_Corr")
                end

            case 'R100'  % Only RAB5 cells

                load([folder filename_LOAD_Traj_R],'cell_ALL_trajectories_R');
                load([folder filename_LOAD_Velocities],'V_fake_R');

                cell_ALL_trajectories_R_Corr = CorrectTrajectories(cell_ALL_trajectories_R,V_fake_R);

                cell_ALL_trajectories_Corr = cell_ALL_trajectories_R_Corr;

                if save_results == true
                    disp('Saving results ..')
                    save([folder filename_SAVE_R], "cell_ALL_trajectories_R_Corr",'V_fake_R')
                    save([folder filename_SAVE_ALL], "cell_ALL_trajectories_Corr")
                end

            otherwise % Mixed samples

                load([folder filename_LOAD_Velocities],'V_fake_C','V_fake_R','V_fake');

                % CTRL

                load([folder filename_LOAD_Traj_C],'cell_ALL_trajectories_C');

                cell_ALL_trajectories_C_Corr = CorrectTrajectories(cell_ALL_trajectories_C,V_fake_C);

                % RAB5

                load([folder filename_LOAD_Traj_R],'cell_ALL_trajectories_R');

                cell_ALL_trajectories_R_Corr = CorrectTrajectories(cell_ALL_trajectories_R,V_fake_R);

                % All trajectories together

                cell_ALL_trajectories = [cell_ALL_trajectories_C; cell_ALL_trajectories_R];
                cell_ALL_trajectories_Corr = CorrectTrajectories(cell_ALL_trajectories,V_fake);


                % Save results
                if save_results == true
                    disp('Saving results ..')
                    save([folder filename_SAVE_C], "cell_ALL_trajectories_C_Corr",'V_fake_C')
                    save([folder filename_SAVE_R], "cell_ALL_trajectories_R_Corr",'V_fake_R')
                    save([folder filename_SAVE_ALL], "cell_ALL_trajectories_Corr",'V_fake')
                end

        end

    end
    toc

end

toc
