% This script builds all trajectories of CTRL and RAB5 cells from TrackMate output CSV files for each FOV and sample
% The trajectories are stored and saved in a cell array called "cell_ALL_trajectories"
% C. Guidolin

clear
close all
clc
tic

DataSet = 2021; % The dataset is mostly identified either by the year (YYYY) or the full date (YYYYMMDD) in case of more than one experiment per year
% Read the function GetDatasetInfo for special cases or if in doubt

HD = 'D'; % Hard disk

save_results = false;


% It saves:
filename_SAVE_C = '\MAT_Traj_C.mat';
filename_SAVE_R = '\MAT_Traj_R.mat';
filename_SAVE_ALL = '\MAT_Traj_ALL.mat';

% It uses the functions:
% GetDatasetInfo.m; GetTrajectories.m
addpath([HD ':\GUIDOLIN_CHIARA\MIXED_POPULATIONS\LATEST_CODE_VERSION\Functions\'])

% Load dataset info
[main_folder, sample_name_cell, fov, ~, ~, ~, ~, ~] = GetDatasetInfo(DataSet,HD);


n_samples = numel(sample_name_cell);

% Loop over the samples
for sss = 1 : n_samples

    sample_name = sample_name_cell{sss};

    disp(['Sample: ' sample_name])

    % Loop over the fields of view
    for FOV = fov

        disp(['FOV=' num2str(FOV)])

        % Folder containing TrackMate output tables
        folder = [main_folder sample_name '\FOV_' num2str(FOV)];

        % Import tracking tables (Trackmate output)

        switch sample_name

            case {'C100','C100_RACDN'} % Only CTRL

                % Tracking tables of CTRL cells + remove first (empty) row
                tab_C_spots = readtable([folder '\TrackMate_C_spots.csv']);     tab_C_spots(1,:) = [];
                tab_C_tracks = readtable([folder '\TrackMate_C_tracks.csv']);    tab_C_tracks(1,:) = [];

                % Add +1 to frame number in TrackMate spot table (as it starts from 0 instead of 1)
                tab_C_spots.FRAME = tab_C_spots.FRAME + 1;


            case 'R100'  % Only RAB5

                % Tracking tables of RAB5 cells + remove first (empty) row
                tab_R_spots = readtable([folder '\TrackMate_R_spots.csv']);      tab_R_spots(1,:) = [];
                tab_R_tracks = readtable([folder '\TrackMate_R_tracks.csv']);    tab_R_tracks(1,:) = [];

                % Add +1 to frame number in TrackMate spot table (as it starts from 0 instead of 1)
                tab_R_spots.FRAME = tab_R_spots.FRAME + 1;


            otherwise  % Mixed samples

                % Tracking tables of CTRL cells + remove first (empty) row
                tab_C_spots = readtable([folder '\TrackMate_C_spots.csv']);      tab_C_spots(1,:) = [];
                tab_C_tracks = readtable([folder '\TrackMate_C_tracks.csv']);    tab_C_tracks(1,:) = [];

                % Tracking tables of RAB5 cells + remove first (empty) row
                tab_R_spots = readtable([folder '\TrackMate_R_spots.csv']);      tab_R_spots(1,:) = [];
                tab_R_tracks = readtable([folder '\TrackMate_R_tracks.csv']);    tab_R_tracks(1,:) = [];

                % Add +1 to frame number in Trackmate spot table (as it starts from 0 instead of 1)
                tab_C_spots.FRAME = tab_C_spots.FRAME + 1;
                tab_R_spots.FRAME = tab_R_spots.FRAME + 1;

        end


        % Find trajectories

        switch sample_name

            case {'C100','C100_RACDN'}

                % Only CTRL cells
                [fraction_trajectories_C, cell_ALL_trajectories_C, tab_C_spots, ID_tracking_C] = GetTrajectories(tab_C_spots, tab_C_tracks);

                tab_ALL_spots = tab_C_spots;
                cell_ALL_trajectories = cell_ALL_trajectories_C;

            case 'R100'

                % Only RAB5 cells
                [fraction_trajectories_R, cell_ALL_trajectories_R, tab_R_spots, ID_tracking_R] = GetTrajectories(tab_R_spots, tab_R_tracks);

                tab_ALL_spots = tab_R_spots;
                cell_ALL_trajectories = cell_ALL_trajectories_R;

            otherwise

                % CTRL cells
                [fraction_trajectories_C, cell_ALL_trajectories_C, tab_C_spots, ID_tracking_C] = GetTrajectories(tab_C_spots, tab_C_tracks);

                % RAB5 cells
                [fraction_trajectories_R, cell_ALL_trajectories_R, tab_R_spots, ID_tracking_R] = GetTrajectories(tab_R_spots, tab_R_tracks);

                % ALL together
                tab_ALL_spots = [tab_C_spots; tab_R_spots];
                cell_ALL_trajectories = [cell_ALL_trajectories_C; cell_ALL_trajectories_R];

        end


        % Save results

        if save_results == true

            disp('Saving results ..')

            switch sample_name

                case {'C100','C100_RACDN'}

                    save([folder filename_SAVE_C],'sample_name','fraction_trajectories_C', 'cell_ALL_trajectories_C', 'tab_C_spots','ID_tracking_C')
                    save([folder filename_SAVE_ALL],'sample_name','cell_ALL_trajectories','tab_ALL_spots')

                case 'R100'

                    save([folder filename_SAVE_R],'sample_name','fraction_trajectories_R', 'cell_ALL_trajectories_R', 'tab_R_spots','ID_tracking_R')
                    save([folder filename_SAVE_ALL],'sample_name','cell_ALL_trajectories' ,'tab_ALL_spots')

                otherwise

                    save([folder filename_SAVE_C],'sample_name','fraction_trajectories_C', 'cell_ALL_trajectories_C', 'tab_C_spots','ID_tracking_C')
                    save([folder filename_SAVE_R],'sample_name','fraction_trajectories_R', 'cell_ALL_trajectories_R', 'tab_R_spots','ID_tracking_R')
                    save([folder filename_SAVE_ALL],'sample_name','cell_ALL_trajectories','tab_ALL_spots' )

            end

        end

        toc

    end

end