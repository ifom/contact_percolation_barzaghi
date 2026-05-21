% This script runs the PIV analysis on phase contrast movies
% C. Guidolin

clear
close all
clc

addpath('D:\CODE_PIV\PIV_Lab\PIVlab_CG') % it uses the function PIV_Labber_CG.m
addpath('D:\CODE_PIV\PIV_Lab\PIVlab')

display = false;


% Folder containing the Phase Contrast movies
main_folder = 'D:\MIXED_POPULATIONS\Leonardo_2023_Analysis\MixedPop_RACDN_20230612\PIV\PhC\';


% Get list of all movies
files = dir([main_folder '\*.tif']);

tic

% Find all filenames
names = cell(numel(files),1);
for gg = 1 : numel(files)
    names{gg} = files(gg).name;
end

% Loop over the movies
for sss = 1 : numel(names)

    folderin = main_folder;
    folderout = folderin;

    filename = names{sss};

    disp(filename)

    % Run PIV
    for dt = 1 % [1, 2]

        for win = 16 % [16, 32, 64]

            disp(['dt = ' num2str(dt) ', w = ' num2str(win)])
            PIV_Labber_CG(filename,folderin,folderout,win,dt,display);
            toc

        end
    end

end

toc
