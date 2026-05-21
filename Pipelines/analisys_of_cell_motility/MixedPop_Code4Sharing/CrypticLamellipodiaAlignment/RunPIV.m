% Script to run PIV on phase contrast movies

clear
close all
clc

addpath('G:\Il mio Drive\UNIMI\MIXED_POPULATIONS_2026\MATLAB_CODE\CODE_PIV\PIV_Lab\PIVlab_CG')
addpath('G:\Il mio Drive\UNIMI\MIXED_POPULATIONS_2026\MATLAB_CODE\CODE_PIV\PIV_Lab\PIVlab\')

display = false;

% List of folders containing movies
subfolder_list = {  'G:\Il mio Drive\UNIMI\Leonardo_CRYPTIC_LAMELLIPODIA\Lamellipodia_per_Chiara_2026\Rac\RAC\20250730\Movies\PhC_PIV'; ...
                    'G:\Il mio Drive\UNIMI\Leonardo_CRYPTIC_LAMELLIPODIA\Lamellipodia_per_Chiara_2026\Rac\RAC\20251008\Movies\PhC_PIV'; ...
                    'G:\Il mio Drive\UNIMI\Leonardo_CRYPTIC_LAMELLIPODIA\Lamellipodia_per_Chiara_2026\Rac\WT\20250723\Movies\PhC_PIV'; ...
                    'G:\Il mio Drive\UNIMI\Leonardo_CRYPTIC_LAMELLIPODIA\Lamellipodia_per_Chiara_2026\Rac\WT\20250730\Movies\PhC_PIV'; ...
                    'G:\Il mio Drive\UNIMI\Leonardo_CRYPTIC_LAMELLIPODIA\Lamellipodia_per_Chiara_2026\Rac\WT\20251008\Movies\PhC_PIV'};

% Loop over the list of folders
for lll = 1 : numel(subfolder_list)

    % Folder containing the movies to analyse
    main_folder = subfolder_list{lll};
    disp(main_folder)

    files = dir([main_folder '\*.tif']);

    tic

    % Find all filenames
    names = cell(numel(files),1);
    for gg = 1 : numel(files)
        names{gg} = files(gg).name;
    end

    % Loop over the movies
    for sss = 1 : numel(names)

        disp(names{sss})

        folderin = main_folder;
        folderout = folderin;

        filename = names{sss};

        for dt = [1 2]
            for win = [64 128]
                disp(['dt' num2str(dt) ', w' num2str(win)])
                PIV_Labber_CG(filename,folderin,folderout,win,dt,display);
                toc
            end
        end

    end

end

toc

