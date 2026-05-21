function [main_folder, sample_name_cell, fov, scale, dt, img_size, dist_edge, n_frames] = GetDatasetInfo(DataSet,HD)
% This function is used to automatically load all the relevant experimental info of the specified dataset

arguments
    DataSet (1,1) double {mustBeInteger}
    HD (1,1) char {mustBeTextScalar}
end

switch DataSet

    case 2016

        main_folder = [HD ':\GUIDOLIN_CHIARA\MIXED_POPULATIONS\MCF10A_2016_Analysis\'];
        sample_name_cell = {'C100'; 'C89R11'; 'C80R20'; 'C67R33'; 'C50R50'; 'C33R67'; 'C20R80'; 'C11R89'; 'R100'};
        fov = 1 : 4;
        scale = 660.48/512; % um / pxl
        dt = 10/60; % hours
        img_size = [672, 512]; % this is reversed compared to size(img) as this gives [height, width]=(n_rows,n_cols), but centroid coordinates (x,y) follow std convention (x horizontal, y vertical)
        dist_edge = 30;
        % n_frames = 48; % numel(img_DS_C.Files);
        n_frames = 72; % numel(img_DS_C.Files);

    case 2021

        main_folder = [HD ':\GUIDOLIN_CHIARA\MIXED_POPULATIONS\Leonardo_2021_Analysis\'];
        sample_name_cell = {'C100'; 'C80R20'; 'C60R40'; 'C50R50'; 'C40R60'; 'C20R80'; 'R100'};
        fov = 1 : 5;
        scale = 1040/800; % um / pxl
        dt = 5/60; % hours
        img_size = [800, 800];
        dist_edge = 30;
        n_frames = 96; % numel(img_DS_C.Files);



    case 2022

        % main_folder = 'D:\MIXED_POPULATIONS\Leonardo_2022_Analysis\MixedPop_RACDN_20221202\';
        % sample_name_cell = {'C100'; 'C80R20'; 'C80R20_RACDN'; 'C70R30'; 'C70R30_RACDN'; 'C50R50'; 'C50R50_RACDN'; 'C30R70'; 'C30R70_RACDN'; 'C20R80'; 'C20R80_RACDN'; 'R100'};
        % fov = 1 : 5;
        % scale = 1342.85/1024; % um / pxl
        % dt = 5/60; % hours
        % img_size = [800, 800];
        % dist_edge = 30;
        % n_frames = 96; % numel(img_DS_C.Files);

        main_folder = [HD ':\GUIDOLIN_CHIARA\MIXED_POPULATIONS\Leonardo_2022_Analysis\MixedPop_RACDN_20221202_newTW\'];
        sample_name_cell = {'C100'; 'C80R20'; 'C80R20_RACDN'; 'C70R30'; 'C70R30_RACDN'; 'C50R50'; 'C50R50_RACDN'; 'C30R70'; 'C30R70_RACDN'; 'C20R80'; 'C20R80_RACDN'; 'R100'};
        fov = 1 : 5;
        scale = 1342.85/1024; % um / pxl
        dt = 5/60; % hours
        img_size = [800, 800];
        dist_edge = 30;
        n_frames = 200; % numel(img_DS_C.Files);

    case 202301

        main_folder = [HD ':\GUIDOLIN_CHIARA\MIXED_POPULATIONS\Leonardo_2023_Analysis\MixCtrlBicolor\'];
        sample_name_cell = {'C100'; 'C80R20'; 'C65R35'; 'C50R50'; 'C35R65'; 'C20R80'; 'R100'};
        fov = 1 : 5;
        scale = 1331.2/1024; % um / pxl
        dt = 5/60; % hours
        img_size = [800, 800];
        dist_edge = 30;
        n_frames = 96;

    case 20230217

%         main_folder = 'D:\MIXED_POPULATIONS\Leonardo_2023_Analysis\MixedPop_RACDN_20230217\';
%         sample_name_cell = {'C100_RACDN'; 'C80R20'; 'C80R20_RACDN'; 'C70R30'; 'C70R30_RACDN'; 'C50R50'; 'C50R50_RACDN'; 'C30R70'; 'C30R70_RACDN'; 'C20R80'; 'C20R80_RACDN'; 'R100'};
%         fov = 1 : 5;
%         scale = 1342.85/1024; % um / pxl
%         dt = 5/60; % hours
%         img_size = [800, 800];
%         dist_edge = 30;
%         n_frames = 96;

        main_folder = [HD ':\GUIDOLIN_CHIARA\MIXED_POPULATIONS\Leonardo_2023_Analysis\MixedPop_RACDN_20230217_newTW\'];
        sample_name_cell = {'C100_RACDN'; 'C80R20'; 'C80R20_RACDN'; 'C70R30'; 'C70R30_RACDN'; 'C50R50'; 'C50R50_RACDN'; 'C30R70'; 'C30R70_RACDN'; 'C20R80'; 'C20R80_RACDN'; 'R100'};
        fov = 1 : 5;
        scale = 1342.85/1024; % um / pxl
        dt = 5/60; % hours
        img_size = [800, 800];
        dist_edge = 30;
        n_frames = 200;

    case 20230504

        main_folder = [HD ':\GUIDOLIN_CHIARA\MIXED_POPULATIONS\Leonardo_2023_Analysis\MixedPop_RACDN_20230504\'];
        sample_name_cell = {'C100_RACDN'; 'C80R20'; 'C80R20_RACDN'; 'C70R30'; 'C70R30_RACDN'; 'C50R50'; 'C50R50_RACDN'; 'C30R70'; 'C30R70_RACDN'; 'C20R80'; 'C20R80_RACDN'; 'R100'};
        fov = 1 : 5;
        scale = 1342.85/1024; % um / pxl
        dt = 5/60; % hours
        img_size = [800, 800];
        dist_edge = 30;
        n_frames = 96;

    case 20230612

        main_folder = [HD ':\GUIDOLIN_CHIARA\MIXED_POPULATIONS\Leonardo_2023_Analysis\MixedPop_RACDN_20230612\'];
        sample_name_cell = {'C100'; 'C80R20'; 'C80R20_RACDN'; 'C70R30'; 'C70R30_RACDN'; 'C50R50'; 'C50R50_RACDN'; 'C30R70'; 'C30R70_RACDN'; 'C20R80'; 'C20R80_RACDN'; 'R100'};
        fov = 1 : 5;
        scale = 1342.85/1024; % um / pxl
        dt = 5/60; % hours
        img_size = [800, 800];
        dist_edge = 30;
        n_frames = 96;

    case 20230802

        main_folder = [HD ':\GUIDOLIN_CHIARA\MIXED_POPULATIONS\Leonardo_2023_Analysis\MixedPop_RACDN_20230802\'];
        sample_name_cell = {'C100_RACDN'; 'C80R20'; 'C80R20_RACDN'; 'C65R35'; 'C65R35_RACDN'; 'C50R50'; 'C50R50_RACDN'; 'C35R65'; 'C35R65_RACDN'; 'C20R80'; 'C20R80_RACDN'; 'R100'};
        fov = 1 : 5;
        scale = 1342.85/1024; % um / pxl
        dt = 5/60; % hours
        img_size = [800, 800];
        dist_edge = 30;
        n_frames = 96;

    case 20230906

        main_folder = [HD ':\GUIDOLIN_CHIARA\MIXED_POPULATIONS\Leonardo_2023_Analysis\MixedPop_RACDN_20230906\'];
        sample_name_cell = {'C100'; 'C80R20'; 'C80R20_RACDN'; 'C70R30'; 'C70R30_RACDN'; 'C50R50'; 'C50R50_RACDN'; 'C30R70'; 'C30R70_RACDN'; 'C20R80'; 'C20R80_RACDN'; 'R100'};
        fov = 1 : 5;
        scale = 1339.38/1024; % um / pxl
        dt = 5/60; % hours
        img_size = [800, 800];
        dist_edge = 30;
        n_frames = 96;


    case 2024020200

        main_folder = [HD ':\GUIDOLIN_CHIARA\MIXED_POPULATIONS\Leonardo_2024_Analysis\MixedPop_RACDN_20240202_TW4h12h\'];
        sample_name_cell = {'C80R20'; 'C80R20_RACDN'; 'C65R35'; 'C65R35_RACDN'; 'C50R50'; 'C50R50_RACDN'; 'C35R65'; 'C35R65_RACDN'; 'C20R80'; 'C20R80_RACDN'; 'R100'};
        fov = 1 : 5;
        scale = 1339.38/1024; % um / pxl
        dt = 5/60; % hours
        img_size = [800, 800];
        dist_edge = 30;
        n_frames = 96;

    case 20240202

        main_folder = [HD ':\GUIDOLIN_CHIARA\MIXED_POPULATIONS\Leonardo_2024_Analysis\MixedPop_RACDN_20240202\'];
        sample_name_cell = {'C80R20'; 'C80R20_RACDN'; 'C65R35'; 'C65R35_RACDN'; 'C50R50'; 'C50R50_RACDN'; 'C35R65'; 'C35R65_RACDN'; 'C20R80'; 'C20R80_RACDN'; 'R100'};
        fov = 1 : 5;
        scale = 1339.38/1024; % um / pxl
        dt = 5/60; % hours
        img_size = [800, 800];
        dist_edge = 30;
        n_frames = 200;

    case 20240626 % Camillo

        % main_folder = [HD ':\GUIDOLIN_CHIARA\Timelapse_mix_RAC_DN_26-6-2024\TIF_Analysis\'];
        % sample_name_cell = {'C100_RACDN'; 'C80R20'; 'C80R20_RACDN'; 'C70R30'; 'C70R30_RACDN'; 'C50R50'; 'C50R50_RACDN'; 'C30R70'; 'C30R70_RACDN'; 'C20R80'; 'C20R80_RACDN'; 'R100'};
        % fov = 1 : 6;
        % scale = []; % um / pxl
        % dt = []; % hours
        % img_size = [800, 800];
        % dist_edge = 30;
        % n_frames = 96;

        main_folder = [HD ':\GUIDOLIN_CHIARA\MIXED_POPULATIONS\Timelapse_mix_RAC_DN_26-6-2024\TIF_Analysis_TW200400\'];
        sample_name_cell = {'C100_RACDN'; 'C80R20'; 'C80R20_RACDN'; 'C70R30'; 'C70R30_RACDN'; 'C50R50'; 'C50R50_RACDN'; 'C30R70'; 'C30R70_RACDN'; 'C20R80'; 'C20R80_RACDN'; 'R100'};
        fov = 1 : 6;
        scale = 1331.2/1024; % um / pxl
        dt = 5/60; % hours
        img_size = [800, 800];
        dist_edge = 30;
        n_frames = 200;


    case 20240703 % Camillo

        main_folder = [HD ':\GUIDOLIN_CHIARA\MIXED_POPULATIONS\TIF_20240703_Analysis\'];
        sample_name_cell = {'C100_RACDN'; 'C80R20'; 'C80R20_RACDN'; 'C70R30'; 'C70R30_RACDN'; 'C50R50'; 'C50R50_RACDN'; 'C30R70'; 'C30R70_RACDN'; 'C20R80'; 'C20R80_RACDN'; 'R100'};
        fov = 1 : 5;
        scale = 1497.6/1152; % um / pxl
        dt = 5/60; % hours
        img_size = [1152, 1152];
        dist_edge = 30;
        n_frames = 200;


    case 20240724 % Camillo

        main_folder = [HD ':\GUIDOLIN_CHIARA\MIXED_POPULATIONS\TIF_20240724_Analysis\'];
        sample_name_cell = {'C100_RACDN'; 'C80R20'; 'C80R20_RACDN'; 'C70R30'; 'C70R30_RACDN'; 'C50R50'; 'C50R50_RACDN'; 'C30R70'; 'C30R70_RACDN'; 'C20R80'; 'C20R80_RACDN'; 'R100'};
        fov = 1 : 5;
        scale = 1497.6/1152; % um / pxl
        dt = 5/60; % hours
        img_size = [1152, 1152];
        dist_edge = 30;
        n_frames = 200;

    case 20240719 % Leo

        main_folder = [HD ':\GUIDOLIN_CHIARA\MIXED_POPULATIONS\Leonardo_2024_Analysis\MixedPop_RACDN_20240719\'];
        sample_name_cell = {'C100_RACDN'; 'C80R20'; 'C80R20_RACDN'; 'C70R30'; 'C70R30_RACDN'; 'C50R50'; 'C50R50_RACDN'; 'C30R70'; 'C30R70_RACDN'; 'C20R80'; 'C20R80_RACDN'; 'R100'};
        fov = 1 : 5;
        scale = 1331.2/1024; % um / pxl
        dt = 5/60; % hours
        img_size = [800, 800];
        dist_edge = 30;
        n_frames = 200;



end