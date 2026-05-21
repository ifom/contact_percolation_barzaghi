% It averages the results (already averaged over the TW) over the different FOVs
% Run after STATIX_02_AvgResults_TW.m
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
filename_LOAD = '\MAT_STATIX_Results_AVG_Time.mat';
% It saves:
filename_SAVE = '\Results_AVG_FOV_STATIX.mat';


addpath([HD ':\GUIDOLIN_CHIARA\MIXED_POPULATIONS\LATEST_CODE_VERSION\Functions\'])

% Load dataset info
[main_folder, sample_name_cell, fov, ~, ~, ~, ~, ~] = GetDatasetInfo(DataSet,HD);


% Loop over the samples
for sss = 1 : numel(sample_name_cell)

    sample_name = sample_name_cell{sss};

    disp(['Sample: ' sample_name])

    folder_save = [main_folder sample_name '\'];


    % Preallocation
    number = zeros(numel(fov),1);
    % Aspect ratio
    AR_FOV = zeros(numel(fov),1);
    AR_FOV_C = zeros(numel(fov),1);
    AR_FOV_R = zeros(numel(fov),1);
    AR_Voronoi_FOV = zeros(numel(fov),1);
    AR_Voronoi_FOV_C = zeros(numel(fov),1);
    AR_Voronoi_FOV_R = zeros(numel(fov),1);
    % Areas
    Area_FOV = zeros(numel(fov),1);
    Area_FOV_C = zeros(numel(fov),1);
    Area_FOV_R = zeros(numel(fov),1);
    Area_Voronoi_FOV = zeros(numel(fov),1);
    Area_Voronoi_FOV_C = zeros(numel(fov),1);
    Area_Voronoi_FOV_R = zeros(numel(fov),1);
    % Demixing parameter
    DP_FOV = zeros(numel(fov),1);
    DP_FOV_C = zeros(numel(fov),1);
    DP_FOV_R = zeros(numel(fov),1);
    % Exp cell fractions
    phi_n_R_exp_FOV = zeros(numel(fov),1);
    phi_A_R_exp_FOV = zeros(numel(fov),1);
    phi_n_C_exp_FOV = zeros(numel(fov),1);
    phi_A_C_exp_FOV = zeros(numel(fov),1);
    % Percolating cluster
    phi_A_PerCluster_C_FOV = zeros(numel(fov),1);
    phi_A_PerCluster_R_FOV = zeros(numel(fov),1);


    % Loop over the fields of view
    for FOV = fov

        disp(['FOV=' num2str(FOV)])

        folder_AVG = [main_folder sample_name '\FOV_' num2str(FOV)];

        load([folder_AVG filename_LOAD],'*_AVG_*','*_AVG','*_binned','Pearson*','pval*','number_nuclei_tot','n','n*','AR*','Area*')

        switch sample_name

            case {'C100','C100_RACDN','R100'}

                number(FOV) = mean(number_nuclei_tot,'omitnan');

                % STRUCTURE

                AR_FOV(FOV) = AR_AVG;

                AR_Voronoi_FOV(FOV) = AR_Voronoi_AVG;

                Area_FOV(FOV) = Area_AVG;

                Area_Voronoi_FOV(FOV) = Area_Voronoi_AVG;



            otherwise

                % Correlation with local environment
                if FOV == 1
                    n_FOV = n;
                    nR_FOV = nR;
                    nC_FOV = nC;
                    nR_norm_FOV = nR_norm;
                    AR_C_FOV = AR_C;
                    Area_C_FOV = Area_C;
                    AR_Voronoi_C_FOV = AR_Voronoi_C;
                    Area_Voronoi_C_FOV = Area_Voronoi_C;
                else
                    n_FOV = [n_FOV; n];
                    nR_FOV = [nR_FOV; nR];
                    nC_FOV = [nC_FOV; nC];
                    nR_norm_FOV = [nR_norm_FOV; nR_norm];
                    AR_C_FOV = [AR_C_FOV; AR_C];
                    Area_C_FOV = [Area_C_FOV; Area_C];
                    AR_Voronoi_C_FOV = [AR_Voronoi_C_FOV; AR_Voronoi_C];
                    Area_Voronoi_C_FOV = [Area_Voronoi_C_FOV; Area_Voronoi_C];
                end

                number(FOV) = mean(number_nuclei_tot,'omitnan');

                % STRUCTURE

                AR_FOV(FOV) = AR_AVG;
                AR_FOV_C(FOV) = AR_AVG_C;
                AR_FOV_R(FOV) = AR_AVG_R;

                AR_Voronoi_FOV(FOV) = AR_Voronoi_AVG;
                AR_Voronoi_FOV_C(FOV) = AR_Voronoi_AVG_C;
                AR_Voronoi_FOV_R(FOV) = AR_Voronoi_AVG_R;

                Area_FOV(FOV) = Area_AVG;
                Area_FOV_C(FOV) = Area_AVG_C;
                Area_FOV_R(FOV) = Area_AVG_R;

                Area_Voronoi_FOV(FOV) = Area_Voronoi_AVG;
                Area_Voronoi_FOV_C(FOV) = Area_Voronoi_AVG_C;
                Area_Voronoi_FOV_R(FOV) = Area_Voronoi_AVG_R;

                DP_FOV(FOV) = DP_AVG;
                DP_FOV_C(FOV) = DP_AVG_C;
                DP_FOV_R(FOV) = DP_AVG_R;

                phi_n_C_exp_FOV(FOV) = phi_n_C_AVG;
                phi_n_R_exp_FOV(FOV) = phi_n_R_AVG;
                phi_A_C_exp_FOV(FOV) = phi_area_C_AVG;
                phi_A_R_exp_FOV(FOV) = phi_area_R_AVG;

                phi_A_PerCluster_C_FOV(FOV) = phi_A_PerCluster_C_AVG;
                phi_A_PerCluster_R_FOV(FOV) = phi_A_PerCluster_R_AVG;


        end

    end


    % Averages over the FOV

    switch sample_name

        case {'C100','C100_RACDN','R100'}

            number_tot = mean(number(fov),'omitnan');
            number_tot_std = std(number(fov));

            AR_FOV_avg = mean(AR_FOV(fov));

            AR_Voronoi_FOV_avg = mean(AR_Voronoi_FOV(fov));

            Area_FOV_avg = mean(Area_FOV(fov));

            Area_Voronoi_FOV_avg = mean(Area_Voronoi_FOV(fov));

        otherwise

            number_tot = mean(number(fov),'omitnan');
            number_tot_std = std(number(fov));

            % Structure
            phi_n_C_exp_FOV_avg = mean(phi_n_C_exp_FOV(fov));
            phi_n_R_exp_FOV_avg = mean(phi_n_R_exp_FOV(fov));
            phi_A_C_exp_FOV_avg = mean(phi_A_C_exp_FOV(fov));
            phi_A_R_exp_FOV_avg = mean(phi_A_R_exp_FOV(fov));
            phi_A_PerCluster_C_FOV_avg = mean(phi_A_PerCluster_C_FOV(fov));
            phi_A_PerCluster_R_FOV_avg = mean(phi_A_PerCluster_R_FOV(fov));

            AR_FOV_avg = mean(AR_FOV(fov));
            AR_FOV_avg_C = mean(AR_FOV_C(fov));
            AR_FOV_avg_R = mean(AR_FOV_R(fov));

            AR_Voronoi_FOV_avg = mean(AR_Voronoi_FOV(fov));
            AR_Voronoi_FOV_avg_C = mean(AR_Voronoi_FOV_C(fov));
            AR_Voronoi_FOV_avg_R = mean(AR_Voronoi_FOV_R(fov));

            Area_FOV_avg = mean(Area_FOV(fov));
            Area_FOV_avg_C = mean(Area_FOV_C(fov));
            Area_FOV_avg_R = mean(Area_FOV_R(fov));

            Area_Voronoi_FOV_avg = mean(Area_Voronoi_FOV(fov));
            Area_Voronoi_FOV_avg_C = mean(Area_Voronoi_FOV_C(fov));
            Area_Voronoi_FOV_avg_R = mean(Area_Voronoi_FOV_R(fov));

            DP_FOV_avg = mean(DP_FOV(fov));
            DP_FOV_avg_C = mean(DP_FOV_C(fov));
            DP_FOV_avg_R = mean(DP_FOV_R(fov));


    end

    % Save results
    if save_results == true
        save([folder_save filename_SAVE], '*_FOV_*','*_FOV','number_tot','number_tot_std')
    end

    toc

end