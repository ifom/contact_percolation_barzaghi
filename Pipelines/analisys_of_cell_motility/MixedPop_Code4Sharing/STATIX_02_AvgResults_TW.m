% This scripts imports the cell array "cell_ALL_tables", evaluates the ensemble average of different parameters for a single frame, and then averages over the time window for each sample and FOV
% *avg* indicates the ensemble average in each frame
% *AVG* indicates the final average over the time window
% Run after STATIX_01_VoronoiClusters.m
% C. Guidolin

clear
close all
clc
tic



DataSet = 2021; % The dataset is mostly identified either by the year (YYYY) or the full date (YYYYMMDD) in case of more than one experiment per year
% Read the function GetDatasetInfo for special cases or if in doubt

HD = 'D'; % Hard disk


save_results = false;

% It uses the function GetDatasetInfo.m
addpath([HD ':\GUIDOLIN_CHIARA\MIXED_POPULATIONS\LATEST_CODE_VERSION\Functions\'])

filename_LOAD = '\MAT_STATIX_Results.mat';
filename_SAVE = '\MAT_STATIX_Results_AVG_Time.mat';


[main_folder, sample_name_cell, fov, ~, ~, ~, ~, n_frames] = GetDatasetInfo(DataSet,HD);

n_samples = numel(sample_name_cell);

for sss =  1 : n_samples

    sample_name = sample_name_cell{sss};

    disp(['Sample: ' sample_name])


    for FOV = fov

        disp(['FOV = ' num2str(FOV)])

        % Folder containing results file
        folder_AVG = [main_folder sample_name '\FOV_' num2str(FOV)];


        switch sample_name

            case {'C100','C100_RACDN','R100'}

                load([folder_AVG filename_LOAD],'cell_ALL_tables')
               
                % Preallocation
                AR_avg = zeros(n_frames,1);
                AR_Voronoi_avg = zeros(n_frames,1);
                Area_avg = zeros(n_frames,1);
                Area_Voronoi_avg = zeros(n_frames,1);
                number_nuclei_tot = zeros(n_frames,1);


                for ggg = 1 : n_frames

                    % Index bulk nuclei
                    idx_ALL = cell_ALL_tables{ggg}.idx_Bulk == 1;

                    % Number of nuclei
                    number_nuclei_tot(ggg) = nnz(idx_ALL);

                    % Aspect ratio
                    AR_avg(ggg) = mean(cell_ALL_tables{ggg}.AspectRatio_ALL(idx_ALL));
                    AR_Voronoi_avg(ggg) = mean(cell_ALL_tables{ggg}.AspectRatio_Voronoi(idx_ALL));

                    % Areas nuclei and Voronoi cells
                    Area_avg(ggg) = mean(cell_ALL_tables{ggg}.Area(idx_ALL));
                    Area_Voronoi_avg(ggg) = mean(cell_ALL_tables{ggg}.area_VoronoiCell(idx_ALL));

                end



                % AVERAGE OVER TIME WINDOW

                % Aspect ratio
                AR_AVG = mean(AR_avg);
                AR_Voronoi_AVG = mean(AR_Voronoi_avg,'omitnan');

                Area_AVG = mean(Area_avg);
                Area_Voronoi_AVG = mean(Area_Voronoi_avg,'omitnan');

            otherwise

                load([folder_AVG filename_LOAD],'cell_ALL_tables','isPercolating_C','isPercolating_R','phi_A_MaxCluster_ctrl','phi_A_MaxCluster_RAB5')


                % Preallocation
                AR_avg = zeros(n_frames,1);
                AR_avg_C = zeros(n_frames,1);
                AR_avg_R = zeros(n_frames,1);
                AR_Voronoi_avg = zeros(n_frames,1);
                AR_Voronoi_avg_C = zeros(n_frames,1);
                AR_Voronoi_avg_R = zeros(n_frames,1);
                Area_avg = zeros(n_frames,1);
                Area_avg_C = zeros(n_frames,1);
                Area_avg_R = zeros(n_frames,1);
                Area_Voronoi_avg = zeros(n_frames,1);
                Area_Voronoi_avg_C = zeros(n_frames,1);
                Area_Voronoi_avg_R = zeros(n_frames,1);
                DP_avg = zeros(n_frames,1);
                DP_avg_C = zeros(n_frames,1);
                DP_avg_R = zeros(n_frames,1);
                number_nuclei_tot = zeros(n_frames,1);


                for ggg = 1 : n_frames

                    % Indices
                    idx_ALL = cell_ALL_tables{ggg}.idx_Bulk == 1;
                    idx_C = cell_ALL_tables{ggg}.idx_Bulk == 1 & cell_ALL_tables{ggg}.CellType == 1;
                    idx_R = cell_ALL_tables{ggg}.idx_Bulk == 1 & cell_ALL_tables{ggg}.CellType == 2;

                    % Number of nuclei
                    number_nuclei_tot(ggg) = nnz(idx_ALL);

                    % Aspect ratio nuclei and Voronoi cells
                    AR_avg(ggg) = mean(cell_ALL_tables{ggg}.AspectRatio_ALL(idx_ALL));
                    AR_avg_C(ggg) = mean(cell_ALL_tables{ggg}.AspectRatio_ALL(idx_C));
                    AR_avg_R(ggg) = mean(cell_ALL_tables{ggg}.AspectRatio_ALL(idx_R));
                    AR_Voronoi_avg(ggg) = mean(cell_ALL_tables{ggg}.AspectRatio_Voronoi(idx_ALL));
                    AR_Voronoi_avg_C(ggg) = mean(cell_ALL_tables{ggg}.AspectRatio_Voronoi(idx_C));
                    AR_Voronoi_avg_R(ggg) = mean(cell_ALL_tables{ggg}.AspectRatio_Voronoi(idx_R));

                    % Areas nuclei and Voronoi cells
                    Area_avg(ggg) = mean(cell_ALL_tables{ggg}.Area(idx_ALL));
                    Area_avg_C(ggg) = mean(cell_ALL_tables{ggg}.Area(idx_C));
                    Area_avg_R(ggg) = mean(cell_ALL_tables{ggg}.Area(idx_R));
                    Area_Voronoi_avg(ggg) = mean(cell_ALL_tables{ggg}.area_VoronoiCell(idx_ALL));
                    Area_Voronoi_avg_C(ggg) = mean(cell_ALL_tables{ggg}.area_VoronoiCell(idx_C));
                    Area_Voronoi_avg_R(ggg) = mean(cell_ALL_tables{ggg}.area_VoronoiCell(idx_R));

                    % Experimental (area and number) cell fraction
                    area_tot_C(ggg) = sum(cell_ALL_tables{ggg}.area_VoronoiCell(idx_C),'omitnan');
                    area_tot_R(ggg) = sum(cell_ALL_tables{ggg}.area_VoronoiCell(idx_R),'omitnan');
                    area_tot_ALL(ggg) = area_tot_C(ggg) + area_tot_R(ggg);  % sum(areas_VoronoiCells_ALL(tab_ALL.idx_Bulk == 1),'omitnan');
                    phi_area_C(ggg) = area_tot_C(ggg) / area_tot_ALL(ggg);
                    phi_area_R(ggg) = area_tot_R(ggg) / area_tot_ALL(ggg);
                    n_C = nnz(idx_C);
                    n_R = nnz(idx_R);
                    phi_n_R(ggg) = n_R / (n_C + n_R);
                    phi_n_C(ggg) = n_C / (n_C + n_R);

                    

                    % Correct the DP definition here:
                    
                    % Correct column Demixing Parameter for CTRL cells
                    cell_ALL_tables{ggg}.DemixingParameter(idx_C) = cell_ALL_tables{ggg}.DemixingParameter(idx_C) *(1/2) /(1 - phi_n_C(ggg));
                    % Correct column Demixing Parameter for RAB5 cells
                    cell_ALL_tables{ggg}.DemixingParameter(idx_R) = cell_ALL_tables{ggg}.DemixingParameter(idx_R) *(1/2) /(1 - phi_n_R(ggg));


                    % Demixing parameter
                    DP_avg(ggg) = mean(cell_ALL_tables{ggg}.DemixingParameter(idx_ALL));
                    DP_avg_C(ggg) = mean(cell_ALL_tables{ggg}.DemixingParameter(idx_C));
                    DP_avg_R(ggg) = mean(cell_ALL_tables{ggg}.DemixingParameter(idx_R));

                    % Correlation between Area, AR for C cells and their number of R neighbours
                    if ggg == 1
                        n = cell_ALL_tables{ggg}.n_neighbours(idx_C); % number of neighbours
                        nR = cell_ALL_tables{ggg}.n_neighbours_R(idx_C); % number of RAB5 neighbours
                        nC = cell_ALL_tables{ggg}.n_neighbours_C(idx_C); % number of CTRL neighbours
                        nR_norm = nR ./ n;
                        AR_C = cell_ALL_tables{ggg}.AspectRatio_ALL(idx_C);
                        Area_C = cell_ALL_tables{ggg}.Area(idx_C);
                        AR_Voronoi_C = cell_ALL_tables{ggg}.AspectRatio_Voronoi(idx_C);
                        Area_Voronoi_C = cell_ALL_tables{ggg}.area_VoronoiCell(idx_C);
                    else
                        n = [n; cell_ALL_tables{ggg}.n_neighbours(idx_C)];
                        nR = [nR; cell_ALL_tables{ggg}.n_neighbours_R(idx_C)];
                        nC = [nC; cell_ALL_tables{ggg}.n_neighbours_C(idx_C)];
                        nR_norm = [nR_norm; cell_ALL_tables{ggg}.n_neighbours_R(idx_C)./(cell_ALL_tables{ggg}.n_neighbours_R(idx_C)+cell_ALL_tables{ggg}.n_neighbours_C(idx_C))];
                        AR_C = [AR_C; cell_ALL_tables{ggg}.AspectRatio_ALL(idx_C)];
                        Area_C = [Area_C; cell_ALL_tables{ggg}.Area(idx_C)];
                        AR_Voronoi_C = [AR_Voronoi_C; cell_ALL_tables{ggg}.AspectRatio_Voronoi(idx_C)];
                        Area_Voronoi_C = [Area_Voronoi_C; cell_ALL_tables{ggg}.area_VoronoiCell(idx_C)];
                    end

                end

                % figure(100*sss)
                % plot(DP_avg_C,'b-','DisplayName',['fov=' num2str(FOV)])
                % hold on
                % plot(DP_avg_R,'r-','DisplayName',['fov=' num2str(FOV)])
                % hold on
                % legend show
                % title(replace(sample_name,'_','-'))


                



                % AVERAGE OVER TIME WINDOW

                % Aspect ratio
                AR_AVG = mean(AR_avg,'omitnan');
                AR_AVG_C = mean(AR_avg_C,'omitnan');
                AR_AVG_R = mean(AR_avg_R,'omitnan');
                AR_Voronoi_AVG = mean(AR_Voronoi_avg,'omitnan');
                AR_Voronoi_AVG_C = mean(AR_Voronoi_avg_C,'omitnan');
                AR_Voronoi_AVG_R = mean(AR_Voronoi_avg_R,'omitnan');

                % Areas
                Area_AVG = mean(Area_avg,'omitnan');
                Area_AVG_C = mean(Area_avg_C,'omitnan');
                Area_AVG_R = mean(Area_avg_R,'omitnan');
                Area_Voronoi_AVG = mean(Area_Voronoi_avg,'omitnan');
                Area_Voronoi_AVG_C = mean(Area_Voronoi_avg_C,'omitnan');
                Area_Voronoi_AVG_R = mean(Area_Voronoi_avg_R,'omitnan');

                % Cell fractions
                phi_n_C_AVG = mean(phi_n_C,'omitnan');
                phi_n_R_AVG = mean(phi_n_R,'omitnan');
                phi_area_C_AVG = mean(phi_area_C,'omitnan');
                phi_area_R_AVG = mean(phi_area_R,'omitnan');

                % Demixing parameter
                DP_AVG = mean(DP_avg);
                DP_AVG_C = mean(DP_avg_C);
                DP_AVG_R = mean(DP_avg_R);


                phi_A_MaxCluster_ctrl_AVG = mean(phi_A_MaxCluster_ctrl,'omitnan');
                phi_A_MaxCluster_RAB5_AVG = mean(phi_A_MaxCluster_RAB5,'omitnan');

                % phi_A Percolating cluster ! DEFINIZIONE CAMBIATA PER TENERE GLI ZERI

                %             phi_A_PerCluster_C_AVG = mean(phi_A_MaxCluster_ctrl(isPercolating_C == 1),'omitnan');
                %             phi_A_PerCluster_R_AVG = mean(phi_A_MaxCluster_RAB5(isPercolating_R == 1),'omitnan');
                idx_perc_C = find(isPercolating_C == 1);
                test_C = zeros(n_frames,1);
                test_C(idx_perc_C) = phi_A_MaxCluster_ctrl(idx_perc_C);
                phi_A_PerCluster_C_AVG = mean(test_C);

                idx_perc_R = find(isPercolating_R == 1);
                test_R = zeros(n_frames,1);
                test_R(idx_perc_R) = phi_A_MaxCluster_RAB5(idx_perc_R);
                phi_A_PerCluster_R_AVG = mean(test_R);



        end

        % % Plot number of nuclei over time
        % figure(sss)
        % plot(number_nuclei_tot,'s','DisplayName',['FOV ' num2str(FOV)]);
        % hold on
        % drawnow
        % xlabel('frame')
        % ylabel('N')
        % title(replace(sample_name,'_','-'))
        % legend show

        % Save results
        if save_results == true
            save([folder_AVG filename_SAVE])
        end

        toc

    end

end

toc