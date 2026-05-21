% This script groups and plots final averaged STATIX results for a given dataset
% Run after STATIX_03_AvgResults_FOV.m
% C. Guidolin


clear
close all
clc
tic



DataSet = 20240202; % The dataset is mostly identified either by the year (YYYY) or the full date (YYYYMMDD) in case of more than one experiment per year
% Read the function GetDatasetInfo for special cases or if in doubt

HD = 'D'; % Hard disk

save_results = false;
save_figures = false;
PlotResults = true;



% It loads:
filename_LOAD = '\Results_AVG_FOV_STATIX.mat';
% It saves:
filename_SAVE_SSS = '\FinalResults_STATIX.mat';


addpath([HD ':\GUIDOLIN_CHIARA\MIXED_POPULATIONS\LATEST_CODE_VERSION\Functions\'])


% Load dataset info
[main_folder, sample_name_cell, ~, scale, ~, ~, ~, ~] = GetDatasetInfo(DataSet,HD);

n_samples = numel(sample_name_cell);

folder_SAVE_fig = [main_folder '\Figures_' num2str(DataSet) '\'];
mkdir(folder_SAVE_fig);

% Experimental fractions
PHI_R_n = zeros(n_samples,1);       PHI_R_n_err = zeros(n_samples,1);
PHI_R_A = zeros(n_samples,1);       PHI_R_A_err = zeros(n_samples,1);

% Percolating Cluster
PHI_Perc_C = zeros(n_samples,1);        PHI_Perc_C_err = zeros(n_samples,1);
PHI_Perc_R = zeros(n_samples,1);        PHI_Perc_R_err = zeros(n_samples,1);

% Loop over the samples
for sss = 1 : n_samples

    sample_name = sample_name_cell{sss};

    switch sample_name

        case {'C100','C100_RACDN'}
            PHI_R_n(sss) = 0;
            PHI_R_A(sss) = 0;
            PHI_Perc_C(sss) = 100;

        case 'R100'
            PHI_R_n(sss) = 100;
            PHI_R_A(sss) = 100;
            PHI_Perc_R(sss) = 100;

        otherwise

            folder_results_AVG = [main_folder sample_name_cell{sss} '\'];
            load([folder_results_AVG filename_LOAD],'phi_n_R_exp_FOV_avg','phi_n_R_exp_FOV','phi_A_R_exp_FOV_avg','phi_A_R_exp_FOV')

            PHI_R_n(sss) = phi_n_R_exp_FOV_avg * 100;
            PHI_R_A(sss) = phi_A_R_exp_FOV_avg * 100;
            PHI_R_n_err(sss) = std(phi_n_R_exp_FOV) / sqrt(numel(phi_n_R_exp_FOV)) * 100;
            PHI_R_A_err(sss) = std(phi_A_R_exp_FOV) / sqrt(numel(phi_A_R_exp_FOV)) * 100;

    end
end

% Comparison between nominal and experimental value
% figure(100)
% plot(0:100,0:100,'k--')
% hold on
% plot(phi_R_nom,phi_R_exp_A,'ms','MarkerFaceColor','m','MarkerSize',20)
% errorbar(phi_R_nom,phi_R_exp_A,phi_R_exp_A_std,'k','LineWidth',1.5)
% plot(phi_R_nom,phi_R_exp_n,'gs','MarkerFaceColor','g','MarkerSize',20)
% errorbar(phi_R_nom,phi_R_exp_n,phi_R_exp_n_std,'k','LineWidth',1.5)
% xlim([0 100])
% xlabel('\phi_R^{nom} [%]')
% ylabel('\phi_R^{A}, \phi_R^{n} [%]')
% set(gca,'PlotBoxAspectRatio',[1,1,1],'FontSize',20,'FontWeight','bold','LineWidth',1)


% Preallocation final STATIX parameters
% Area nuclei
AREA = zeros(n_samples,1);              AREA_err = zeros(n_samples,1);
AREA_C = zeros(n_samples,1);            AREA_C_err = zeros(n_samples,1);
AREA_R = zeros(n_samples,1);            AREA_R_err = zeros(n_samples,1);

% Area Voronoi
AREA_VOR = zeros(n_samples,1);          AREA_VOR_err = zeros(n_samples,1);
AREA_VOR_C = zeros(n_samples,1);        AREA_VOR_C_err = zeros(n_samples,1);
AREA_VOR_R = zeros(n_samples,1);        AREA_VOR_R_err = zeros(n_samples,1);

% Aspect Ratio nuclei
AR = zeros(n_samples,1);                AR_err = zeros(n_samples,1);
AR_C = zeros(n_samples,1);              AR_C_err = zeros(n_samples,1);
AR_R = zeros(n_samples,1);              AR_R_err = zeros(n_samples,1);

% Aspect Ratio Voronoi
AR_VOR = zeros(n_samples,1);            AR_VOR_err = zeros(n_samples,1);
AR_VOR_C = zeros(n_samples,1);          AR_VOR_C_err = zeros(n_samples,1);
AR_VOR_R = zeros(n_samples,1);          AR_VOR_R_err = zeros(n_samples,1);


% Demixing Parameter
DP = zeros(n_samples,1);                DP_err = zeros(n_samples,1);
DP_C = zeros(n_samples,1);              DP_C_err = zeros(n_samples,1);
DP_R = zeros(n_samples,1);              DP_R_err = zeros(n_samples,1);




%% AVG results

for sss = 1 : n_samples

    sample_name = sample_name_cell{sss};

    folder_results_AVG = [main_folder sample_name_cell{sss} '\'];

    load([folder_results_AVG filename_LOAD])


    switch sample_name

        case {'C100','C100_RACDN','R100'}

            % Number nuclei
            N_nuclei(sss) = number_tot;
            N_nuclei_err(sss) = number_tot_std / sqrt(numel(N_nuclei));

            % Area nuclei
            AREA(sss) = Area_FOV_avg;
            AREA_C(sss) = NaN;
            AREA_R(sss) = NaN;
            % Error bars
            AREA_err(sss) = std(Area_FOV) / sqrt(numel(Area_FOV));
            AREA_R_err(sss) = NaN;
            AREA_C_err(sss) = NaN;

            % Area Voronoi
            AREA_VOR(sss) = Area_Voronoi_FOV_avg;
            AREA_VOR_C(sss) = NaN;
            AREA_VOR_R(sss) = NaN;
            % Error bars
            AREA_VOR_err(sss) = std(Area_Voronoi_FOV) / sqrt(numel(Area_Voronoi_FOV));
            AREA_VOR_C_err(sss) = NaN;
            AREA_VOR_R_err(sss) = NaN;

            % Aspect Ratio nuclei
            AR(sss) = AR_FOV_avg;
            AR_C(sss) = NaN;
            AR_R(sss) = NaN;
            % Error bars
            AR_err(sss) = std(AR_FOV) / sqrt(numel(AR_FOV));
            AR_C_err(sss) = NaN;
            AR_R_err(sss) = NaN;

            % Aspect Ratio Voronoi
            AR_VOR(sss) = AR_Voronoi_FOV_avg;
            AR_VOR_C(sss) = NaN;
            AR_VOR_R(sss) = NaN;
            % Error bars
            AR_VOR_err(sss) = std(AR_Voronoi_FOV) / sqrt(numel(AR_Voronoi_FOV));
            AR_VOR_C_err(sss) = NaN;
            AR_VOR_R_err(sss) = NaN;

            % Demixing parameter
            DP_R(sss) = NaN;
            DP_R_err(sss) = NaN;
            DP_C(sss) = NaN;
            DP_C_err(sss) = NaN;

        otherwise

            N_nuclei(sss) = number_tot;
            N_nuclei_err(sss) = number_tot_std / sqrt(numel(N_nuclei));

            % Area nuclei
            AREA(sss) = Area_FOV_avg;
            AREA_C(sss) = Area_FOV_avg_C;
            AREA_R(sss) = Area_FOV_avg_R;
            % Error bars
            AREA_err(sss) = std(Area_FOV) / sqrt(numel(Area_FOV));
            AREA_C_err(sss) = std(Area_FOV_C) / sqrt(numel(Area_FOV_C));
            AREA_R_err(sss) = std(Area_FOV_R) / sqrt(numel(Area_FOV_R));

            % Area Voronoi
            AREA_VOR(sss) = Area_Voronoi_FOV_avg;
            AREA_VOR_C(sss) = Area_Voronoi_FOV_avg_C;
            AREA_VOR_R(sss) = Area_Voronoi_FOV_avg_R;
            % Error bars
            AREA_VOR_err(sss) = std(Area_Voronoi_FOV) / sqrt(numel(Area_Voronoi_FOV));
            AREA_VOR_C_err(sss) = std(Area_Voronoi_FOV_C) / sqrt(numel(Area_Voronoi_FOV_C));
            AREA_VOR_R_err(sss) = std(Area_Voronoi_FOV_R) / sqrt(numel(Area_Voronoi_FOV_R));

            % Aspect Ratio nuclei
            AR(sss) = AR_FOV_avg;
            AR_C(sss) = AR_FOV_avg_C;
            AR_R(sss) = AR_FOV_avg_R;
            % Error bars
            AR_err(sss) = std(AR_FOV) / sqrt(numel(AR_FOV));
            AR_C_err(sss) = std(AR_FOV_C) / sqrt(numel(AR_FOV_C));
            AR_R_err(sss) = std(AR_FOV_R) / sqrt(numel(AR_FOV_R));

            % Aspect Ratio Voronoi
            AR_VOR(sss) = AR_Voronoi_FOV_avg;
            AR_VOR_C(sss) = AR_Voronoi_FOV_avg_C;
            AR_VOR_R(sss) = AR_Voronoi_FOV_avg_R;
            % Error bars
            AR_VOR_err(sss) = std(AR_Voronoi_FOV) / sqrt(numel(AR_Voronoi_FOV));
            AR_VOR_C_err(sss) = std(AR_Voronoi_FOV_C) / sqrt(numel(AR_Voronoi_FOV_C));
            AR_VOR_R_err(sss) = std(AR_Voronoi_FOV_R) / sqrt(numel(AR_Voronoi_FOV_R));

            % Demixing Parameter
            DP(sss) = DP_FOV_avg;
            DP_C(sss) = DP_FOV_avg_C;
            DP_R(sss) = DP_FOV_avg_R;
            % Error bars
            DP_err(sss) = std(DP_FOV) / sqrt(numel(DP_FOV));
            DP_C_err(sss) = std(DP_FOV_C) / sqrt(numel(DP_FOV_C));
            DP_R_err(sss) = std(DP_FOV_R) / sqrt(numel(DP_FOV_R));


            % Percolating Cluster
            PHI_Perc_C(sss) = phi_A_PerCluster_C_FOV_avg * 100;
            PHI_Perc_R(sss) = phi_A_PerCluster_R_FOV_avg * 100;
            % Error bars
            PHI_Perc_C_err(sss) = std(phi_A_PerCluster_C_FOV) * 100 / sqrt(numel(phi_A_PerCluster_C_FOV));
            PHI_Perc_R_err(sss) = std(phi_A_PerCluster_R_FOV) * 100 / sqrt(numel(phi_A_PerCluster_R_FOV));

    end

    if save_results == true
        disp('Saving results ..')
        % Save final results for each sample
        save([main_folder filename_SAVE_SSS],'PHI_R_n','sample_name_cell')
    end

end


%% PLOT RESULTS

if PlotResults == true

    % idx NORMAL / RACDN samples
    idx_RACDN = contains(sample_name_cell,'RACDN');
    idx_NORMAL = ~idx_RACDN;

    % Convert areas to real units (um^2)
    AREA = AREA * scale^2;
    AREA_C = AREA_C * scale^2;
    AREA_R = AREA_R * scale^2;
    AREA_err = AREA_err * scale^2;
    AREA_C_err = AREA_C_err * scale^2;
    AREA_R_err = AREA_R_err * scale^2;

    AREA_VOR = AREA_VOR * scale^2;
    AREA_VOR_C = AREA_VOR_C * scale^2;
    AREA_VOR_R = AREA_VOR_R * scale^2;
    AREA_VOR_err = AREA_VOR_err * scale^2;
    AREA_VOR_C_err = AREA_VOR_C_err * scale^2;
    AREA_VOR_R_err = AREA_VOR_R_err * scale^2;

    % Create final results tables
    T_AREAS = table(sample_name_cell,PHI_R_n,AREA,AREA_err,AREA_C,AREA_C_err,AREA_R,AREA_R_err,AREA_VOR,AREA_VOR_err,AREA_VOR_C,AREA_VOR_C_err,AREA_VOR_R,AREA_VOR_R_err,'VariableNames',{'sample','phi_R','Area','Area_err','Area_C','Area_C_err','Area_R','Area_R_err','Area_Voronoi','Area_Voronoi_err','Area_Voronoi_C','Area_Voronoi_C_err','Area_Voronoi_R','Area_Voronoi_R_err'});
    T_AR = table(sample_name_cell,PHI_R_n,AR,AR_err,AR_C,AR_C_err,AR_R,AR_R_err,AR_VOR,AR_VOR_err,AR_VOR_C,AR_VOR_C_err,AR_VOR_R,AR_VOR_R_err,'VariableNames',{'sample','phi_R','AR','AR_err','AR_C','AR_C_err','AR_R','AR_R_err','AR_Voronoi','AR_Voronoi_err','AR_Voronoi_C','AR_Voronoi_C_err','AR_Voronoi_R','AR_Voronoi_R_err'});
    T_DP = table(sample_name_cell,PHI_R_n,DP,DP_err,DP_C,DP_C_err,DP_R,DP_R_err,'VariableNames',{'sample','phi_R','DP','DP_err','DP_C','DP_C_err','DP_R','DP_R_err'});
    T_PhiPerc = table(sample_name_cell,PHI_R_n,PHI_Perc_R,PHI_Perc_R_err,PHI_Perc_C,PHI_Perc_C_err,'VariableNames',{'sample','phi_R','PhiPerc_R','PhiPerc_R_err','PhiPérc_C','PhiPerc_C_err'});

    if save_results == true
        % Save tables
        writetable(T_AREAS,[main_folder '\Areas_' num2str(DataSet) '.xlsx']);
        writetable(T_AR,[main_folder '\AR_' num2str(DataSet) '.xlsx']);
        writetable(T_DP,[main_folder '\DP_' num2str(DataSet) '.xlsx']);
        writetable(T_PhiPerc,[main_folder '\PhiPercolatingCluster_' num2str(DataSet) '.xlsx']);
    end




    % Area nuclei and Voronoi cells

    f1 = figure(1);%('Position',[1 100 2000 1000]);
    sgtitle('Area','FontSize',25,'FontWeight','bold')

    subplot(1,2,1)
    errorbar(PHI_R_n(idx_NORMAL),AREA(idx_NORMAL),AREA_err(idx_NORMAL),'ks','MarkerFaceColor','k','MarkerSize',20,'LineWidth',1.5,'DisplayName','ALL')
    hold on
    errorbar(PHI_R_n(idx_NORMAL),AREA_C(idx_NORMAL),AREA_C_err(idx_NORMAL),'bs','MarkerFaceColor','b','MarkerSize',20,'LineWidth',1.5,'DisplayName','CTRL')
    errorbar(PHI_R_n(idx_NORMAL),AREA_R(idx_NORMAL),AREA_R_err(idx_NORMAL),'rs','MarkerFaceColor','r','MarkerSize',20,'LineWidth',1.5,'DisplayName','RAB5')
    xlim([0 100])
    xlabel('\phi_R [%]')
    ylabel('\langleA\rangle [\mum^2]')
    title('Nuclei')
    set(gca,'PlotBoxAspectRatio',[1,1,1],'FontSize',20,'FontWeight','bold','LineWidth',1)
    legend('show','Location','bestoutside')


    subplot(1,2,2)
    errorbar(PHI_R_n(idx_NORMAL),AREA_VOR(idx_NORMAL),AREA_VOR_err(idx_NORMAL),'kv','MarkerFaceColor','k','MarkerSize',20,'LineWidth',1.5,'DisplayName','ALL')
    hold on
    errorbar(PHI_R_n(idx_NORMAL),AREA_VOR_C(idx_NORMAL),AREA_VOR_C_err(idx_NORMAL),'bv','MarkerFaceColor','b','MarkerSize',20,'LineWidth',1.5,'DisplayName','CTRL')
    errorbar(PHI_R_n(idx_NORMAL),AREA_VOR_R(idx_NORMAL),AREA_VOR_R_err(idx_NORMAL),'rv','MarkerFaceColor','r','MarkerSize',20,'LineWidth',1.5,'DisplayName','RAB5')
    xlim([0 100])
    xlabel('\phi_R [%]')
    ylabel('\langleA\rangle [\mum^2]')
    title('Voronoi cells')
    set(gca,'PlotBoxAspectRatio',[1,1,1],'FontSize',20,'FontWeight','bold','LineWidth',1)
    legend('show','Location','bestoutside')

    set(gcf,'WindowState','maximized')



    % AR nuclei and Voronoi cells

    f2 = figure(2);%('Position',[1 100 2000 1000]);
    sgtitle('Aspect Ratio','FontSize',25,'FontWeight','bold')

    subplot(1,2,1)  % nuclei
    errorbar(PHI_R_n(idx_NORMAL),AR(idx_NORMAL),AR_err(idx_NORMAL),'ks','MarkerFaceColor','k','MarkerSize',20,'LineWidth',1.5,'DisplayName','ALL')
    hold on
    errorbar(PHI_R_n(idx_NORMAL),AR_C(idx_NORMAL),AR_C_err(idx_NORMAL),'bs','MarkerFaceColor','b','MarkerSize',20,'LineWidth',1.5,'DisplayName','CTRL')
    errorbar(PHI_R_n(idx_NORMAL),AR_R(idx_NORMAL),AR_R_err(idx_NORMAL),'rs','MarkerFaceColor','r','MarkerSize',20,'LineWidth',1.5,'DisplayName','RAB5')
    xlim([0 100])
    xlabel('\phi_R [%]')
    ylabel('\langleAR\rangle')
    title('Nuclei')
    set(gca,'PlotBoxAspectRatio',[1,1,1],'FontSize',20,'FontWeight','bold','LineWidth',1)
    legend('show','Location','bestoutside')

    subplot(1,2,2)  % Voronoi
    errorbar(PHI_R_n(idx_NORMAL),AR_VOR(idx_NORMAL),AR_VOR_err(idx_NORMAL),'kv','MarkerFaceColor','k','MarkerSize',20,'LineWidth',1.5,'DisplayName','ALL')
    hold on
    errorbar(PHI_R_n(idx_NORMAL),AR_VOR_C(idx_NORMAL),AR_VOR_C_err(idx_NORMAL),'bv','MarkerFaceColor','b','MarkerSize',20,'LineWidth',1.5,'DisplayName','CTRL')
    errorbar(PHI_R_n(idx_NORMAL),AR_VOR_R(idx_NORMAL),AR_VOR_R_err(idx_NORMAL),'rv','MarkerFaceColor','r','MarkerSize',20,'LineWidth',1.5,'DisplayName','RAB5')
    xlim([0 100])
    xlabel('\phi_R [%]')
    ylabel('\langleAR\rangle')
    title('Voronoi cells')
    set(gca,'PlotBoxAspectRatio',[1,1,1],'FontSize',20,'FontWeight','bold','LineWidth',1)
    legend('show','Location','bestoutside')

    set(gcf,'WindowState','maximized')



    % Demixing Parameter
    f3 = figure(3);%('Position',[100 100 1000 1000]);
    errorbar(PHI_R_n(idx_NORMAL),DP(idx_NORMAL),DP_err(idx_NORMAL),'ks','MarkerFaceColor','k','MarkerSize',20,'LineWidth',1.5,'DisplayName','ALL')
    hold on
    errorbar(PHI_R_n(idx_NORMAL),DP_C(idx_NORMAL),DP_C_err(idx_NORMAL),'bs','MarkerFaceColor','b','MarkerSize',20,'LineWidth',1.5,'DisplayName','CTRL')
    errorbar(PHI_R_n(idx_NORMAL),DP_R(idx_NORMAL),DP_R_err(idx_NORMAL),'rs','MarkerFaceColor','r','MarkerSize',20,'LineWidth',1.5,'DisplayName','RAB5')
    xlim([0 100])
    ylim([0 1])
    xlabel('\phi_R [%]')
    ylabel('\langleDP\rangle')
    title('Demixing Parameter')
    set(gca,'PlotBoxAspectRatio',[1,1,1],'FontSize',20,'FontWeight','bold','LineWidth',1)
    legend('show','Location','bestoutside')
    set(gcf,'WindowState','maximized')


    % Area fraction of percolating cluster
    f4 = figure(4);%('Position',[100 100 1000 1000]);
    errorbar(PHI_R_n(idx_NORMAL),PHI_Perc_C(idx_NORMAL),PHI_Perc_C_err(idx_NORMAL),'bs-','MarkerFaceColor','b','MarkerSize',20,'LineWidth',1.5,'DisplayName','\phi^A_C')
    hold on
    errorbar(PHI_R_n(idx_NORMAL),PHI_Perc_R(idx_NORMAL),PHI_Perc_R_err(idx_NORMAL),'rs-','MarkerFaceColor','r','MarkerSize',20,'LineWidth',1.5,'DisplayName','\phi^A_R')
    xline(50,'k--','DisplayName','\phi=50%')
    xlim([0 100])
    xlabel('\phi_R [%]')
    ylabel('\phi^A_{PC} [%]')
    title('Percolating Cluster')
    set(gca,'PlotBoxAspectRatio',[1,1,1],'FontSize',20,'FontWeight','bold','LineWidth',1)
    legend('show','Location','bestoutside')
    set(gcf,'WindowState','maximized')

    % Total number of nuclei
    f5 = figure(5);
    errorbar(PHI_R_n(idx_NORMAL),N_nuclei(idx_NORMAL),N_nuclei_err(idx_NORMAL),'ks','MarkerFaceColor','k','MarkerSize',20,'LineWidth',1.5,'DisplayName','NORMAL')
    hold on
    xlim([0 100])
    xlabel('\phi_R [%]')
    ylabel('N')
    set(gca,'PlotBoxAspectRatio',[1,1,1],'FontSize',20,'FontWeight','bold','LineWidth',1)
    legend('show','Location','bestoutside')
    set(gcf,'WindowState','maximized')


    if DataSet == 20240202 || DataSet == 2022 || DataSet == 20230217 || DataSet == 20230504 || DataSet == 20230612 || DataSet == 20230802 || DataSet == 20230906 || DataSet == 20240703 || DataSet == 20240626 || DataSet == 20240724 || DataSet == 20240719

        %         Total number of nuclei
        f5 = figure(5);
        errorbar(PHI_R_n(idx_RACDN),N_nuclei(idx_RACDN),N_nuclei_err(idx_RACDN),'s','MarkerFaceColor',[1 0.5 0.5],'MarkerEdgeColor',[1 0.5 0.5],'MarkerSize',20,'LineWidth',1.5,'DisplayName','RACDN')
        hold on
        xlim([0 100])
        xlabel('\phi_R [%]')
        ylabel('N')
        set(gca,'PlotBoxAspectRatio',[1,1,1],'FontSize',20,'FontWeight','bold','LineWidth',1)
        legend('show','Location','bestoutside')
        set(gcf,'WindowState','maximized')

        % Area nuclei and Voronoi cells

        f11 = figure(11);%('Position',[1 100 2000 1000]);
        sgtitle('Area - RACDN','FontSize',25,'FontWeight','bold')

        subplot(1,2,1)
        errorbar(PHI_R_n(idx_RACDN),AREA(idx_RACDN),AREA_err(idx_RACDN),'ks','MarkerFaceColor','k','MarkerSize',20,'LineWidth',1.5,'DisplayName','ALL')
        hold on
        errorbar(PHI_R_n(idx_RACDN),AREA_C(idx_RACDN),AREA_C_err(idx_RACDN),'bs','MarkerFaceColor','b','MarkerSize',20,'LineWidth',1.5,'DisplayName','CTRL')
        errorbar(PHI_R_n(idx_RACDN),AREA_R(idx_RACDN),AREA_R_err(idx_RACDN),'rs','MarkerFaceColor','r','MarkerSize',20,'LineWidth',1.5,'DisplayName','RAB5')
        xlim([0 100])
        xlabel('\phi_R [%]')
        ylabel('\langleA\rangle [\mum^2]')
        title('Nuclei')
        set(gca,'PlotBoxAspectRatio',[1,1,1],'FontSize',20,'FontWeight','bold','LineWidth',1)
        legend('show','Location','bestoutside')
        
        subplot(1,2,2)
        errorbar(PHI_R_n(idx_RACDN),AREA_VOR(idx_RACDN),AREA_VOR_err(idx_RACDN),'kv','MarkerFaceColor','k','MarkerSize',20,'LineWidth',1.5,'DisplayName','ALL')
        hold on
        errorbar(PHI_R_n(idx_RACDN),AREA_VOR_C(idx_RACDN),AREA_VOR_C_err(idx_RACDN),'bv','MarkerFaceColor','b','MarkerSize',20,'LineWidth',1.5,'DisplayName','CTRL')
        errorbar(PHI_R_n(idx_RACDN),AREA_VOR_R(idx_RACDN),AREA_VOR_R_err(idx_RACDN),'rv','MarkerFaceColor','r','MarkerSize',20,'LineWidth',1.5,'DisplayName','RAB5')
        xlim([0 100])
        xlabel('\phi_R [%]')
        ylabel('\langleA\rangle [\mum^2]')
        title('Voronoi cells')
        set(gca,'PlotBoxAspectRatio',[1,1,1],'FontSize',20,'FontWeight','bold','LineWidth',1)
        legend('show','Location','bestoutside')

        set(gcf,'WindowState','maximized')


        % AR nuclei and Voronoi cells

        f22 = figure(22);%figure('Position',[1 100 2000 1000]);
        sgtitle('Aspect Ratio - RACDN','FontSize',25,'FontWeight','bold')

        subplot(1,2,1)  % nuclei
        errorbar(PHI_R_n(idx_RACDN),AR(idx_RACDN),AR_err(idx_RACDN),'ks','MarkerFaceColor','k','MarkerSize',20,'LineWidth',1.5,'DisplayName','ALL')
        hold on
        errorbar(PHI_R_n(idx_RACDN),AR_C(idx_RACDN),AR_C_err(idx_RACDN),'bs','MarkerFaceColor','b','MarkerSize',20,'LineWidth',1.5,'DisplayName','CTRL')
        errorbar(PHI_R_n(idx_RACDN),AR_R(idx_RACDN),AR_R_err(idx_RACDN),'rs','MarkerFaceColor','r','MarkerSize',20,'LineWidth',1.5,'DisplayName','RAB5')
        xlim([0 100])
        xlabel('\phi_R [%]')
        ylabel('\langleAR\rangle')
        title('Nuclei')
        set(gca,'PlotBoxAspectRatio',[1,1,1],'FontSize',20,'FontWeight','bold','LineWidth',1)
        legend('show','Location','bestoutside')

        subplot(1,2,2)  % Voronoi
        errorbar(PHI_R_n(idx_RACDN),AR_VOR(idx_RACDN),AR_VOR_err(idx_RACDN),'kv','MarkerFaceColor','k','MarkerSize',20,'LineWidth',1.5,'DisplayName','ALL')
        hold on
        errorbar(PHI_R_n(idx_RACDN),AR_VOR_C(idx_RACDN),AR_VOR_C_err(idx_RACDN),'bv','MarkerFaceColor','b','MarkerSize',20,'LineWidth',1.5,'DisplayName','CTRL')
        errorbar(PHI_R_n(idx_RACDN),AR_VOR_R(idx_RACDN),AR_VOR_R_err(idx_RACDN),'rv','MarkerFaceColor','r','MarkerSize',20,'LineWidth',1.5,'DisplayName','RAB5')
        xlim([0 100])
        xlabel('\phi_R [%]')
        ylabel('\langleAR\rangle')
        title('Voronoi cells')
        set(gca,'PlotBoxAspectRatio',[1,1,1],'FontSize',20,'FontWeight','bold','LineWidth',1)
        legend('show','Location','bestoutside')

        set(gcf,'WindowState','maximized')




        % Demixing Parameter
        f33 = figure(33);%('Position',[100 100 1000 1000]);
        errorbar(PHI_R_n(idx_RACDN),DP(idx_RACDN),DP_err(idx_RACDN),'ks','MarkerFaceColor','k','MarkerSize',20,'LineWidth',1.5,'DisplayName','ALL')
        hold on
        errorbar(PHI_R_n(idx_RACDN),DP_C(idx_RACDN),DP_C_err(idx_RACDN),'bs','MarkerFaceColor','b','MarkerSize',20,'LineWidth',1.5,'DisplayName','CTRL')
        errorbar(PHI_R_n(idx_RACDN),DP_R(idx_RACDN),DP_R_err(idx_RACDN),'rs','MarkerFaceColor','r','MarkerSize',20,'LineWidth',1.5,'DisplayName','RAB5')
        xlim([0 100])
        ylim([0 1])
        xlabel('\phi_R [%]')
        ylabel('\langleDP\rangle')
        title('Demixing Parameter - RACDN')
        set(gca,'PlotBoxAspectRatio',[1,1,1],'FontSize',20,'FontWeight','bold','LineWidth',1)
        legend('show','Location','bestoutside')
        set(gcf,'WindowState','maximized')


        % Area fraction of percolating cluster
        f44 = figure(44);%('Position',[100 100 1000 1000]);
        errorbar(PHI_R_n(idx_RACDN),PHI_Perc_C(idx_RACDN),PHI_Perc_C_err(idx_RACDN),'bs-','MarkerFaceColor','b','MarkerSize',20,'LineWidth',1.5,'DisplayName','\phi^A_C')
        hold on
        errorbar(PHI_R_n(idx_RACDN),PHI_Perc_R(idx_RACDN),PHI_Perc_R_err(idx_RACDN),'rs-','MarkerFaceColor','r','MarkerSize',20,'LineWidth',1.5,'DisplayName','\phi^A_R')
        xline(50,'k--','DisplayName','\phi=50%')
        xlim([0 100])
        xlabel('\phi_R [%]')
        ylabel('\phi^A_{PC} [%]')
        title('Percolating Cluster - RACDN')
        set(gca,'PlotBoxAspectRatio',[1,1,1],'FontSize',20,'FontWeight','bold','LineWidth',1)
        legend('show','Location','bestoutside')
        set(gcf,'WindowState','maximized')

    end




    % Save Figures
    if save_figures == true

        savefig(f1,[folder_SAVE_fig 'Areas.fig'])
        print(f1,[folder_SAVE_fig 'Areas'],'-dpng','-r0')

        savefig(f2,[folder_SAVE_fig 'AspectRatio.fig'])
        print(f2,[folder_SAVE_fig 'AspectRatio'],'-dpng','-r0')

        savefig(f3,[folder_SAVE_fig 'DemixingParameter_' num2str(DataSet) '.fig'])
        print(f3,[folder_SAVE_fig 'DemixingParameter_' num2str(DataSet)],'-dpng','-r0')

        savefig(f4,[folder_SAVE_fig 'PhiPercolatingCluster_' num2str(DataSet) '.fig'])
        print(f4,[folder_SAVE_fig 'PhiPercolatingCluster_' num2str(DataSet)],'-dpng','-r0')

        if DataSet == 20240202 || DataSet == 2022 || DataSet == 20230217 || DataSet == 20230504 || DataSet == 20230612 || DataSet == 20230802 || DataSet == 20230906 || DataSet == 20240703 || DataSet == 20240626 || DataSet == 20240724 || DataSet == 20240719

            savefig(f11,[folder_SAVE_fig 'Areas_RACDN.fig'])
            print(f11,[folder_SAVE_fig 'Areas_RACDN'],'-dpng','-r0')

            savefig(f22,[folder_SAVE_fig 'AspectRatio_RACDN.fig'])
            print(f22,[folder_SAVE_fig 'AspectRatio_RACDN'],'-dpng','-r0')

            savefig(f33,[folder_SAVE_fig 'DemixingParameter_RACDN_' num2str(DataSet) '.fig'])
            print(f33,[folder_SAVE_fig 'DemixingParameter_RACDN_' num2str(DataSet)],'-dpng','-r0')

            savefig(f44,[folder_SAVE_fig 'PhiPercolatingCluster_RACDN_' num2str(DataSet) '.fig'])
            print(f44,[folder_SAVE_fig 'PhiPercolatingCluster_RACDN_' num2str(DataSet)],'-dpng','-r0')
        end

    end

end

