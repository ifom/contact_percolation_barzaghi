% This script fits the MSD(dt) curves to get the persistence length
% Run after DYNAMIX_05_MSD.m
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
filename_LOAD_STATIX = '\FinalResults_STATIX.mat';
filename_LOAD_MSD = '\MAT_MSD_AllTraj.mat';     % '\MAT_MSD_LminHalfTW.mat';  ['\MAT_MSD_L' num2str(n_frames) '.mat'];
% It saves:
filename_SAVE_Lpers_FOV = '\MAT_MSD_Lpers_.mat';
filename_SAVE_Lpers_sss = '\Global_MSD_Lpers.mat';


% It uses the functions:
% GetDatasetInfo.m
addpath([HD ':\GUIDOLIN_CHIARA\MIXED_POPULATIONS\LATEST_CODE_VERSION\Functions\'])

% Load dataset info
[main_folder, sample_name_cell, fov, scale, dt, img_size, dist_edge, n_frames] = GetDatasetInfo(DataSet,HD);





% Delta_t = dt*(1 : n_frames/2); % hours
Delta_t = dt*unique(round(logspace(0,log10(round(n_frames/2)),50)));

n_samples = numel(sample_name_cell);
n_fov = numel(fov);

% Load experimental phi_R
load([main_folder filename_LOAD_STATIX],'PHI_R_n')





% idx NORMAL / RACDN samples
idx_racdn = contains(sample_name_cell,'RACDN');



% Preallocation
L_PERS = NaN(n_samples,1); % Final Lpers - one length for each sample
L_PERS_err = NaN(n_samples,1); % Error of the average
L_PERS_avgMSD = NaN(n_samples,1); % Final Lpers - one length for each sample
D_avgMSD = NaN(n_samples,1);
gamma_avgMSD = NaN(n_samples,1);
% MSD_sss = zeros(n_frames-1,n_samples); % Final MSD - one column for each sample
% MSD_sss_err = zeros(n_frames-1,n_samples); % Error of the average - one column for each sample

% Loop over the samples
for sss = 1 : n_samples     %[1 6 n_samples]

    sample_name = sample_name_cell{sss};

    disp(['Sample: ' sample_name])

    folder_SAMPLE = [main_folder sample_name];



    % Preallocation FOV
    % MSD_fov = NaN(n_frames-1,n_fov);
    D = NaN(n_fov,1);
    gamma = NaN(n_fov,1);
    v_0 = NaN(n_fov,1);
    L_pers = NaN(n_fov,1);
    

    % Loop over the fields of view
    for FOV = fov

        disp(['FOV = ' num2str(FOV)])

        folder_FOV = [main_folder sample_name '\FOV_' num2str(FOV)];
        
        load([folder_FOV filename_LOAD_MSD],'MSD_CORR_ALL');%,'sel_delays');
        
        % Delta_t = dt * sel_delays;
        
        MSD_fit = nonzeros(MSD_CORR_ALL);
        MSD_fov(:,FOV) = nonzeros(MSD_CORR_ALL);
        
        



        % FIT MSD curves
                
        % New fit:
        new_fit = fittype('2*A*x + 2*A/B * (exp(-B*x)-1) + C');
        opts = fitoptions(new_fit);
        opts.StartPoint = [1e2, 1, 0.5]; %[StartPoint_A(sss), StartPoint_B(sss), 1];
        opts.Lower = [10 .1 0.3];
        opts.Upper = [1e4 10 1];
        opts.Weight = 1 ./ Delta_t.^3; % 1 ./ Delta_t;
        
        [fit_MSD,~] = fit(Delta_t',MSD_fit(1:numel(Delta_t))*scale^2,new_fit,opts);
        D(FOV) = fit_MSD.A;
        gamma(FOV) = fit_MSD.B;
        % baseline = fit_MSD.C;
        
        v_0(FOV) = sqrt(D(FOV) * gamma(FOV));
        L_pers(FOV) = v_0(FOV) / gamma(FOV);
        
        tt = 0.05:0.001:10;
        
    end

    

    MSD = mean(MSD_fov(:,fov),2);
    MSD_err = std(MSD_fov(:,fov),0,2) / sqrt(numel(fov));

    % Average MSD over FOVs
    MSD_sss(:,sss) = MSD * scale^2;
    MSD_sss_err(:,sss) = MSD_err * scale^2;

    L_PERS(sss) = mean(L_pers);

    

    % Fit averaged MSD
    [fit_MSD_allfov,~] = fit(Delta_t',MSD(1:numel(Delta_t))*scale^2,new_fit,opts);
    

    D_avgMSD(sss) = fit_MSD_allfov.A;
    gamma_avgMSD(sss) = fit_MSD_allfov.B;
    baseline_avgMSD(sss) = fit_MSD_allfov.C;

    % Persistence length obtained from the fit of the averaged MSD
    L_PERS_avgMSD(sss) =  sqrt(fit_MSD_allfov.A * fit_MSD_allfov.B) / fit_MSD_allfov.B;

    L_PERS_err(sss) = max( std(L_pers) / sqrt(numel(fov)), 0.10 * L_PERS_avgMSD(sss) );



    if ~contains(sample_name,'RACDN')
        figure(1)
        errorbar(Delta_t,MSD(1:numel(Delta_t))*scale^2,MSD_err(1:numel(Delta_t)) * scale^2,'s','DisplayName',sample_name)
        hold on
        plot(tt,fit_MSD_allfov(tt),'k-','DisplayName',[sample_name ' fit'])
        xlabel('\Deltat [h]')
        ylabel('MSD(\Deltat) [\mum^2]')
        legend show
        set(gca,'PlotBoxAspectRatio',[1,1,1],'FontSize',20,'FontWeight','bold','LineWidth',1)
        set(gca, 'YScale', 'log','XScale','log')
    elseif contains(sample_name,'RACDN')
        figure(2)
        errorbar(Delta_t,MSD(1:numel(Delta_t))*scale^2,MSD_err(1:numel(Delta_t)) * scale^2,'s','DisplayName',replace(sample_name,'_','-'))
        hold on
        plot(tt,fit_MSD_allfov(tt),'k-','DisplayName',[replace(sample_name,'_','-') ' fit'])
        xlabel('\Deltat [h]')
        ylabel('MSD(\Deltat) [\mum^2]')
        legend show
        set(gca,'PlotBoxAspectRatio',[1,1,1],'FontSize',20,'FontWeight','bold','LineWidth',1)
        set(gca, 'YScale', 'log','XScale','log')
    end


    % Save results
    if save_results == true
        save([folder_SAMPLE filename_SAVE_Lpers_FOV],'MSD_fov','L_pers','v_0','gamma','D','MSD','MSD_err')
        save([main_folder filename_SAVE_Lpers_sss],'L_PERS','L_PERS_err','L_PERS_avgMSD','MSD_sss','MSD_sss_err')
    end

end



%% Plot results

% idx NORMAL / RACDN samples
idx_RACDN = contains(sample_name_cell,'RACDN');
idx_NORMAL = ~idx_RACDN;
idx_R100 = contains(sample_name_cell,'R100');

idx_sss = idx_NORMAL;


if all(contains(sample_name_cell(idx_sss),'RACDN'))
    plotNAME = ['EXP ' num2str(DataSet) ' - RACDN'];
    COLOUR = 'm';
else
    plotNAME = ['EXP ' num2str(DataSet) ' - NORMAL'];
    COLOUR = 'k';
end

figure(100)
subplot(1,2,1)
errorbar(PHI_R_n(idx_sss),L_PERS_avgMSD(idx_sss),L_PERS_err(idx_sss),'s','Color',COLOUR,'MarkerEdgeColor',COLOUR,'MarkerFaceColor',COLOUR,'MarkerSize',20,'DisplayName',plotNAME)
hold on
xlim([0 100])
ylim([0 300])
xlabel('\phi_R [%]')
ylabel('L_{pers} [\mum]')
set(gca,'PlotBoxAspectRatio',[1,1,1],'FontSize',20,'FontWeight','bold','LineWidth',1)
% set(gca, 'YScale', 'log')
legend show


subplot(1,2,2)
% errorbar(PHI_R_n(idx_sss),L_PERS(idx_sss)/L_PERS(idx_R100),L_PERS_err(idx_sss),'s','Color',COLOUR,'MarkerEdgeColor',COLOUR,'MarkerFaceColor',COLOUR,'MarkerSize',20,'DisplayName','NORMAL')
plot(PHI_R_n(idx_sss),L_PERS_avgMSD(idx_sss)/L_PERS_avgMSD(idx_R100),'s','Color',COLOUR,'MarkerEdgeColor',COLOUR,'MarkerFaceColor',COLOUR,'MarkerSize',20,'DisplayName',plotNAME)
hold on
xlim([0 100])
ylim([0 1])
xlabel('\phi_R [%]')
ylabel('L_{pers}/L_{pers}^{R100}')
set(gca,'PlotBoxAspectRatio',[1,1,1],'FontSize',20,'FontWeight','bold','LineWidth',1)
legend show

% %%
% % 
% % Create final results tables
% T_Lpers = table(sample_name_cell,PHI_R_n,L_PERS_avgMSD,L_PERS_err,'VariableNames',{'sample','phi_R','L_pers','L_pers_err'});
% aaa = append('MSD -',sample_name_cell);
% bbb = append('MSD_err -',sample_name_cell);
% VarNames = [{'dt_h'};aaa;bbb]';
% T_MSD = table(Delta_t',MSD_sss(1:length(Delta_t),:),MSD_sss_err(1:length(Delta_t),:));%,'VariableNames',VarNames);
% writetable(T_Lpers,[main_folder '\L_pers_' num2str(DataSet) '.xlsx']);
% writetable(T_MSD,[main_folder '\MSD_' num2str(DataSet) '.xlsx']); % need to add var names manually


