% This script evaluates spatial velocity correlations and fits the curves to get the correlation length
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
CalculateVnVm = false;

% It loads:
filename_LOAD = '\MAT_Traj_Corr_ALL.mat';
filename_LOAD_STATIX = '\FinalResults_STATIX.mat'; % to load real phi_R
% It saves:
filename_SAVE = '\VelocityCorrelation.mat';
filename_SAVE_AVG = '\VelocityCorrelation_AVG';

% It uses the functions:
% GetDatasetInfo.m, SelectTrajectories.m, GetVelocityCorrelation.m, BinVelocityCorr.m
addpath([HD ':\GUIDOLIN_CHIARA\MIXED_POPULATIONS\LATEST_CODE_VERSION\Functions\'])

% Load dataset info
[main_folder, sample_name_cell, fov, scale, dt, img_size, dist_edge, n_frames] = GetDatasetInfo(DataSet,HD);



n_samples = numel(sample_name_cell);
n_fov = numel(fov);

% Bin details
delta_R = 10; % pxl
MinR = delta_R + 1;
MaxR = 500;
R_bins = (MinR : delta_R : MaxR)'; % pxl


% % Bin details only for MCF10A
% delta_R = 5; % pxl
% MinR = 10;%delta_R + 1;
% MaxR = 400;
% R_bins = (MinR : delta_R : MaxR)'; % pxl




if CalculateVnVm == true

    % Loop over the samples
    for sss = 1 : n_samples

        sample_name = sample_name_cell{sss};

        disp(['Sample: ' sample_name])

        % Loop over the fields of view
        for FOV = fov

            disp(['FOV = ' num2str(FOV)])

            % Folder containing TrackMate output tables
            folder = [main_folder sample_name '\FOV_' num2str(FOV)];

            disp('Uploading data ..')

            load([folder filename_LOAD],'cell_ALL_trajectories_Corr')

            cell_ALL_traj = cell_ALL_trajectories_Corr;

            min_traj_length = n_frames/2;
            
            % Select trajectories
            disp('Selecting trajectories ..')
            cell_ALL_traj = SelectTrajectories(cell_ALL_traj,min_traj_length);

            % Calculate correlations
            disp('Calculating correlations ..')
            [~, cell_Rnm, cell_VnVm] = GetVelocityCorrelation(cell_ALL_traj,n_frames);
            % Now it is calculated in the lab sdr from corrected trajectories

            toc


            % Bin results

            disp('Binning results ..')

            [V_binned, V_binned_err, V_dR] = BinVelocityCorr(n_frames,R_bins,delta_R,cell_Rnm,cell_VnVm); %,cell_VnVm_X,cell_VnVm_Y);
            
            toc

            


            % Save results
            if save_results == true
                disp('Saving results ..')
                save([folder filename_SAVE],'V_binned', 'V_binned_err', 'V_dR','R_bins')
                disp('Results saved.')
            end
            toc

        end

    end

end







%% AVERAGE OVER TIME AND FOVs

% Load experimental phi_R
load([main_folder filename_LOAD_STATIX],'PHI_R_n')


tic

disp('AVERAGE over time and FOVs')

% Preallocation sss
L_CORR_FIN = zeros(n_samples,1);
L_CORR_FIN_ERR = zeros(n_samples,1);
VCORR_TAB = zeros(numel(R_bins),n_samples);
VCORR_ERR_TAB = zeros(numel(R_bins),n_samples);


for sss = 1 : n_samples

    disp(['Sample: ' sample_name_cell{sss}])

    sample_name = sample_name_cell{sss};
    sname = erase(sample_name,'_');

    folder_AVG_FOV = [main_folder sample_name '\'];


    for FOV = fov
        
        folder = [main_folder sample_name '\FOV_' num2str(FOV)];
        load([folder filename_SAVE],'V_binned')
        

        % Average over the time window for each FOV
        V_binned_TW_avg = mean(V_binned,2);
        V_binned_TW_std = std(V_binned,0,2);
        V_allfov(:,FOV) = V_binned_TW_avg;
        V_allfov_err(:,FOV) = V_binned_TW_std;
        
               
        if save_results == true
            save([folder filename_SAVE_AVG], 'V_binned_TW_avg','V_binned_TW_std','R_bins')
            disp('Results saved.')
        end
        
        % Plot each FOV for each sample
        figure(sss)
        tt = 0 : 0.1 : 500*scale;
        errorbar(R_bins * scale,V_binned_TW_avg * scale^2 / dt^2,V_binned_TW_std * scale^2 / dt^2,'s','DisplayName',[ 'fov ' num2str(FOV)])
        hold on
                
    end

    % Weighted average
    V_binned_FOV_avg = sum(V_allfov .* (1./V_allfov_err.^2), 2) ./ sum(1./V_allfov_err.^2,2);
    V_binned_FOV_err = sqrt( 1./ sum(1./V_allfov_err.^2,2));

    % % Average
    % V_binned_FOV_avg = mean(V_allfov,2);
    % V_binned_FOV_err = std(V_allfov,0,2) / sqrt(numel(fov));

    

    % Stop at smaller R
    idx_fit = R_bins * scale < 250;

    % FIT OF THE AVERAGED CORR FUNCTION

    fit_Lcorr = fittype('A*(exp(-(x/B)^C)) + D'); % stretched exp

    opts = fitoptions(fit_Lcorr);
    opts.StartPoint = [300, 10, .5, 0];
    opts.Lower = [1 0.1 .5 0];
    opts.Upper = [1e3 Inf .5 0];
    opts.Weights = 1./R_bins(idx_fit);% 1./sqrt(R_bins);
   

    
    [fit_V,~] = fit(R_bins(idx_fit) * scale, V_binned_FOV_avg(idx_fit) * scale^2 / dt^2, fit_Lcorr,  opts);
    V_rms_av = sqrt(fit_V.A);
    L_corr = fit_V.B; 
    beta(sss) = fit_V.C; 
    plateau = fit_V.D;
    rr = 1 : 0.1 : 500*scale;


    L_corr_err = 0.2 * L_corr;

    
    figure(100)
    subplot(3,4,sss)
    errorbar(R_bins(idx_fit) * scale,V_binned_FOV_avg(idx_fit) * scale^2 / dt^2,V_binned_FOV_err(idx_fit)*scale^2/dt^2, 's','Color',[.5 .5 .5],'DisplayName',sname,'LineWidth',1.5)
    hold on
    plot(rr,fit_V(rr),'k-','LineWidth',1.5)
    title(sname)
    xlabel('R [\mum]')
    ylabel('\langlev_n\cdotv_m\rangle [\mum^2/h^2]')
    set(gca,'FontSize',10,'FontWeight','bold','LineWidth',1)
    % set(gca,'XScale','log')

    figure(200)
    plot(rr,fit_V(rr),'k-','LineWidth',1.5,'DisplayName',sname)
    hold on
    errorbar(R_bins * scale,V_binned_FOV_avg * scale^2 / dt^2,V_binned_FOV_err*scale^2/dt^2, 's','DisplayName',sname,'LineWidth',1.5)
    xlabel('R [\mum]')
    ylabel('\langlev_n\cdotv_m\rangle [\mum^2/h^2]')
    set(gca,'PlotBoxAspectRatio',[1,1,1],'FontSize',20,'FontWeight','bold','LineWidth',1)

    L_CORR_FIN(sss) = L_corr; % um
    L_CORR_FIN_ERR(sss) = L_corr_err; % um

    VCORR_TAB(:,sss) = V_binned_FOV_avg * scale^2 / dt^2;
    VCORR_ERR_TAB(:,sss) = V_binned_FOV_err * scale^2 / dt^2;

    R_TAB = R_bins * scale;
   
    if save_results == true
        save([folder_AVG_FOV filename_SAVE_AVG],'L_corr','L_corr_err', 'V_binned_FOV_avg','V_binned_FOV_err','R_bins')
        disp('Results saved.')
    end

end
toc

%% Save results

if save_results == true
    % Create and save final results tables
    T_Lcorr = table(sample_name_cell,PHI_R_n,L_CORR_FIN,L_CORR_FIN_ERR,'VariableNames',{'sample','phi_R','L_corr','L_corr_err'});
    aaa = append('VnVm -',sample_name_cell);
    bbb = append('VnVm_err -',sample_name_cell);
    VarNames = [{'dr_um'};aaa;bbb]';
    T_Vcorr = table(R_TAB,VCORR_TAB,VCORR_ERR_TAB);  %,'VariableNames',[{'dt_h'};aaa;bbb]);
    writetable(T_Lcorr,[main_folder '\L_corr_' num2str(DataSet) '.xlsx']);
    writetable(T_Vcorr,[main_folder '\Vcorr_' num2str(DataSet) '.xlsx']); % need to add var names manually
end

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

figure(1000)
errorbar(PHI_R_n(idx_sss),L_CORR_FIN(idx_sss),L_CORR_FIN_ERR(idx_sss),'o','MarkerFaceColor',COLOUR,'MarkerEdgeColor',COLOUR,'MarkerSize',20,'DisplayName',plotNAME)
hold on
xlim([0 100])
% ylim([1 1000])
xlabel('\phi_R [%]')
ylabel('L_{corr} [\mum]')
set(gca,'PlotBoxAspectRatio',[1,1,1],'FontSize',20,'FontWeight','bold','LineWidth',1)
set(gca, 'YScale', 'log')




