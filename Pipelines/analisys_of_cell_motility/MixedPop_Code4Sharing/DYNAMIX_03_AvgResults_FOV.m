% This script averages the velocities over the different FOVs
% It saves the final results for each sample and plots them
% Run after DYNAMIX_02_Velocities.m
% C. Guidolin

clear
close all
clc
tic

DataSet = 20240703; % The dataset is mostly identified either by the year (YYYY) or the full date (YYYYMMDD) in case of more than one experiment per year
% Read the function GetDatasetInfo for special cases or if in doubt

HD = 'D'; % Hard disk

addpath([HD ':\GUIDOLIN_CHIARA\MIXED_POPULATIONS\LATEST_CODE_VERSION\Functions\'])

save_results = false;
save_figures = false;
PlotResults = true;
CalculateAverages = false;


% It loads:
filename_LOAD = '\MAT_Vrms_Vcm_psi.mat';%'\MAT_Vrms_Vcm_psi_new.mat';
filename_LOAD_STATIX = '\FinalResults_STATIX.mat';
% It saves:
filename_SAVE_AVG = '\Results_AVG_FOV_DYNAMIX.mat';
filename_SAVE_SSS = '\FinalResults_DYNAMIX.mat';%'\FinalResults_DYNAMIX_new.mat';





% Load dataset info
[main_folder, sample_name_cell, fov, scale, dt, img_size, dist_edge, n_frames] = GetDatasetInfo(DataSet,HD);


n_samples = numel(sample_name_cell);
n_fov = numel(fov);


if CalculateAverages == true

    % Preallocation sss

    VCM = NaN(n_samples,1);
    VCM_C = NaN(n_samples,1);
    VCM_R = NaN(n_samples,1);
    VCM_err = NaN(n_samples,1);
    VCM_C_err = NaN(n_samples,1);
    VCM_R_err = NaN(n_samples,1);

    VRMS = NaN(n_samples,1);
    VRMS_C = NaN(n_samples,1);
    VRMS_R = NaN(n_samples,1);
    VRMS_err = NaN(n_samples,1);
    VRMS_C_err = NaN(n_samples,1);
    VRMS_R_err = NaN(n_samples,1);

    VRMS_LAB = NaN(n_samples,1);
    VRMS_LAB_C = NaN(n_samples,1);
    VRMS_LAB_R = NaN(n_samples,1);
    VRMS_LAB_err = NaN(n_samples,1);
    VRMS_LAB_C_err = NaN(n_samples,1);
    VRMS_LAB_R_err = NaN(n_samples,1);

    VCM_SMOOTH = NaN(n_samples,1);
    VCM_SMOOTH_C = NaN(n_samples,1);
    VCM_SMOOTH_R = NaN(n_samples,1);
    VCM_SMOOTH_err = NaN(n_samples,1);
    VCM_SMOOTH_C_err = NaN(n_samples,1);
    VCM_SMOOTH_R_err = NaN(n_samples,1);

    PSI = NaN(n_samples,1);
    PSI_C = NaN(n_samples,1);
    PSI_R = NaN(n_samples,1);
    PSI_err = NaN(n_samples,1);
    PSI_C_err = NaN(n_samples,1);
    PSI_R_err = NaN(n_samples,1);


    % Loop over the samples
    for sss = 1 : n_samples

        sample_name = sample_name_cell{sss};

        disp(['Sample: ' sample_name])

        folder_SAMPLE = [main_folder sample_name];


        % Preallocation FOV
        
        V_CM_mod_FOV = NaN(n_fov,1);
        V_CM_mod_FOV_C = NaN(n_fov,1);
        V_CM_mod_FOV_R = NaN(n_fov,1);

        V_rms_FOV = NaN(n_fov,1);
        V_rms_FOV_C = NaN(n_fov,1);
        V_rms_FOV_R = NaN(n_fov,1);

        V_rms_lab_CORR_FOV = NaN(n_fov,1);
        V_rms_lab_CORR_FOV_C = NaN(n_fov,1);
        V_rms_lab_CORR_FOV_R = NaN(n_fov,1);

        V_CM_mod_smooth_FOV = NaN(n_fov,1);
        V_CM_mod_smooth_FOV_C = NaN(n_fov,1);
        V_CM_mod_smooth_FOV_R = NaN(n_fov,1);

        psi_FOV = NaN(n_fov,1);
        psi_FOV_C = NaN(n_fov,1);
        psi_FOV_R = NaN(n_fov,1);


        % Loop over the fields of view
        for FOV = fov

            disp(['FOV=' num2str(FOV)])

            folder_FOV = [main_folder sample_name '\FOV_' num2str(FOV)];

            switch sample_name

                case {'C100','C100_RACDN','R100'}

                    load([folder_FOV filename_LOAD],'V*','psi_*')

                    V_CM_mod_FOV(FOV) = V_CM_mod_AVG;

                    V_CM_mod_smooth_FOV(FOV) = V_CM_mod_AVG_smooth;

                    V_rms_FOV(FOV) = V_rms_AVG; % in the CM sdr

                    V_rms_lab_CORR_FOV(FOV) = V_rms_lab_CORR_AVG;

                    psi_FOV(FOV) = psi_CORR_AVG;

                otherwise % Mixed samples

                    load([folder_FOV filename_LOAD],'V*','psi_*')

                    V_CM_mod_FOV(FOV) = V_CM_mod_AVG;
                    V_CM_mod_FOV_C(FOV) = V_CM_mod_AVG_C;
                    V_CM_mod_FOV_R(FOV) = V_CM_mod_AVG_R;

                    V_CM_mod_smooth_FOV(FOV) = V_CM_mod_AVG_smooth;
                    V_CM_mod_smooth_FOV_C(FOV) = V_CM_mod_AVG_smooth_C;
                    V_CM_mod_smooth_FOV_R(FOV) = V_CM_mod_AVG_smooth_R;

                    V_rms_FOV(FOV) = V_rms_AVG; % in the CM sdr
                    V_rms_FOV_C(FOV) = V_rms_AVG_C;
                    V_rms_FOV_R(FOV) = V_rms_AVG_R;

                    V_rms_lab_CORR_FOV(FOV) = V_rms_lab_CORR_AVG;
                    V_rms_lab_CORR_FOV_C(FOV) = V_rms_lab_CORR_AVG_C;
                    V_rms_lab_CORR_FOV_R(FOV) = V_rms_lab_CORR_AVG_R;

                    psi_FOV(FOV) = psi_CORR_AVG;
                    psi_FOV_C(FOV) = psi_CORR_AVG_C;
                    psi_FOV_R(FOV) = psi_CORR_AVG_R;

            end

        end


        % Average over the FOV
        switch sample_name
            case {'C100','C100_RACDN','R100'}


                % V_CM (raw)
                VCM(sss) = mean(V_CM_mod_FOV);
                % Error bars
                VCM_err(sss) = std(V_CM_mod_FOV) / sqrt(n_fov);

                % V_CM (smooth)
                VCM_SMOOTH(sss) = mean(V_CM_mod_smooth_FOV);
                % Error bars
                VCM_SMOOTH_err(sss) = std(V_CM_mod_smooth_FOV) / sqrt(n_fov);

                % V_RMS (CM sdr)
                VRMS(sss) = mean(V_rms_FOV);
                % Error bars
                VRMS_err(sss) = std(V_rms_FOV) / sqrt(n_fov);

                % V_RMS (LAB sdr)
                VRMS_LAB(sss) = mean(V_rms_lab_CORR_FOV);
                % Error bars
                VRMS_LAB_err(sss) = std(V_rms_lab_CORR_FOV) / sqrt(n_fov);

                % PSI
                PSI(sss) = mean(psi_FOV);
                % Error bars
                PSI_err(sss) = std(psi_FOV) / sqrt(n_fov);

            otherwise % Mixed samples

                % V_CM (raw)
                VCM(sss) = mean(V_CM_mod_FOV);
                VCM_C(sss) = mean(V_CM_mod_FOV_C);
                VCM_R(sss) = mean(V_CM_mod_FOV_R);
                % Error bars
                VCM_err(sss) = std(V_CM_mod_FOV) / sqrt(n_fov);
                VCM_C_err(sss) = std(V_CM_mod_FOV_C) / sqrt(n_fov);
                VCM_R_err(sss) = std(V_CM_mod_FOV_R) / sqrt(n_fov);

                % V_CM (smooth)
                VCM_SMOOTH(sss) = mean(V_CM_mod_smooth_FOV);
                VCM_SMOOTH_C(sss) = mean(V_CM_mod_smooth_FOV_C);
                VCM_SMOOTH_R(sss) = mean(V_CM_mod_smooth_FOV_R);
                % Error bars
                VCM_SMOOTH_err(sss) = std(V_CM_mod_smooth_FOV) / sqrt(n_fov);
                VCM_SMOOTH_C_err(sss) = std(V_CM_mod_smooth_FOV_C) / sqrt(n_fov);
                VCM_SMOOTH_R_err(sss) = std(V_CM_mod_smooth_FOV_R) / sqrt(n_fov);


                % V_RMS (CM sdr)
                VRMS(sss) = mean(V_rms_FOV);
                VRMS_C(sss) = mean(V_rms_FOV_C);
                VRMS_R(sss) = mean(V_rms_FOV_R);
                % Error bars
                VRMS_err(sss) = max( std(V_rms_FOV) / sqrt(n_fov), 0.05*VRMS(sss) );
                VRMS_C_err(sss) = max( std(V_rms_FOV_C) / sqrt(n_fov) , 0.05*VRMS_C(sss) );
                VRMS_R_err(sss) = max( std(V_rms_FOV_R) / sqrt(n_fov), 0.05*VRMS_R(sss) );

                % V_RMS (LAB sdr)
                VRMS_LAB(sss) = mean(V_rms_lab_CORR_FOV);
                VRMS_LAB_C(sss) = mean(V_rms_lab_CORR_FOV_C);
                VRMS_LAB_R(sss) = mean(V_rms_lab_CORR_FOV_R);
                % Error bars
                VRMS_LAB_err(sss) = std(V_rms_lab_CORR_FOV) / sqrt(n_fov);
                VRMS_LAB_C_err(sss) = std(V_rms_lab_CORR_FOV_C) / sqrt(n_fov);
                VRMS_LAB_R_err(sss) = std(V_rms_lab_CORR_FOV_R) / sqrt(n_fov);

                % PSI
                PSI(sss) = mean(psi_FOV);
                PSI_C(sss) = mean(psi_FOV_C);
                PSI_R(sss) = mean(psi_FOV_R);
                % Error bars
                PSI_err(sss) = max( std(psi_FOV) / sqrt(n_fov), 0.05*PSI(sss) );
                PSI_C_err(sss) = max( std(psi_FOV_C) / sqrt(n_fov), 0.05*PSI_C(sss) );
                PSI_R_err(sss) = max( std(psi_FOV_R) / sqrt(n_fov), 0.05*PSI_R(sss) );

        end

        Note = 'Velocities are in pxl/frame';

        if save_results == true
            disp('Saving results ..')
            % Results for each FOV
            save([folder_SAMPLE filename_SAVE_AVG],'*_FOV_*','*_FOV','Note')
            % Final results for each sample
            save([main_folder filename_SAVE_SSS],'VC*','VRM*','PS*','sample_name_cell','Note')
        end

        toc

    end

end


%% PLOT RESULTS

if PlotResults == true

    % Load experimental phi_R
    load([main_folder filename_LOAD_STATIX],'PHI_R_n')
    % Load DYNAMIX results
    load([main_folder filename_SAVE_SSS])


    % idx NORMAL / RACDN samples
    idx_RACDN = contains(sample_name_cell,'RACDN');
    idx_NORMAL = ~idx_RACDN;
    idx_R100 = contains(sample_name_cell,'R100');


    % Convert velocities to real units (micrometers/hour)

    % V_CM (raw)
    VCM = VCM * scale / dt;
    VCM_C = VCM_C * scale / dt;
    VCM_R = VCM_R * scale / dt;
    % Error bars
    VCM_err = VCM_err * scale / dt;
    VCM_C_err  = VCM_C_err * scale / dt;
    VCM_R_err = VCM_R_err * scale / dt;

    % V_CM (smooth)
    VCM_SMOOTH = VCM_SMOOTH * scale / dt;
    VCM_SMOOTH_C = VCM_SMOOTH_C * scale / dt;
    VCM_SMOOTH_R = VCM_SMOOTH_R * scale / dt;
    % Error bars
    VCM_SMOOTH_err = VCM_SMOOTH_err * scale / dt;
    VCM_SMOOTH_C_err = VCM_SMOOTH_C_err * scale / dt;
    VCM_SMOOTH_R_err = VCM_SMOOTH_R_err * scale / dt;


    % V_RMS (CM sdr)
    VRMS = VRMS * scale / dt;
    VRMS_C = VRMS_C * scale / dt;
    VRMS_R = VRMS_R * scale / dt;
    % Error bars
    VRMS_err = VRMS_err * scale / dt;
    VRMS_C_err = VRMS_C_err * scale / dt;
    VRMS_R_err = VRMS_R_err * scale / dt;

    % V_RMS (LAB sdr)
    VRMS_LAB = VRMS_LAB * scale / dt;
    VRMS_LAB_C = VRMS_LAB_C * scale / dt;
    VRMS_LAB_R = VRMS_LAB_R * scale / dt;
    % Error bars
    VRMS_LAB_err = VRMS_LAB_err * scale / dt;
    VRMS_LAB_C_err = VRMS_LAB_C_err * scale / dt;
    VRMS_LAB_R_err = VRMS_LAB_R_err * scale / dt;


    % Create final results tables
    T_Velocities = table(sample_name_cell,PHI_R_n,VCM_SMOOTH,VCM_SMOOTH_err,VCM_SMOOTH_C,VCM_SMOOTH_C_err,VCM_SMOOTH_R,VCM_SMOOTH_R_err,VRMS,VRMS_err,VRMS_C,VRMS_C_err,VRMS_R,VRMS_R_err,VRMS_LAB,VRMS_LAB_err,VRMS_LAB_C,VRMS_LAB_C_err,VRMS_LAB_R,VRMS_LAB_R_err,PSI,PSI_err,PSI_C,PSI_C_err,PSI_R,PSI_R_err,'VariableNames',{'sample','phi_R','V_CM','V_CM_err','V_CM_C','V_CM_C_err','V_CM_R','V_CM_R_err','V_rms','V_rms_err','V_rms_C','V_rms_C_err','V_rms_R','V_rms_R_err','V_rms_LAB','V_rms_LAB_err','V_rms_LAB_C','V_rms_LAB_C_err','V_rms_LAB_R','V_rms_LAB_R_err','PSI','PSI_err','PSI_C','PSI_C_err','PSI_R','PSI_R_err'});
    writetable(T_Velocities,[main_folder '\Velocities_Psi_' num2str(DataSet) '.xlsx']);



    % Plot

    f1 = figure(1); % NORMAL

    sgtitle(['EXP ' num2str(DataSet) ' - NORMAL'])

    subplot(1,3,1) % V_CM smoothed
    errorbar(PHI_R_n(idx_NORMAL),VCM_SMOOTH(idx_NORMAL),VCM_SMOOTH_err(idx_NORMAL),'ks','MarkerFaceColor','k','MarkerSize',15)
    hold on
    errorbar(PHI_R_n(idx_NORMAL),VCM_SMOOTH_C(idx_NORMAL),VCM_SMOOTH_C_err(idx_NORMAL),'bs','MarkerFaceColor','b','MarkerSize',15)
    errorbar(PHI_R_n(idx_NORMAL),VCM_SMOOTH_R(idx_NORMAL),VCM_SMOOTH_R_err(idx_NORMAL),'rs','MarkerFaceColor','r','MarkerSize',15)
    xlabel('\phi_R [%]')
    ylabel('V_{CM} [\mum/h]')
    xlim([0 100])
    ylim([0 40])
    set(gca,'PlotBoxAspectRatio',[1,1,1],'FontSize',20,'FontWeight','bold','LineWidth',1)

    subplot(1,3,2) % V_rms in the CM sdr
    errorbar(PHI_R_n(idx_NORMAL),VRMS(idx_NORMAL),VRMS_err(idx_NORMAL),'ks','MarkerFaceColor','k','MarkerSize',15)
    hold on
    errorbar(PHI_R_n(idx_NORMAL),VRMS_C(idx_NORMAL),VRMS_C_err(idx_NORMAL),'bs','MarkerFaceColor','b','MarkerSize',15)
    errorbar(PHI_R_n(idx_NORMAL),VRMS_R(idx_NORMAL),VRMS_R_err(idx_NORMAL),'rs','MarkerFaceColor','r','MarkerSize',15)
    xlabel('\phi_R [%]')
    ylabel('V_{rms} [\mum/h]')
    xlim([0 100])
    ylim([0 40])
    set(gca,'PlotBoxAspectRatio',[1,1,1],'FontSize',20,'FontWeight','bold','LineWidth',1)


    subplot(1,3,3) % PSI
    errorbar(PHI_R_n(idx_NORMAL),PSI(idx_NORMAL),PSI_err(idx_NORMAL),'ks','MarkerFaceColor','k','MarkerSize',15)
    hold on
    errorbar(PHI_R_n(idx_NORMAL),PSI_C(idx_NORMAL),PSI_C_err(idx_NORMAL),'bs','MarkerFaceColor','b','MarkerSize',15)
    errorbar(PHI_R_n(idx_NORMAL),PSI_R(idx_NORMAL),PSI_R_err(idx_NORMAL),'rs','MarkerFaceColor','r','MarkerSize',15)
    xlabel('\phi_R [%]')
    ylabel('\psi')
    xlim([0 100])
    ylim([0 1])
    set(gca,'PlotBoxAspectRatio',[1,1,1],'FontSize',20,'FontWeight','bold','LineWidth',1)

    set(gcf,'WindowState','maximized')



    f10 = figure(10); % NORMAL

    subplot(1,3,1) % V_CM smoothed
    errorbar(PHI_R_n(idx_NORMAL),VCM_SMOOTH(idx_NORMAL),VCM_SMOOTH_err(idx_NORMAL),'ks','MarkerFaceColor','k','MarkerSize',15,'DisplayName',['EXP ' num2str(DataSet) ' - NORMAL'])
    hold on
    xlabel('\phi_R [%]')
    ylabel('V_{CM} [\mum/h]')
    xlim([0 100])
    ylim([0 40])
    set(gca,'PlotBoxAspectRatio',[1,1,1],'FontSize',20,'FontWeight','bold','LineWidth',1)

    subplot(1,3,2) % V_rms in the CM sdr
    errorbar(PHI_R_n(idx_NORMAL),VRMS(idx_NORMAL),VRMS_err(idx_NORMAL),'ks','MarkerFaceColor','k','MarkerSize',15,'DisplayName',['EXP ' num2str(DataSet) ' - NORMAL'])
    hold on
    xlabel('\phi_R [%]')
    ylabel('V_{rms} [\mum/h]')
    xlim([0 100])
    ylim([0 40])
    set(gca,'PlotBoxAspectRatio',[1,1,1],'FontSize',20,'FontWeight','bold','LineWidth',1)

    subplot(1,3,3) % PSI
    errorbar(PHI_R_n(idx_NORMAL),PSI(idx_NORMAL),PSI_err(idx_NORMAL),'ks','MarkerFaceColor','k','MarkerSize',15,'DisplayName',['EXP ' num2str(DataSet) ' - NORMAL'])
    hold on
    xlabel('\phi_R [%]')
    ylabel('\psi')
    xlim([0 100])
    ylim([0 1])
    set(gca,'PlotBoxAspectRatio',[1,1,1],'FontSize',20,'FontWeight','bold','LineWidth',1)

    set(gcf,'WindowState','maximized')



    
    if DataSet == 20240202 || DataSet == 2022 || DataSet == 20230217 || DataSet == 20230504 || DataSet == 20230612 || DataSet == 20230802 || DataSet == 20230906 || DataSet == 20240626 || DataSet == 20240703

        % RACDN

        f2 = figure(2);
        
        sgtitle(['EXP ' num2str(DataSet) ' - RACDN'])

        subplot(1,3,1) % V_CM smoothed
        errorbar(PHI_R_n(idx_RACDN),VCM_SMOOTH(idx_RACDN),VCM_SMOOTH_err(idx_RACDN),'ko','MarkerFaceColor','k','MarkerSize',15)
        hold on
        errorbar(PHI_R_n(idx_RACDN),VCM_SMOOTH_C(idx_RACDN),VCM_SMOOTH_C_err(idx_RACDN),'bo','MarkerFaceColor','b','MarkerSize',15)
        errorbar(PHI_R_n(idx_RACDN),VCM_SMOOTH_R(idx_RACDN),VCM_SMOOTH_R_err(idx_RACDN),'ro','MarkerFaceColor','r','MarkerSize',15)
        xlabel('\phi_R [%]')
        ylabel('V_{CM} [\mum/h]')
        xlim([0 100])
        ylim([0 40])
        set(gca,'PlotBoxAspectRatio',[1,1,1],'FontSize',20,'FontWeight','bold','LineWidth',1)

        subplot(1,3,2) % V_rms in the CM sdr
        errorbar(PHI_R_n(idx_RACDN),VRMS(idx_RACDN),VRMS_err(idx_RACDN),'ko','MarkerFaceColor','k','MarkerSize',15)
        hold on
        errorbar(PHI_R_n(idx_RACDN),VRMS_C(idx_RACDN),VRMS_C_err(idx_RACDN),'bo','MarkerFaceColor','b','MarkerSize',15)
        errorbar(PHI_R_n(idx_RACDN),VRMS_R(idx_RACDN),VRMS_R_err(idx_RACDN),'ro','MarkerFaceColor','r','MarkerSize',15)
        xlabel('\phi_R [%]')
        ylabel('V_{rms}^{(CM)} [\mum/h]')
        xlim([0 100])
        ylim([0 40])
        set(gca,'PlotBoxAspectRatio',[1,1,1],'FontSize',20,'FontWeight','bold','LineWidth',1)

       
        subplot(1,3,3) % PSI
        errorbar(PHI_R_n(idx_RACDN),PSI(idx_RACDN),PSI_err(idx_RACDN),'ko','MarkerFaceColor','k','MarkerSize',15)
        hold on
        errorbar(PHI_R_n(idx_RACDN),PSI_C(idx_RACDN),PSI_C_err(idx_RACDN),'bo','MarkerFaceColor','b','MarkerSize',15)
        errorbar(PHI_R_n(idx_RACDN),PSI_R(idx_RACDN),PSI_R_err(idx_RACDN),'ro','MarkerFaceColor','r','MarkerSize',15)
        xlabel('\phi_R [%]')
        ylabel('\psi')
        xlim([0 100])
        ylim([0 1])
        set(gca,'PlotBoxAspectRatio',[1,1,1],'FontSize',20,'FontWeight','bold','LineWidth',1)

        set(gcf,'WindowState','maximized')



        f10 = figure(10); % RACDN

        subplot(1,3,1) % V_CM smoothed
        errorbar(PHI_R_n(idx_RACDN),VCM_SMOOTH(idx_RACDN),VCM_SMOOTH_err(idx_RACDN),'ms','MarkerFaceColor','m','MarkerSize',15,'DisplayName',['EXP ' num2str(DataSet) ' - RACDN'])
        hold on
        xlabel('\phi_R [%]')
        ylabel('V_{CM} [\mum/h]')
        xlim([0 100])
        ylim([0 40])
        set(gca,'PlotBoxAspectRatio',[1,1,1],'FontSize',20,'FontWeight','bold','LineWidth',1)

        subplot(1,3,2) % V_rms in the CM sdr
        errorbar(PHI_R_n(idx_RACDN),VRMS(idx_RACDN),VRMS_err(idx_RACDN),'ms','MarkerFaceColor','m','MarkerSize',15,'DisplayName',['EXP ' num2str(DataSet) ' - RACDN'])
        hold on
        xlabel('\phi_R [%]')
        ylabel('V_{rms} [\mum/h]')
        xlim([0 100])
        ylim([0 40])
        set(gca,'PlotBoxAspectRatio',[1,1,1],'FontSize',20,'FontWeight','bold','LineWidth',1)

        subplot(1,3,3) % PSI
        errorbar(PHI_R_n(idx_RACDN),PSI(idx_RACDN),PSI_err(idx_RACDN),'ms','MarkerFaceColor','m','MarkerSize',15,'DisplayName',['EXP ' num2str(DataSet) ' - RACDN'])
        hold on
        xlabel('\phi_R [%]')
        ylabel('\psi')
        xlim([0 100])
        ylim([0 1])
        set(gca,'PlotBoxAspectRatio',[1,1,1],'FontSize',20,'FontWeight','bold','LineWidth',1)

        set(gcf,'WindowState','maximized')


    end

    % Save Figures
    if save_figures == true

        folder_SAVE_fig = [main_folder '\Figures_' num2str(DataSet) '\'];
        mkdir(folder_SAVE_fig);

        savefig(f1,[folder_SAVE_fig 'Vcm_Vrms_Psi.fig'])
        print(f1,[folder_SAVE_fig 'Vcm_Vrms_Psi'],'-dpng','-r0')

        savefig(f1,[folder_SAVE_fig 'Velocities_NORMAL_' num2str(DataSet) '.fig'])
        print(f1,[folder_SAVE_fig 'Velocities_NORMAL_' num2str(DataSet)],'-dpng','-r0')
        

        if DataSet == 20240202 || DataSet == 2022 || DataSet == 20230217 || DataSet == 20230504 || DataSet == 20230612 || DataSet == 20230802 || DataSet == 20230906 || DataSet == 20240703 || DataSet == 20240626 || DataSet == 20240724 || DataSet == 20240719

            savefig(f2,[folder_SAVE_fig 'Vcm_Vrms_Psi_RACDN.fig'])
            print(f2,[folder_SAVE_fig 'Vcm_Vrms_Psi_RACDN'],'-dpng','-r0')

            savefig(f2,[folder_SAVE_fig 'Velocities_RACDN_' num2str(DataSet) '.fig'])
            print(f2,[folder_SAVE_fig 'Velocities_RACDN_' num2str(DataSet)],'-dpng','-r0')

            savefig(f10,[folder_SAVE_fig 'Velocities_COMPARISON_' num2str(DataSet) '.fig'])
            print(f10,[folder_SAVE_fig 'Velocities_COMPARISON_' num2str(DataSet)],'-dpng','-r0')

            savefig(f100,[folder_SAVE_fig 'Velocities_COMPARISON_NORM_' num2str(DataSet) '.fig'])
            print(f100,[folder_SAVE_fig 'Velocities_COMPARISON_NORM_' num2str(DataSet)],'-dpng','-r0')
        end

    end

end