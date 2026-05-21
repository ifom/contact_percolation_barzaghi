% This script calculates the CM and RMS velocities from PIV results and plots them for a dataset
% C. Guidolin


close all
clear
clc

main_folder = 'D:\MIXED_POPULATIONS\Leonardo_2023_Analysis\MixedPop_RACDN_20230612\PIV\PhC\';


dt_PIV = 1; % frames
window = 32; % pixels

motility_peak = 144:240;%192:288;%144:240;%49:144; % frames

scale = 1342.85/1024; % um / pxl
dt = dt_PIV*5/60; % hours

nfov = 5;

files = dir([main_folder '\*.tif']);

% Find all filenames
names = cell(numel(files),1);
sample_names = cell(numel(files),1);
for gg = 1 : numel(files)
    names{gg} = files(gg).name;
    sample_names{gg} = names{gg}(1:end-9);
end

sample_names = unique(sample_names);
phi = [0 80 80 70 70 50 50 30 30 20 20 100]; % Nominal values


% Dataset 20230504
phi_exp = [0 78.6730947252809 78.7395316368063 74.3408075025222 73.5794510359851 49.5993269057908 51.7975506177464 24.4457304478182 25.1648838211218 21.0738569670689 22.2833534894069 100];
% phi_exp is the list of experimental values of phi calculated from STATIX scripts


idx_RACDN = contains(sample_names,'RACDN');
idx_NORMAL = ~idx_RACDN;

% Loop over the samples
for sss = 1 : numel(sample_names)

    for fov = 1 : nfov

        load([main_folder sample_names{sss} '_FOV' num2str(fov) '_PIV_window' num2str(window) '_dt' num2str(dt_PIV) '_res.mat'])

        % U and V are the velocity components in the lab sdr (raw)

        % Remove datapoints corresponding to first and last frame
        U = U(2:end-1,2:end-1,:);
        V = V(2:end-1,2:end-1,:);
        U(isnan(U)) = .1*randn;
        V(isnan(V)) = .1*randn;


        % Velocity of the CM (raw):
        u = squeeze(mean(mean(U))); % ('squeeze' just removes all singleton dimensions)
        v = squeeze(mean(mean(V)));

        % Velocity of the CM (smoothed):
        u_sm = smooth(u,15);
        v_sm = smooth(v,15);
        % Delta correction:
        du = u - u_sm;
        dv = v - v_sm;

        % Preallocation
        UU = U*0; Ucm = UU; nU = UU;
        VV = V*0; Vcm = VV; nV = VV;

        for kk = 1 : length(U) % time
            UU(:,:,kk) = U(:,:,kk) - du(kk); % Velocity in the lab (corrected)
            VV(:,:,kk) = V(:,:,kk) - dv(kk);
            Ucm(:,:,kk) = U(:,:,kk) - u(kk); % Velocity in the CM
            Vcm(:,:,kk) = V(:,:,kk) - v(kk);
            VU = sqrt(UU(:,:,kk).^2 + VV(:,:,kk).^2);
            nU(:,:,kk) = UU(:,:,kk) ./ VU;
            nV(:,:,kk) = VV(:,:,kk) ./ VU;
        end

        us = squeeze(mean(mean(UU)));
        vs = squeeze(mean(mean(VV)));

        h0 = sqrt(us.^2+vs.^2)./(squeeze(mean(mean(sqrt(UU.^2+VV.^2)))));

        h3 = sqrt(squeeze(mean(mean(Ucm.^2+Vcm.^2)))); % Velocity rms in the CM
        h1 = sqrt(squeeze(mean(mean(UU.^2+VV.^2)))); % Velocity rms in the lab

        h4 = sqrt(u_sm.^2 + v_sm.^2); % Velocity of the CM (smoothed)

        h = sqrt(u_sm.^2 + v_sm.^2) ./ h1; % Order Parameter (V_CM modulus / Vrms lab)

        

        figure(sss)
        sgtitle(replace(sample_names{sss},'_','-'),'FontSize',30,'FontWeight','bold')

        subplot(1,2,1) % Vcm
        plot(h4(2:end-1)*scale/dt,'-','LineWidth',1,'DisplayName',['FOV ' num2str(fov)])
        hold on
        xlim([0 600])
        % ylim([0 30])
        xlabel('frame')
        ylabel('V_{CM} [\mum/h]')
        legend show
        set(gca,'PlotBoxAspectRatio',[1,1,1],'FontSize',20,'FontWeight','bold','LineWidth',1)


        subplot(1,2,2) % Vrms
        plot(h3*scale/dt,'-','LineWidth',1,'DisplayName',['FOV ' num2str(fov)])
        hold on
        xlim([0 600])
        % ylim([0 30])
        xlabel('frame')
        ylabel('V_{rms} [\mum/h]')
        set(gca,'PlotBoxAspectRatio',[1,1,1],'FontSize',20,'FontWeight','bold','LineWidth',1)


        VRMS_fov_time(fov,:) = h3;   % Vrms
        Vcm_fov_time(fov,:) = h4; % V_CM



        OP_fov(fov,sss) = mean(h(motility_peak));      % Order Parameter
        VRMS_fov(fov,sss) = mean(h3(motility_peak));   % Vrms

        Vcm_fov(fov,sss) = mean(h4(motility_peak)); % V_CM



    end

    VRMS_time = mean(VRMS_fov_time);
    VCM_time = mean(Vcm_fov_time);

    VRMS_time_err = std(VRMS_fov_time) / sqrt(nfov);
    VCM_time_err = std(Vcm_fov_time) / sqrt(nfov);


    figure(sss)
    subplot(1,2,1)
    plot(VCM_time*scale/dt,'r-','LineWidth',4,'DisplayName','AVG')
    hold on
    xline(motility_peak(1),'k--','DisplayName','t_i')
    xline(motility_peak(end),'k--','DisplayName','t_f')
    hold on

    subplot(1,2,2)
    plot(VRMS_time*scale/dt,'r-','LineWidth',4,'DisplayName','AVG')
    hold on
    xline(motility_peak(1),'k--','DisplayName','t_i')
    xline(motility_peak(end),'k--','DisplayName','t_f')
    hold on
    

end



OP = mean(OP_fov);      OP_err = std(OP_fov)/sqrt(nfov);
VRMS = mean(VRMS_fov);  VRMS_err = std(VRMS_fov)/sqrt(nfov);
VCM = mean(Vcm_fov);    VCM_err = std(Vcm_fov)/sqrt(nfov);



figure(100)

subplot(1,2,1)
errorbar(phi_exp(idx_NORMAL),VCM(idx_NORMAL)*scale/dt,VCM_err(idx_NORMAL)*scale/dt,'ks','MarkerSize',20,'LineWidth',1.5,'DisplayName','NORMAL')
hold on
errorbar(phi_exp(idx_RACDN),VCM(idx_RACDN)*scale/dt,VCM_err(idx_RACDN)*scale/dt,'ms','MarkerSize',20,'LineWidth',1.5,'DisplayName','RACDN')
xlabel('\phi_R [%]')
ylabel('V_{CM} [\mum/h]')
xlim([0 100])
ylim([0 30])
legend show
set(gca,'PlotBoxAspectRatio',[1,1,1],'FontSize',20,'FontWeight','bold','LineWidth',1)


subplot(1,2,2)
errorbar(phi_exp(idx_NORMAL),VRMS(idx_NORMAL)*scale/dt,VRMS_err(idx_NORMAL)*scale/dt,'ks','MarkerSize',20,'LineWidth',1.5,'DisplayName','NORMAL')
hold on
errorbar(phi_exp(idx_RACDN),VRMS(idx_RACDN)*scale/dt,VRMS_err(idx_RACDN)*scale/dt,'ms','MarkerSize',20,'LineWidth',1.5,'DisplayName','RACDN')
xlabel('\phi_R [%]')
ylabel('V_{RMS} [\mum/h]')
xlim([0 100])
ylim([0 35])
legend show
set(gca,'PlotBoxAspectRatio',[1,1,1],'FontSize',20,'FontWeight','bold','LineWidth',1)


