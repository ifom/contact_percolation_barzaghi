% This script allows visualizing the movie with cryptic lamellipodia orientation together with PIV velocity field
% Plus it quantifies the alignment of lamellipodia with the local velocity field

clear
close all
clc
tic

ShowMovie = false;
EXP = 20260121;
cell_type = 'KD';
correct_CM_velocity = true;

switch EXP
    case {20260113,20250723,20250730}
        winPIV = 128;
        dtPIV = 2; % If dtPIV = 2 the velocity is taken in the previous frame (as difference between next and previous)
        RR = 150; % blob radius in pxl
    case {20260114,20260120,20260121,20251008}
        winPIV = 64;
        dtPIV = 2; % If dtPIV = 2 the velocity is taken in the previous frame (as difference between next and previous)
        RR = 150/2; % blob radius in pxl
end


main_folder = ['G:\Il mio Drive\UNIMI\Leonardo_CRYPTIC_LAMELLIPODIA\Lamellipodia_per_Chiara_2026\ECAD\' cell_type '\' num2str(EXP) '\'];


% Find filenames
DS = imageDatastore([main_folder '\Movies\Fluo']);
SAMPLE_NAMES_Fluo = extractBetween(DS.Files,[main_folder 'Movies\Fluo\'],'.tif');
DS = imageDatastore([main_folder '\Movies\PhC_PIV']);
SAMPLE_NAMES_PhC = extractBetween(DS.Files,[main_folder 'Movies\PhC_PIV\'],'.tif');
clear DS


% Loop over the different samples (positions):
for sss = 1 : numel(SAMPLE_NAMES_PhC)

    SAMPLE_PhC = SAMPLE_NAMES_PhC{sss};
    SAMPLE_Fluo = SAMPLE_NAMES_Fluo{sss};
    
    % Load table with arrows information
    TAB = readtable([main_folder 'ROI\' SAMPLE_Fluo '_roi.csv']);


    arr_Length = TAB.Length;
    arr_Frame = TAB.Slice;
    arr_Angle = TAB.Angle; % Angle from -180 to +180 degrees, it is the angle between the arrow axis and a horizontal line
    arr_X = TAB.X;
    arr_Y = TAB.Y;

    % Load PIV results
    cd([main_folder 'Movies\PhC_PIV\'])
    load([SAMPLE_PhC '_PIV_window' num2str(winPIV) '_dt' num2str(dtPIV) '_res.mat'],'x','y','U','V')

    % Some of the experimental movies contain stage oscillations that result in global shifts
    % that need to be corrected before quantifying the alignment
    if correct_CM_velocity == true
        % Velocity of the CM (raw):
        u = squeeze(mean(mean(U))); % ('squeeze' just removes all singleton dimensions)
        v = squeeze(mean(mean(V)));

        % Velocity of the CM (smoothed):
        u_sm = smooth(u,5);
        v_sm = smooth(v,5);
        % Delta correction:
        du = u - u_sm;
        dv = v - v_sm;

        % Correct velocities to reduce stage oscillation artefacts
        for kk = 1 : length(U) % time
            U(:,:,kk) = U(:,:,kk) - du(kk); % Velocities in the lab (corrected)
            V(:,:,kk) = V(:,:,kk) - dv(kk);
        end
    end

    % Preallocation
    ALIGNMENT = NaN(length(U),10); % There are never more than 10 arrows in the same frame


    % Loop over the frames:
    for ggg = 1 : length(U)
       
        
        % Check if in the frame ggg there is a hand drawn arrow:
        if ismember(ggg,arr_Frame)

            disp(ggg)

            if ShowMovie == true
                img = imread([SAMPLE_PhC '.tif'],ggg);

                img_fluo = imread([main_folder 'Movies\Fluo\' SAMPLE_Fluo '.tif'],ggg);

                if EXP == 20260113
                    img_fluo = imresize(img_fluo,2);
                end

                % Combine Fluo and PhC in an RGB image
                img_overlay = cat(3,rescale(uint16(img)),rescale(uint16(img)+img_fluo),rescale(uint16(img)));

                figure(1)
                imshow(img_overlay)
                hold on
                quiver(x,y,U(:,:,ggg)*10,V(:,:,ggg)*10,'k','AutoScale', 'off');
                drawnow
            end

            % Find index corresponding to the arrow in frame ggg
            idx_ggg = arr_Frame == ggg;

            theta = deg2rad(arr_Angle(idx_ggg));  % angle of the hand draw arrow
            L = arr_Length(idx_ggg); % length of the hand draw arrow

            % Arrow centroid coordinates
            x_C = arr_X(idx_ggg);
            y_C = arr_Y(idx_ggg);
            % Coordinates of the first end, I use them to plot the arrow in the right position
            xx = arr_X(idx_ggg) - L/2.*cos(theta); % I need the dot .* as sometimes there are more than one arrow in the same frame
            yy = arr_Y(idx_ggg) + L/2.*sin(theta);

            if EXP == 20260113
                xx = 2*xx;
                yy = 2*yy;
            end

            % Versor lamellipodia arrow
            V_lam_x = cos(theta);
            V_lam_y = -sin(theta); % the y-axis is directed downwards

            % Note 21/08/2021: We first evaluate the average velocity in the blob and
            % then we calculate its versor and its scalar product with the hand drawn arrow
            % as otherwise we might be too affected by small artifacts, if any inside the blob.
            

            % Loop over the number of arrows in the same frame
            for ttt = 1 : numel(x_C)    % numel(x_C) > 1 in case there are more than one arrow in the same frame

                DD = sqrt((x - xx(ttt)).^2 + (y - yy(ttt)).^2); % distance matrix

                % Find the PIV arrows within a distance RR from the lamellipodia arrow
                idx_RR = DD <= RR;

                
                if dtPIV == 2 && ggg > 1
                    ggg_PIV = ggg - 1; % This is to consider the velocity betweeen the previous and the next frame of the frame with the hand drawn error (ggg)
                else
                    ggg_PIV = ggg;
                end

                uu = U(:,:,ggg_PIV); % all x components
                vv = V(:,:,ggg_PIV); % all y components
                Vmoduli = sqrt(uu.^2 + vv.^2); % all moduli
                % restrict to the PIV arrows inside the blob (blob = circle of radius RR, centered around the hand drawn arrow position)
                uu_blob = uu(idx_RR);
                vv_blob = vv(idx_RR);
                x_blob = x(idx_RR);
                y_blob = y(idx_RR);

                % Find average PIV velocity inside the blob
                uu_blob_avg = mean(uu_blob); % x component
                vv_blob_avg = mean(vv_blob); % y component
                x_blob_avg = mean(x_blob);
                y_blob_avg = mean(y_blob);
                Vavg_modulus = sqrt(uu_blob_avg.^2 + vv_blob_avg.^2); % modulus

                if ShowMovie == true
                    quiver(x(idx_RR),y(idx_RR),uu_blob*10,vv_blob*10,'r','AutoScale','off');
                    viscircles([xx,yy],RR,'Color','r','EnhanceVisibility',false)
                    quiver(xx(ttt),yy(ttt),uu_blob_avg*10,vv_blob_avg*10,'b','LineWidth',2,'AutoScale','off');
                    drawnow
                end

                % Versor PIV mean velocity in the blob
                uuu = uu_blob_avg / Vavg_modulus;
                vvv = vv_blob_avg / Vavg_modulus;

                % Evaluate alignment = scalar product between the versors
                ALIGNMENT(ggg,ttt) = V_lam_x(ttt) * uuu + V_lam_y(ttt) * vvv;
                
            end

            ALIGNMENT_sss(:,sss) = ALIGNMENT(:); % each column contains the data of a different sample (position)

            
            if ShowMovie == true
                quiver(xx,yy,L.*cos(theta),-L.*sin(theta),'g','LineWidth',2,'AutoScale', 'off');
                drawnow                
%                 quiver(x_C,y_C,V_lam_x*5,V_lam_y*5,'b','AutoScale', 'off'); % Plot versor
            end
            
        end
              
    end

end

toc

%% Plot angle and alignment distributions

ANGLES_sss = acos(ALIGNMENT_sss);

% Generate random uniform distribution of angles
angles_rand = rand(10000,1)*pi;
alignment_rand = cos(angles_rand);

figure(100)
subplot(1,2,2)
h1 = polarhistogram(angles_rand,8,'Normalization','pdf','DisplayName','RND');
hold on
polarhistogram(ANGLES_sss(~isnan(ANGLES_sss)),'BinEdges',h1.BinEdges,'Normalization','pdf','DisplayName','EXP');
box on
set(gca,'Layer','top','FontSize',15,'FontWeight','bold','LineWidth',1)
title('Angle \theta')

subplot(1,2,1)
h2 = histogram(alignment_rand,8,'Normalization','pdf','DisplayName','RND');
hold on
histogram(ALIGNMENT_sss(~isnan(ALIGNMENT_sss)),'BinEdges',h2.BinEdges,'Normalization','pdf','DisplayName','EXP');
box on
set(gca,'PlotBoxAspectRatio',[1,1,1],'Layer','top','FontSize',15,'FontWeight','bold','LineWidth',1)
legend show
xlabel('cos(\theta)')
ylabel('pdf')

% Test Statistical Significance with Kolmogorov-Smirnov test
[h,p_value] = kstest2(ANGLES_sss(:),angles_rand,'Alpha',0.005);

n_arrows = numel(ANGLES_sss(~isnan(ANGLES_sss)));
disp(['# arrows = ' num2str(n_arrows)])

% Save results
save([main_folder 'RESULTS.mat'], 'p_value', 'n_arrows', 'ALIGNMENT_sss', 'ANGLES_sss', 'alignment_rand','angles_rand','winPIV','dtPIV','RR')

