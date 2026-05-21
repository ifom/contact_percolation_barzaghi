function PIV_Labber_CG(filename,folderin,folderout,intarea2,dt,display)

cudi = cd;
cd(folderin);

% Analyze the image with piv_FFTmulti
%disp('Performing PIV analysis with deforming windows and 4 passes...')
% Standard PIV Settings
s = cell(15,2); % To make it more readable, let's create a "settings table"
% Parameter                         % Setting           % Options
s{1,1}= 'Int. area 1';              s{1,2}=64;          % window size of first pass
s{2,1}= 'Step size 1';              s{2,2}=32;          % step of first pass
s{3,1}= 'Subpix. finder';           s{3,2}=1;           % 1 = 3point Gauss, 2 = 2D Gauss
s{4,1}= 'Mask';                     s{4,2}=[];          % If needed, generate via: imagesc(image); [temp,Mask{1,1},Mask{1,2}]=roipoly;
s{5,1}= 'ROI';                      s{5,2}=[];          % Region of interest: [x,y,width,height] in pixels, may be left empty
s{6,1}= 'Nr. of passes';            s{6,2}=2;           % 1-4 nr. of passes
s{7,1}= 'Int. area 2';              s{7,2}=intarea2;    % second pass window size
s{8,1}= 'Int. area 3';              s{8,2}=32;          % third pass window size
s{9,1}= 'Int. area 4';              s{9,2}=32;          % fourth pass window size
s{10,1}='Window deformation';       s{10,2}='*spline';  % '*spline' is more accurate, but slower
s{11,1}='Repeated Correlation';     s{11,2}=0;          % 0 or 1 : Repeat the correlation four times and multiply the correlation matrices.
s{12,1}='Disable Autocorrelation';  s{12,2}=0;          % 0 or 1 : Disable Autocorrelation in the first pass. 
s{13,1}='Correlation style';        s{13,2}=0;          % 0 or 1 : Use circular correlation (0) or linear correlation (1).
s{14,1}='Repeat last pass';         s{14,2}=0;          % 0 or 1 : Repeat the last pass of a multipass analyis
s{15,1}='Last pass quality slope';  s{15,2}=0.025;      % Repetitions of last pass will stop when the average difference to the previous pass is less than this number.

% Standard image preprocessing settings
p = cell(8,1);
% Parameter                      % Setting          % Options
p{1,1}= 'ROI';                   p{1,2}=s{5,2};     % same as in PIV settings
p{2,1}= 'CLAHE';                 p{2,2}=1;          % 1 = enable CLAHE (contrast enhancement), 0 = disable
p{3,1}= 'CLAHE size';            p{3,2}=50;         % CLAHE window size
p{4,1}= 'Highpass';              p{4,2}=0;          % 1 = enable highpass, 0 = disable
p{5,1}= 'Highpass size';         p{5,2}=15;         % highpass size
p{6,1}= 'Clipping';              p{6,2}=0;          % 1 = enable clipping, 0 = disable
p{7,1}= 'Wiener';                p{7,2}=0;          % 1 = enable Wiener2 adaptive denaoise filter, 0 = disable
p{8,1}= 'Wiener size';           p{8,2}=3;          % Wiener2 window size
p{9,1}= 'Minimum intensity';     p{9,2}=0.0;        % Minimum intensity of input image (0 = no change) 
p{10,1}='Maximum intensity';     p{10,2}=1.0;       % Maximum intensity on input image (1 = no change)



%addpath 'C:\Users\Utente\Documents\MATLAB\PIVlab'

T=length(imfinfo(filename));
time=1:1:(T-dt);
length_time=length(time);
counter=0;
for k0=time
    counter=counter+1;
    %k0
    %filename
    %tic
    A = imread(filename,k0);
    B = imread(filename,k0+dt);
    clearvars -except A B k0 filename U V counter time length_time s p folderin folderout cudi intarea2 dt display

    %fprintf('\n\n');
    % PIV analysis:

    image1 = PIVlab_preproc (A,p{1,2},p{2,2},p{3,2},p{4,2},p{5,2},p{6,2},p{7,2},p{8,2},p{9,2},p{10,2}); %preprocess images
    image2 = PIVlab_preproc (B,p{1,2},p{2,2},p{3,2},p{4,2},p{5,2},p{6,2},p{7,2},p{8,2},p{9,2},p{10,2});
    % start timer for PIV analysis only
    [x,y,u,v,~,~,~] = piv_FFTmulti(image1,image2,s{1,2},s{2,2},s{3,2},s{4,2},s{5,2},s{6,2},s{7,2},s{8,2},s{9,2},s{10,2},s{11,2},s{12,2},s{13,2},0,s{14,2},s{15,2});
    %toc

    if counter == 1
        [a,b]=size(x);
        U=zeros(a,b,length_time);
        V=zeros(a,b,length_time);
    end

    clearvars -except x y u v typevector image1 image2 A B filename k0 U V counter time length_time s p folderin folderout cudi  intarea2 dt display
%     toc

    U(:,:,counter) = u;
    V(:,:,counter) = v;

    if display == true
        hold off
        image((double(image1)+double(image2))/2);
        colormap('gray');
        hold on
        quiver(x,y,u,v,'g','AutoScaleFactor', 1.5);
        %hold off;
        axis image;
        set(gca,'xtick',[],'ytick',[])
        title('Vector map of PIV analysis')
        drawnow
        whos u
    end

end

cd(folderout);
disp('Saving results ..')
save([filename(1:(end-4)) '_PIV_window' num2str(intarea2) '_dt' num2str(dt) '_res.mat'],'U','V','x','y','s','p','time');
cd(cudi)

% return
% figure;imagesc(sqrt(u_real_reduced.^2+v_real_reduced.^2));title('Real displacement magnitude');
% figure;imagesc(sqrt(u.^2+v.^2));title('Calculated displacement magnitude');
% 
% figure;scatter(reshape(u_real_reduced,size(x,1)*size(x,2),1),reshape(u,size(x,1)*size(x,2),1),'g.')% plots real vs calculated u displacements
% xlabel('Real displacement in x-direction [px]');ylabel('Measured displacement in x-direction [px]');title('Real vs. calculated displacements in x-direction')
% figure;scatter(reshape(v_real_reduced,size(x,1)*size(x,2),1),reshape(v,size(x,1)*size(x,2),1),'b.')% plots real vs calculated v displacements
% xlabel('Real displacement in y-direction [px]');ylabel('Measured displacement in y-direction [px]');title('Real vs. calculated displacements in y-direction')
% 
% fprintf('\n\n');
% disp(['Accuracy tests finished. Elapsed time: ' num2str(elapsedtime) ' seconds.'])
% disp([ 'Mean (n = ' num2str(numel(x)) ') error u displacement: ' num2str(abs(mean2(u-u_real_reduced))) ' +- ' num2str(std2(u-u_real_reduced)) ' px'])
% disp([ 'Mean (n = ' num2str(numel(x)) ') error v displacement: ' num2str(abs(mean2(v-v_real_reduced))) ' +- ' num2str(std2(v-v_real_reduced)) ' px'])
% disp(['See figures for the detailed results.'])
% clear i j typevector
 
