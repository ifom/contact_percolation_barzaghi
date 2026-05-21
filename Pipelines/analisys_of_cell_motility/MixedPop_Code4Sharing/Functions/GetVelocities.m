function [CM, V_CM, V_CM_smooth, V_rms_CM, V_rms_lab, V_rms_lab_CORR, V_fake] = GetVelocities(n_frames,tab_spots,SmoothingWindow)   % dist_edge,img_size)

arguments
    n_frames (1,1) double {mustBeInteger}
    tab_spots table
    SmoothingWindow (1,1) double {mustBeInteger}
end

% Note: All velocities are calculated between two consecutive frames (dt=1)
% The imported table tab_spots has already only the spots corresponding to
% trajectories selected according to their length

% Image size
% img_size_x = img_size(1); % img_size_x corresponds to POSITION_X: checked
% img_size_y = img_size(2); % img_size_y corresponds to POSITION_Y: checked

% I put the velocity at the end of each trajectory equal to NaN
% instead of 0 to avoid them when averaging
tab_spots.dX(tab_spots.dX == 0) = NaN;
tab_spots.dY(tab_spots.dY == 0) = NaN;

max_disp = 5; % pxl

% Preallocation
V_CM_x = zeros(n_frames-1,1);
V_CM_y = zeros(n_frames-1,1);
V_rms_CM = zeros(n_frames-1,1);
V_rms_lab = zeros(n_frames-1,1);
V_rms_lab_CORR = zeros(n_frames-1,1);
x_CM = zeros(n_frames,1);
y_CM = zeros(n_frames,1);

% Calculate raw V_CM and V_rms (the latter both in the CM and in the lab sdr)

for jjj = 1 : n_frames    % time

    % Select nuclei in the same frame
    idx_frame = tab_spots.FRAME == jjj;


    % Instantaneous coordinates (lab sdr)
    x_i = tab_spots.POSITION_X(idx_frame);
    y_i = tab_spots.POSITION_Y(idx_frame);

    % Instantaneous velocity components (lab sdr)
    v_x_i = tab_spots.dX(idx_frame);
    v_y_i = tab_spots.dY(idx_frame);

    % CM coordinates
    x_CM(jjj) = mean( x_i, 'omitnan' );
    y_CM(jjj) = mean( y_i, 'omitnan' );

%     % Change sdr (lab to CM)
%     x_i_cm = x_i - x_CM(jjj);
%     y_i_cm = y_i - y_CM(jjj);

    % Instantaneous CM velocity components
    V_CM_x(jjj) = mean( v_x_i, 'omitnan' );
    V_CM_y(jjj) = mean( v_y_i, 'omitnan' );
    
    % V_rms (CM sdr)
    V_rms_CM(jjj) =  mean( ((v_x_i - V_CM_x(jjj)).^2 + (v_y_i - V_CM_y(jjj)).^2), 'omitnan') ; % the sqrt is taken afterwards, after averaging over the time window
    
    % V_rms (lab sdr)
    V_rms_lab(jjj) =  mean( ((v_x_i).^2 + (v_y_i).^2), 'omitnan') ; % this will not be used anywhere as the raw velocities have fake oscillations coming from the microscope stage movements

end

% CM coordinates and raw CM velocity
CM = [x_CM, y_CM];
V_CM = [V_CM_x, V_CM_y];

% Smooth the CM velocity to reduce fake movements
V_CM_smooth_x = smooth(V_CM_x,SmoothingWindow,'moving');
V_CM_smooth_y = smooth(V_CM_y,SmoothingWindow,'moving');

V_CM_smooth = [V_CM_smooth_x, V_CM_smooth_y];

% Calculate fake velocity as the difference between the raw and the smooth V_CM
V_fake_x = V_CM_x - V_CM_smooth_x;
V_fake_y = V_CM_y - V_CM_smooth_y;

V_fake = [V_fake_x, V_fake_y];


% Calculate again the Vrms from the corrected velocities

for jjj = 1 : n_frames   % time

    % Select nuclei in the same frame
    idx_frame = tab_spots.FRAME == jjj;


    % Instantaneous velocity components - CORRECTED
    v_x_i = tab_spots.dX(idx_frame) - V_fake_x(jjj);
    v_y_i = tab_spots.dY(idx_frame) - V_fake_y(jjj);
    
    % V_rms (lab sdr) - CORRECTED
    V_rms_lab_CORR(jjj) = sqrt( mean( ((v_x_i).^2 + (v_y_i).^2), 'omitnan') ); % here we already take the sqrt
    
end
