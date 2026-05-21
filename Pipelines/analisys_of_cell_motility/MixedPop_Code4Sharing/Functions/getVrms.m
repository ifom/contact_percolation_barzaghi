function [Vrms,vCM] = getVrms(cell_traj,n_frames,img_size,dist_edge)

% Get Vrms from the trajectories

% img_size_x = img_size(1);
% img_size_y = img_size(2);

dx_i = zeros(numel(cell_traj),n_frames-1);
dy_i = zeros(numel(cell_traj),n_frames-1);
x_i = zeros(numel(cell_traj),n_frames);
y_i = zeros(numel(cell_traj),n_frames);
weight = zeros(numel(cell_traj),1);

for jjj = 1 : numel(cell_traj)
    weight(jjj) = height(cell_traj{jjj});
    for ttt = min(cell_traj{jjj}.Frame) : max(cell_traj{jjj}.Frame)   % n_frames-1

        % Like this I'm keeping ALL trajectories
        idx_frame = cell_traj{jjj}.Frame == ttt;
        dx_i(jjj,ttt) = cell_traj{jjj}.dX(idx_frame);
        dy_i(jjj,ttt) = cell_traj{jjj}.dY(idx_frame);
        x_i(jjj,ttt) = cell_traj{jjj}.X(idx_frame);
        y_i(jjj,ttt) = cell_traj{jjj}.Y(idx_frame);

        % Like this I'm discarding the whole trajectory if only one point is on the borders
%         if all(cell_traj{jjj}.X > dist_edge & cell_traj{jjj}.X < (img_size_x-dist_edge) & cell_traj{jjj}.Y > dist_edge & cell_traj{jjj}.Y < (img_size_y-dist_edge))
%             idx_frame = cell_traj{jjj}.Frame == ttt;
%             dx_i(jjj,ttt) = cell_traj{jjj}.dX(idx_frame);
%             dy_i(jjj,ttt) = cell_traj{jjj}.dY(idx_frame);
%             x_i(jjj,ttt) = cell_traj{jjj}.X(idx_frame);
%             y_i(jjj,ttt) = cell_traj{jjj}.Y(idx_frame);
%         end
        
        % Like this I'm discarding only the single point on the border
        % which makes the CM oscillate much more as there are particles that appear and disappear between consecutive frames        
%         idx_frame = cell_traj{ggg}.Frame == jjj;
%         if cell_traj{ggg}.X(idx_frame) > dist_edge && cell_traj{ggg}.X(idx_frame) < (img_size_x-dist_edge) && cell_traj{ggg}.Y(idx_frame) > dist_edge && cell_traj{ggg}.Y(idx_frame) < (img_size_y-dist_edge)
%             dx_i(ggg,jjj) = cell_traj{ggg}.dX(idx_frame);
%             dy_i(ggg,jjj) = cell_traj{ggg}.dY(idx_frame);
%             x_i(ggg,jjj) = cell_traj{ggg}.X(idx_frame);
%             y_i(ggg,jjj) = cell_traj{ggg}.Y(idx_frame);
%         end

    end
    
end

dx_i(dx_i==0)=NaN;
dy_i(dy_i==0)=NaN;
% x_i(x_i==0)=NaN;
% y_i(y_i==0)=NaN;

% Discard jumps
dx_i(abs(dx_i) > (3 * std(dx_i,0,'all','omitnan')) )= NaN;
dy_i(abs(dy_i) > (3 * std(dy_i,0,'all','omitnan')) )= NaN;

% These are a row vectors with vCM components for each frame
vCM_x = mean(dx_i,'omitnan');
vCM_y = mean(dy_i,'omitnan');

% These are a row vectors with CM components for each frame
% CM_x = mean(x_i,'omitnan');
% CM_y = mean(y_i,'omitnan');

vCM = [vCM_x', vCM_y'];
% CM = [CM_x', CM_y'];

weight = weight / sum(weight);

% Here I use implicit expansion to subtract the row vector vCM* to the matrix d*_i

Vrms = sqrt( mean( ((dx_i - vCM_x).^2 + (dy_i - vCM_y).^2),'all','omitnan') ); % average over all elements

% Vrms = sqrt( sum( weight .* mean( ((dx_i - vCM_x).^2 + (dy_i - vCM_y).^2),2,'omitnan') ) ); % weighted average