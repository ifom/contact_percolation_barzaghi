function [cell_TAB, cell_Rnm, cell_VnVm] = GetVelocityCorrelation(cell_ALL_traj,n_frames)

cell_TAB = cell(n_frames-1,1);
cell_Rnm = cell(n_frames-1,1);
cell_VnVm = cell(n_frames-1,1);


for jj = 1 : n_frames-1

    TAB = NaN(numel(cell_ALL_traj),4); % TAB contains X, Y, Vx, Vy as columns in this order

    % This for loop extracts the matrix TAB containing the positions and
    % instant velocities for each particle in one frame from the list of trajectories
    for kk = 1 : length(cell_ALL_traj)
        
        idx_kk = cell_ALL_traj{kk}.Frame == jj;

        if any(idx_kk)

            % If I want to consider the raw trajectories
%             TAB(kk,1) = cell_ALL_traj{kk}.X(idx_kk);
%             TAB(kk,2) = cell_ALL_traj{kk}.Y(idx_kk);
%             TAB(kk,3) = cell_ALL_traj{kk}.dX(idx_kk);
%             TAB(kk,4) = cell_ALL_traj{kk}.dY(idx_kk);

            % If I want to consider the corrected trajectories in the lab sdr
            TAB(kk,1) = cell_ALL_traj{kk}.X_new(idx_kk);
            TAB(kk,2) = cell_ALL_traj{kk}.Y_new(idx_kk);
            TAB(kk,3) = cell_ALL_traj{kk}.dX_new(idx_kk);
            TAB(kk,4) = cell_ALL_traj{kk}.dY_new(idx_kk);
            
        end

    end

    % Switch to the CM sdr
%     TAB = TAB - mean(TAB,'omitnan');

    cell_TAB{jj} = TAB;

    XX = TAB(:,1) * ones(1,length(TAB)); % this is the column of X coord repeated n times (n = n_particles)
    YY = TAB(:,2) * ones(1,length(TAB)); % this is the column of Y coord repeated n times (n = n_particles)
    
    Rnm = sqrt( (XX - XX').^2 + (YY - YY').^2 ); % this is the distance matrix

    VnVm_X = TAB(:,3) * TAB(:,3)';
    VnVm_Y = TAB(:,4) * TAB(:,4)';

    VnVm =  VnVm_X + VnVm_Y;

    cell_Rnm{jj} = Rnm;
    cell_VnVm{jj} = VnVm;

end
    