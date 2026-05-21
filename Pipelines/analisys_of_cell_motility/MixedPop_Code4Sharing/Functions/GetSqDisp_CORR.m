function SqDisp_CORR = GetSqDisp_CORR(cell_trajectories,n_frames,dt)

arguments
    cell_trajectories cell
    n_frames (1,1) double {mustBeInteger}
    dt (1,1) double {mustBeInteger}
end

% It calculates the squared displacement for each cell for a given delay dt

% Preallocation
dX_new = NaN(numel(cell_trajectories),n_frames);
dY_new = NaN(numel(cell_trajectories),n_frames);

for lll = 1 : numel(cell_trajectories)  % trajectories
    for uuu = 1 : height(cell_trajectories{lll})-dt   % time
        dX_new(lll,uuu) = cell_trajectories{lll}.X_new(uuu+dt) - cell_trajectories{lll}.X_new(uuu);
        dY_new(lll,uuu) = cell_trajectories{lll}.Y_new(uuu+dt) - cell_trajectories{lll}.Y_new(uuu);
        
    end
end

SqDisp_CORR = dX_new.^2 + dY_new.^2;
