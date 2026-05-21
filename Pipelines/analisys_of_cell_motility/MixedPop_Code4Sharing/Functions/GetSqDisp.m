function SqDisp = GetSqDisp(cell_trajectories,n_frames,dt)

arguments
    cell_trajectories cell
    n_frames (1,1) double {mustBeInteger}
    dt (1,1) double {mustBeInteger}
end

% It calculates the squared displacement for each cell for a given delay dt

% Preallocation
dX = NaN(numel(cell_trajectories),n_frames);
dY = NaN(numel(cell_trajectories),n_frames);

for lll = 1 : numel(cell_trajectories)  % trajectories
    for uuu = 1 : height(cell_trajectories{lll})-dt %n_frames-dt   % time     % (numel(cell_trajectories{lll}.X)-dt)
        dX(lll,uuu) = cell_trajectories{lll}.X(uuu + dt) - cell_trajectories{lll}.X(uuu);
        dY(lll,uuu) = cell_trajectories{lll}.Y(uuu + dt) - cell_trajectories{lll}.Y(uuu);
    end
end

SqDisp = dX.^2 + dY.^2;
