function [SqDisp_CM,weight] = GetSqDisp_CM(cell_trajectories,n_frames,dt,dCM)

% Preallocation
dX_CM = NaN(numel(cell_trajectories),n_frames);
dY_CM = NaN(numel(cell_trajectories),n_frames);
weight = NaN(numel(cell_trajectories),1);

for jjj = 1 : numel(cell_trajectories)  % trajectories
    weight(jjj) = height(cell_trajectories{jjj}); % trajectory length
    for uuu = 1 : height(cell_trajectories{jjj})-dt % 1 : n_frames-dt   % time
    dX_CM(jjj,uuu) = cell_trajectories{jjj}.X(uuu+dt) - dCM(uuu+dt,1) - cell_trajectories{jjj}.X(uuu) + dCM(uuu,1);
    dY_CM(jjj,uuu) = cell_trajectories{jjj}.Y(uuu+dt) - dCM(uuu+dt,2) - cell_trajectories{jjj}.Y(uuu) + dCM(uuu,2);

    end
end

SqDisp_CM = dX_CM.^2 + dY_CM.^2;


weight = weight / sum(weight);
