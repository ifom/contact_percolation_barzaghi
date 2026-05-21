function cell_trajectories_Corr = CorrectTrajectories(cell_trajectories, V_fake)

arguments
    cell_trajectories cell
    V_fake (:,2) double
end

V_fake_x = V_fake(:,1);
V_fake_y = V_fake(:,2);

V_fake2_x = cumsum(V_fake(:,1));
V_fake2_y = cumsum(V_fake(:,2));

for k = 1 : length(cell_trajectories)

    cell_trajectories{k}.X_new(1) = cell_trajectories{k}.X(1);
    cell_trajectories{k}.Y_new(1) = cell_trajectories{k}.Y(1);

    cell_trajectories{k}.dX_new(1) = cell_trajectories{k}.dX(1);
    cell_trajectories{k}.dY_new(1) = cell_trajectories{k}.dY(1);

    idx_frame = cell_trajectories{k}.Frame(2:end);

    % Correct positions
    cell_trajectories{k}.X_new(2:end) = cell_trajectories{k}.X(2:end) - V_fake2_x(idx_frame-1);
    cell_trajectories{k}.Y_new(2:end) = cell_trajectories{k}.Y(2:end) - V_fake2_y(idx_frame-1);

    % Correct instant velocities
    cell_trajectories{k}.dX_new(2:end) = cell_trajectories{k}.dX(2:end) - V_fake_x(idx_frame);
    cell_trajectories{k}.dY_new(2:end) = cell_trajectories{k}.dY(2:end) - V_fake_y(idx_frame);

    cell_trajectories{k}.dX_new(end) = NaN;
    cell_trajectories{k}.dY_new(end) = NaN;

end


cell_trajectories_Corr = cell_trajectories;
