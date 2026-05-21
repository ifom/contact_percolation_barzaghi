function cell_ALL_traj_selected = SelectTrajectories(cell_ALL_traj,min_length_traj)

arguments
    cell_ALL_traj cell
    min_length_traj (1,1) double
end

% img_size_x = img_size(1);
% img_size_y = img_size(2);

cell_ALL_traj_selected = cell(numel(cell_ALL_traj),1);

% Remove trajectories shorter than min_length_traj and touching the borders (second if condition)
for ggg = 1 : numel(cell_ALL_traj)
    if height(cell_ALL_traj{ggg}) >= min_length_traj %&& all(cell_ALL_traj{ggg}.X > dist_edge & cell_ALL_traj{ggg}.X < (img_size_x-dist_edge) & cell_ALL_traj{ggg}.Y > dist_edge & cell_ALL_traj{ggg}.Y < (img_size_y-dist_edge))
        cell_ALL_traj_selected{ggg} = cell_ALL_traj{ggg};
    end
end

cell_ALL_traj_selected = cell_ALL_traj_selected(~cellfun('isempty',cell_ALL_traj_selected));
