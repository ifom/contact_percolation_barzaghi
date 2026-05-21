function [fraction_trajectories, cell_ALL_trajectories, tab_spots, ID_tracking] = GetTrajectories(tab_spots, tab_tracks)

arguments
    tab_spots table
    tab_tracks table
end

% Sort rows of tab_spots according to time for each trajectory
tab_spots = sortrows(tab_spots,{'TRACK_ID','FRAME'});

% Preallocation
tab_spots.dX = zeros(height(tab_spots),1);
tab_spots.dY = zeros(height(tab_spots),1);

% Select trajectories with no gaps and lasting at least 5 frames
ID_tracking = tab_tracks.TRACK_ID(tab_tracks.NUMBER_GAPS == 0 & tab_tracks.NUMBER_SPOTS >= 5);%tab_tracks.NUMBER_SPOTS >= 5);

% Remove discarded tracks from tab_spots
idx_PreSelection = ismember(tab_spots.TRACK_ID,ID_tracking);
tab_spots = tab_spots(idx_PreSelection,:);

fraction_trajectories = numel(ID_tracking) / numel(tab_tracks.TRACK_ID);

% Cell array that will contain all trajectories
% (Coordinates X(t), Y(t), Frame and Area stored in a table for each trajectory)
cell_ALL_trajectories = cell(numel(nonzeros(ID_tracking)),1);

disp(['n. trajectories = ' num2str(numel(ID_tracking))])

% Loop over the selected trajectories
for kkk = ID_tracking'

    idx_TR = tab_spots.TRACK_ID == kkk;

    X = tab_spots.POSITION_X( idx_TR );
    Y = tab_spots.POSITION_Y( idx_TR );

    dX = [diff(X); 0];
    dY = [diff(Y); 0];

    tab_spots.dX( idx_TR ) = dX;
    tab_spots.dY( idx_TR ) = dY;    

    Frame = tab_spots.FRAME( idx_TR );
    Area = tab_spots.AREA( idx_TR );

    tab_ALL_trajectories = table(X,Y,Frame,Area,dX,dY,'VariableNames',{'X','Y','Frame','Area','dX','dY'});

    cell_ALL_trajectories{kkk+1} = tab_ALL_trajectories; % I put kkk+1 so that if TRACK_ID=0 (kkk=0) it is registered

end

cell_ALL_trajectories = cell_ALL_trajectories(~cellfun('isempty',cell_ALL_trajectories));
