% This script analyses the images of labelled nuclei (obtained with StarDist, in the time window of interest)
% After performing a Delaunay triangulation of the centroids, it gets the Voronoi tessellation and gets different information on cell neighbours and clusters
% All info are stored and then saved in a cell array called "cell_ALL_tables"
% This script can also be used to directly visualize the Voronoi tessellation
% C. Guidolin

clear
close all
clc

tic

DataSet = 2021; % The dataset is mostly identified either by the year (YYYY) or the full date (YYYYMMDD) in case of more than one experiment per year
% Read the function GetDatasetInfo for special cases or if in doubt

HD = 'D'; % Hard disk

save_results = false;
displayTriangulation = false;
displayColouredPattern = false;

% It uses the function GetDatasetInfo.m, FindClusters.m
addpath([HD ':\GUIDOLIN_CHIARA\MIXED_POPULATIONS\LATEST_CODE_VERSION\Functions\'])

dist_edge_percolation = 30; % pxl

% It saves:
filename_SAVE = '\MAT_STATIX_Results.mat';

% Load dataset info
[main_folder, sample_name_cell, fov, ~, ~, img_size, dist_edge, n_frames] = GetDatasetInfo(DataSet,HD);


n_samples = numel(sample_name_cell);
img_size_x = img_size(1);   img_size_y = img_size(2);

% Loop over samples
for sss = 1 : n_samples

    sample_name = sample_name_cell{sss};

    disp(['Sample: ' sample_name])



    % Loop over fields of view
    for FOV = fov

        disp(['FOV=' num2str(FOV)])

        % Folder containing subfolders of labelled images
        folder_FOV = [main_folder sample_name '\FOV_' num2str(FOV)];

        % Create image datastore for both cell types
        switch sample_name
            case {'C100','C100_RACDN'}
                img_DS_C = imageDatastore([folder_FOV '\C\*.tif']);
            case 'R100'
                img_DS_R = imageDatastore([folder_FOV '\R\*.tif']);
            otherwise
                img_DS_C = imageDatastore([folder_FOV '\C\*.tif']);
                img_DS_R = imageDatastore([folder_FOV '\R\*.tif']);
        end


        % Preallocation
        cell_ALL_tables = cell(n_frames,1);
        cell_ClusterSizes = cell(n_frames,2); % first column for Ctrl, second column for RAB5
        A_MaxCluster_sum_ctrl = zeros(n_frames,1);
        A_MaxCluster_sum_RAB5 = zeros(n_frames,1);
        phi_A_MaxCluster_ctrl = zeros(n_frames,1);
        phi_A_MaxCluster_RAB5 = zeros(n_frames,1);
        isPercolating_C = zeros(n_frames,1);
        isPercolating_R = zeros(n_frames,1);
        n_max_cluster_ctrl = zeros(n_frames,1);
        n_max_cluster_RAB5 = zeros(n_frames,1);


        % Loop over frames
        for k = 1 : n_frames

            disp(k)

            switch sample_name

                case {'C100','C100_RACDN'}

                    % CTRL cell nuclei
                    img_C = readimage(img_DS_C,k); %    imread([folder '\C\' num2str(k) '.tif']);

                    tab_C = regionprops('table',img_C,double(img_C),'Area','Perimeter','Centroid','MinorAxisLength','MajorAxisLength','MinIntensity','PixelIdxList');
                    tab_C = rmmissing(tab_C);

                    if iscell(tab_C.MinIntensity)   % isa(tab_C.MinIntensity,'cell')
                        tab_C.MinIntensity = cell2mat(tab_C.MinIntensity);
                    end

                    tab_C.CellType = ones(numel( tab_C.Area ),1);

                case 'R100'

                    % RAB5 cell nuclei
                    img_R = readimage(img_DS_R,k);

                    tab_R = regionprops('table',img_R,double(img_R),'Area','Perimeter','Centroid','MinorAxisLength','MajorAxisLength','MinIntensity','PixelIdxList');
                    tab_R = rmmissing(tab_R);

                    if iscell(tab_R.MinIntensity)   % isa(tab_R.MinIntensity,'cell')
                        tab_R.MinIntensity = cell2mat(tab_R.MinIntensity);
                    end

                    tab_R.CellType = 2 * ones(numel( tab_R.Area ),1);

                otherwise

                    % CTRL cell nuclei
                    img_C = readimage(img_DS_C,k); %    imread([folder '\C\' num2str(k) '.tif']);

                    tab_C = regionprops('table',img_C,double(img_C),'Area','Perimeter','Centroid','MinorAxisLength','MajorAxisLength','MinIntensity','PixelIdxList');
                    tab_C = rmmissing(tab_C);

                    if iscell(tab_C.MinIntensity)   % isa(tab_C.MinIntensity,'cell')
                        tab_C.MinIntensity = cell2mat(tab_C.MinIntensity);
                    end

                    tab_C.CellType = ones(numel( tab_C.Area ),1);


                    % RAB5 cell nuclei
                    img_R = readimage(img_DS_R,k);

                    tab_R = regionprops('table',img_R,double(img_R),'Area','Perimeter','Centroid','MinorAxisLength','MajorAxisLength','MinIntensity','PixelIdxList');
                    tab_R = rmmissing(tab_R);

                    if iscell(tab_R.MinIntensity)   % isa(tab_R.MinIntensity,'cell')
                        tab_R.MinIntensity = cell2mat(tab_R.MinIntensity);
                    end

                    tab_R.CellType = 2 * ones(numel( tab_R.Area ),1);

            end

            % ALL cells together
            switch sample_name
                case {'C100','C100_RACDN'}
                    tab_ALL = tab_C;
                case 'R100'
                    tab_ALL = tab_R;
                otherwise
                    tab_ALL = [tab_C; tab_R];
            end


            % Check if there are doubles between R and C
            % very unlikely but it can happen at the edges
            % e.g. see dataset 2021, C80R20, FOV 2, frame 49, nucleus (x,y)=(798.5152, 288.8788)
            if length(tab_ALL.Centroid) > length(unique(tab_ALL.Centroid,'rows','stable'))
                [~,irows] = unique(tab_ALL.Centroid,'rows','stable');
                tab_ALL = tab_ALL(irows,:);
            end

            tab_ALL.AspectRatio_ALL = tab_ALL.MajorAxisLength ./ tab_ALL.MinorAxisLength;

            x_ALL = tab_ALL.Centroid(:,1);
            y_ALL = tab_ALL.Centroid(:,2);


            tab_ALL.FrameNumber = k * ones(numel(x_ALL),1);

            % Remove nuclei on the edges: idx_Bulk = 1 if cell is not within 'dist_edge' (pxl) from the frame edge
            tab_ALL.idx_Bulk = x_ALL > dist_edge & x_ALL < (img_size_x-dist_edge) & y_ALL > dist_edge & y_ALL < (img_size_y-dist_edge);



            % Show Delaunay triangulation
            if displayTriangulation == true
                edges_C = boundarymask(img_C,4);
                skel_C = bwskel(edges_C);
                edges_R = boundarymask(img_R,4);
                skel_R = bwskel(edges_R);
                tri = delaunay(x_ALL,y_ALL);

                figure(10)
                imshowpair(skel_C,skel_R,'falsecolor','ColorChannels',[2 1 0])
                hold on
                triplot(tri,x_ALL,y_ALL,'b')
                voronoi(x_ALL,y_ALL,tri,'r')
                drawnow;
                hold off
            end

            % Delaunay triangulation and Voronoi diagram
            TRI = delaunayTriangulation(x_ALL,y_ALL);
            [V, C] = voronoiDiagram(TRI);   % Equivalent to: [V, C] = voronoin([x_ALL y_ALL]);

            T = TRI.ConnectivityList;

            % Make graph and find ALL clusters
            [GGG,subGGG_ctrl,subGGG_RAB5,bins_ctrl,binsizes_ctrl,bins_RAB5,binsizes_RAB5] = FindClusters(T,tab_ALL);

            % Find largest cluster of both cell types
            idx_max_ctrl = binsizes_ctrl(bins_ctrl) == max(binsizes_ctrl);
            max_cluster_ctrl = subgraph(subGGG_ctrl, idx_max_ctrl);
            idx_max_RAB5 = binsizes_RAB5(bins_RAB5) == max(binsizes_RAB5);
            max_cluster_RAB5 = subgraph(subGGG_RAB5, idx_max_RAB5);

            % Store cluster sizes in the cell array
            cell_ClusterSizes{k,1} = binsizes_ctrl;
            cell_ClusterSizes{k,2} = binsizes_RAB5;

            % Number of cell nuclei in the largest clusters
            [~, size_MAX_C] = conncomp(max_cluster_ctrl);
            [~, size_MAX_R] = conncomp(max_cluster_RAB5);

            % The following lines account for the case of more than one largest cluster having the same size
            if numel(size_MAX_C) > 1
                n_max_cluster_ctrl(k) = nnz(idx_max_ctrl) / length(size_MAX_C);
            else
                n_max_cluster_ctrl(k) = nnz(idx_max_ctrl); % equivalent to: max(binsizes_ctrl) or size(max_cluster_ctrl.Nodes,1);
            end

            if numel(size_MAX_R) > 1
                n_max_cluster_RAB5(k) = nnz(idx_max_RAB5) / length(size_MAX_R);
            else
                n_max_cluster_RAB5(k) = nnz(idx_max_RAB5); % equivalent to: max(binsizes_ctrl) or size(max_cluster_ctrl.Nodes,1);
            end


            % Visualisation
            if displayColouredPattern == true
                figure(2)
                imshowpair(skel_C,skel_R,'falsecolor','ColorChannels',[2 0 1])
            end


            %% ALL CELLS

            % Loop over ALL the Voronoi cells
            for kk = str2double(GGG.Nodes.Name(tab_ALL.idx_Bulk==1))'

                xxx = V(C{kk},1); % x coordinate of cell vertices
                yyy = V(C{kk},2); % y coordinate of cell vertices
                if all( xxx < Inf & yyy < Inf)
                    areas_VoronoiCells_ALL(kk) = polyarea(xxx,yyy);
                    BW_SingleCell = poly2mask(xxx,yyy,img_size_y,img_size_x); % it creates an image img_size_y-by-img_size_x
                    SingleCell_tab = regionprops(BW_SingleCell,'MinorAxisLength','MajorAxisLength');
                    AR_Voronoi(kk) = SingleCell_tab(1).MajorAxisLength / SingleCell_tab(1).MinorAxisLength;
                end
                %         ppp = polyshape(xxx,yyy);
                %         perim_VoronoiCells_ALL(kk) = perimeter(ppp);
                %         areas_VoronoiCells_ALL(kk) = area(ppp);
            end

            % Add Area and AspectRatio of Voronoi cells to the main table
            tab_ALL.area_VoronoiCell = zeros(height(tab_ALL),1);
            tab_ALL.area_VoronoiCell(tab_ALL.idx_Bulk==1) = areas_VoronoiCells_ALL(tab_ALL.idx_Bulk==1);
            tab_ALL.AspectRatio_Voronoi(tab_ALL.idx_Bulk==1) = AR_Voronoi(tab_ALL.idx_Bulk==1);

            % Visualisation
            if displayColouredPattern == true
                figure(2)
            end


            switch sample_name

                case {'C100','C100_RACDN','R100'}
                    
                    % do nothing
                    

                otherwise

                    % MAX CLUSTER CTRL

                    xxx_AllCluster_C = [];
                    yyy_AllCluster_C = [];
                    % Loop over the Voronoi cells in the max ctrl cluster
                    for ccc = str2double(subGGG_ctrl.Nodes.Name(idx_max_ctrl))'  % write "str2double(subGGG_ctrl.Nodes.Name)' " if I want to color all cell clusters % 

                        xxx_C = V(C{ccc},1);
                        yyy_C = V(C{ccc},2);
                        xxx_AllCluster_C = [xxx_AllCluster_C; xxx_C];
                        yyy_AllCluster_C = [yyy_AllCluster_C; yyy_C];
                        areas_VoronoiCells_MaxCluster_C(ccc) = polyarea(xxx_C,yyy_C);

                        % Colour CTRL cells
                        if displayColouredPattern == true
                            hold on
                            fill(xxx_C,yyy_C,[.7 .7 1],'FaceAlpha',0.7)
                        end

                    end

                    % Check whether it's percolating or not
                    isPercolating_C(k) = ( isscalar(size_MAX_C) && min(xxx_AllCluster_C) < dist_edge_percolation && max(xxx_AllCluster_C) > (img_size_x - dist_edge_percolation)) || length(size_MAX_C) == 1 && (min(yyy_AllCluster_C) < dist_edge_percolation && max(yyy_AllCluster_C) > (img_size_y - dist_edge_percolation));


                    % MAX CLUSTER RAB5

                    xxx_AllCluster_R = [];
                    yyy_AllCluster_R = [];
                    % Loop over the nuclei in the max RAB5 cluster
                    for rrr = str2double(subGGG_RAB5.Nodes.Name(idx_max_RAB5))'  % find(tab_ALL.idx_Bulk == 1)'    %   1 : numel(x_ALL)

                        xxx_R = V(C{rrr},1);
                        yyy_R = V(C{rrr},2);
                        xxx_AllCluster_R = [xxx_AllCluster_R; xxx_R];
                        yyy_AllCluster_R = [yyy_AllCluster_R; yyy_R];
                        areas_VoronoiCells_MaxCluster_R(rrr) = polyarea(xxx_R,yyy_R);

                        % Colour RAB5 cells
                        if displayColouredPattern == true
                            hold on
                            fill(xxx_R,yyy_R,'r','FaceAlpha',0.7)
                        end

                    end

                    % Check whether it's percolating or not
                    isPercolating_R(k) = ( isscalar(size_MAX_R) && min(xxx_AllCluster_R) < dist_edge_percolation && max(xxx_AllCluster_R) > (img_size_x - dist_edge_percolation)) || length(size_MAX_R) == 1 && (min(yyy_AllCluster_R) < dist_edge_percolation && max(yyy_AllCluster_R) > (img_size_y - dist_edge_percolation));


                    % Experimental number fraction
                    n_C = nnz(tab_ALL.idx_Bulk == 1 & tab_ALL.CellType == 1);
                    n_R = nnz(tab_ALL.idx_Bulk == 1 & tab_ALL.CellType == 2);
                    phi_n_C = n_C / (n_C + n_R);
                    phi_n_R = n_R / (n_C + n_R);

                    % Area max clusters
                    A_MaxCluster_sum_ctrl(k) = sum(areas_VoronoiCells_MaxCluster_C);  % sum(nonzeros(areas_VoronoiCells_ctrl_cluster));
                    A_MaxCluster_sum_RAB5(k) = sum(areas_VoronoiCells_MaxCluster_R);  % sum(nonzeros(areas_VoronoiCells_RAB5_cluster));

                    % Area fraction covered by max cluster
                    area_tot_ALL = sum(areas_VoronoiCells_ALL(tab_ALL.idx_Bulk == 1),'omitnan');
                    phi_A_MaxCluster_ctrl(k) = A_MaxCluster_sum_ctrl(k)./area_tot_ALL;
                    phi_A_MaxCluster_RAB5(k) = A_MaxCluster_sum_RAB5(k)./area_tot_ALL;

                    clear areas_VoronoiCells_MaxCluster_C areas_VoronoiCells_MaxCluster_R
                    clear rrr ccc

            end

            clear areas_VoronoiCells_ALL AR_Voronoi



            %% NEIGHBOURS AND DEMIXING PARAMETER

            % Preallocation
            n_neigh = zeros(length(TRI.Points),1);
            n_neigh_C = zeros(length(TRI.Points),1);
            n_neigh_R = zeros(length(TRI.Points),1);
            DP_ctrl = zeros(length(TRI.Points),1);
            DP_RAB5 = zeros(length(TRI.Points),1);

            for gg = 1 : length(TRI.Points) % Loop over all the nuclei in the frame
                [r, c] = find(T == gg);
                idx_neigh = unique(T(r,:));         % rmv repetitions
                idx_neigh(idx_neigh == gg) = [];    % rmv central cell
                type_center = tab_ALL.CellType(gg); % cell type of central cell
                type_neigh = tab_ALL.CellType(idx_neigh);   % cell type of neighbours
                n_neigh(gg) = numel(idx_neigh);     % count neighbours for each cell

                switch sample_name
                    case {'C100','C100_RACDN'}
                        continue;
                    case 'R100'
                        continue;
                    otherwise
                        n_neigh_C(gg) = nnz(type_neigh == 1);   % count neighbours of type 1 (CTRL) for each cell
                        n_neigh_R(gg) = nnz(type_neigh == 2);   % count neighbours of type 2 (RAB5) for each cell
                        DP_ctrl(gg) = 2*(n_neigh_C(gg)/n_neigh(gg) - phi_n_C); % Demixing parameter CTRL cells
                        DP_RAB5(gg) = 2*(n_neigh_R(gg)/n_neigh(gg) - phi_n_R); % Demixing parameter RAB5 cells
                        % Note: the correction according to the new (correct) DP definition is done afterwards
                end
            end

            % Demixing parameter: if DP>0 cells tend to cluster, if DP<0 repulsion between cells of same type
            DP = zeros(length(TRI.Points),1);
            DP(tab_ALL.CellType == 1) = DP_ctrl(tab_ALL.CellType == 1);
            DP(tab_ALL.CellType == 2) = DP_RAB5(tab_ALL.CellType == 2);
            tab_ALL.DemixingParameter = DP;
            tab_ALL.n_neighbours = n_neigh;
            tab_ALL.n_neighbours_C = n_neigh_C;
            tab_ALL.n_neighbours_R = n_neigh_R;

            clear n_neigh n_neigh_C n_neigh_R DP_ctrl DP_RAB5 DP n_C n_R
            clear area_tot_ALL phi_n_R phi_n_C

            cell_ALL_tables{k} = tab_ALL;

        end

        % Save results
        if save_results == true
            save([folder_FOV filename_SAVE])
        end

        toc

    end

end

toc
