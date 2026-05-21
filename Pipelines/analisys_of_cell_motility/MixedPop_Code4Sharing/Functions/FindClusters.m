function [GGG,subGGG_C,subGGG_R,bins_C,binsizes_C,bins_R, binsizes_R] = FindClusters(ConnList,tab_AllCells)

% Finds the graph of all cell nuclei (with subgraphs for the two CellTypes) and all clusters
% starting from ConnList which is the connectivity list of the Delaunay Triangulation

% Make graph of all cells
TT = [ConnList(:,1) ConnList(:,2); ConnList(:,2) ConnList(:,3); ConnList(:,3) ConnList(:,1)];
GGG = graph(TT(:,1),TT(:,2)); % G = graph(s,t) specifies graph edges (s,t) in node pairs. s and t are indices of connected nodes.
nnn = 1 : numnodes(GGG);
nn = num2str(nnn');
NodeNames = cellstr(strtrim(nn)); % names correspond to node idx in the main graph

% Remake graph giving a name to each node
GGG = graph(TT(:,1),TT(:,2),[],NodeNames);

% Remove double edges
GGG = simplify(GGG);

% Find subgraphs of CTRL and RAB5 cells ("subgraph" preserves node names, not their order)
subGGG_C = subgraph(GGG,tab_AllCells.CellType==1 & tab_AllCells.idx_Bulk==1);
subGGG_R = subgraph(GGG,tab_AllCells.CellType==2 & tab_AllCells.idx_Bulk==1);

% Find connected components (clusters) and their size
[bins_C, binsizes_C] = conncomp(subGGG_C);
[bins_R, binsizes_R] = conncomp(subGGG_R);
