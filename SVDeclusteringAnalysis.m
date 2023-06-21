% written by [JW 2022]

clear; clc

% This code does the following three calculations:
% 1 - Finds the distance between SVs and the closest point on the AZ
% 2 - Finds the distance between SVs and their nearest neighbor SV
% 3 - Bins data from #1 into binsz (can be set below)

% User Variables: ---------------------------------------------------------
svx = 6; % the column number containing the x coordinates for the SVs
svy = 7; % the column number containing the y coordinates for the SVs
binsz = 50; % what value in nms you'd like to bin your data
% -------------------------------------------------------------------------

% -------------------------------------------------------------------------
% Code starts here: 
fpath = uigetdir('','Choose a folder:'); % Opens Windows explorer to grab folder of interest
fpath = [fpath,filesep]; 
paths = {[fpath,'Control',filesep],[fpath,'High',filesep]}; % separates by condition
% ,[fpath,'Low',filesep]
for f=1:numel(paths) % Iterate through conditions
    % 1 is Control
    % 2 is Low concentration
    % 3 is High concentration
    files = dir([paths{f},'*.txt']); 
    q=1;
    dist = {};
    nneighbor = {};
    bindist = {};
    binedges = {};
    for p=1:2:numel(files) % iterates through .txt files in a given condition
        SV = []; % clear previous SV matrix
        AZ = []; % clear previous AZ matrix
        SV = readmatrix([paths{f},files(p+1).name]); % reads in SV coordinate matrix
        AZ = readmatrix([paths{f},files(p).name]);   % reads in AZ coordinate matrix
        
        % If you have unnecessary letters in the files this will delete them.
        [SVNaN,~] = find(isnan(SV));
        SV(SVNaN,:) = [];
        [AZNaN,~] = find(isnan(AZ));
        AZ(AZNaN,:) = [];

        d = [];dmin = []; % clear previous distance matrix
        n = [];nmin = []; % clear previous distance matrix
        bindmin = []; % binned counts
        edges = [];   % Left edge of binned regions
        for i = 1:1:length(SV(:,1)) % loops through SV coordinates
            x1 = SV(i,svx);
            y1 = SV(i,svy);
            for j = 1:length(AZ(:,1)) % loops through AZ coordinates
                x2 = AZ(j,1);
                y2 = AZ(j,2);
                d(j) = sqrt((x1-x2).^2+(y1-y2).^2);
            end
            dmin(i) = min(d); % finds shortest distance to AZ
            for k = 1:1:length(SV(:,1)) % loop through SV coordinates
                x3 = SV(k,svx);
                y3 = SV(k,svy);
                n(k) = sqrt((x1-x3).^2+(y1-y3).^2);
            end
            nmin(i) = min(n(n>0)); % finds distance to nearest neighbor SV
        end
        % Bin SV distance to AZ data
        [bindmin,edges] = histcounts(dmin,0:binsz:max(dmin)+binsz);
        edges(end)=[];
        % Convert important data into cell arrays:
        dist(q) = {dmin'};       % Save distance to AZ in cell array for each synapse 
        nneighbor(q) = {nmin'};  % Save distance to nearest neighbor in cell array for each synapse
        bindist(q) = {bindmin'}; % Save binned values for dist
        binedges(q) = {edges'};  % Save edges for dist
        q=q+1;
    end
    bindist = bindist';
    % Concatinate NN and SV cell arrays for saving:
    dist = vertcat(dist{:});
    nneighbor = vertcat(nneighbor{:});
    % Save data:
    data = {dist;nneighbor};
    writecell(data,[paths{f},'DistanceData.csv'])
    % Save binned data:
    writecell(bindist,[paths{f},'BinnedDisttoAZ.csv'])
end
% -------------------------------------------------------------------------