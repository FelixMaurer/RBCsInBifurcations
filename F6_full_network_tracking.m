%--------------------------------------------------------------------------
% Script Name : F6_full_network_tracking.m
% Author      : Felix Maurer
% Institution : Saarland University
% Email       : mail@felixmilanmaurer.com
% Date        : 2024
%
% Description :
%   This script performs an initial coarse cell tracking to enhance the
%   prediction together with the Lucas Kanade flow estimation in the final
%   more precise tracking.
%
% Usage :
%   - The script will work with the detected peaks from the peak detection
%   algorithm. Those will be analyzes as a point cloud.
%   - When testing make sure to choose the characteristic velocity that
%   converts a frame number (time) to a spacial coordinate (x,y,k)->(x,y,z)
%   such that the trajectories have an angle of around 45° to the
%   x-y-plane.
%   - this conversion factor is stored in the variable 'dist' as a
%   reference pixel distance
%
% Dependencies :
%
% Reference :
%   This script is associated with the publication
%   Impact of Red Blood Cell Rigidity on in vivo Flow Dynamics and Lingering in Bifurcations
%   by Rashidi et al. 2025
% License :
%   MIT
%% settings
% testing the code and settings
flag_test = false;
% optionally select one ROI to test
testROI = 0; % zero means take all ROIs
% ploting, useful for debugging / parameter adjustment
flag_plot = false;
% optionally skip already analyzed files
flag_overwrite = true;
if flag_test
    % when testing also enable plotting
    flag_plot = true;
end
%% source
addpath('src');
addpath('src\SLMtools\');
%% File Loop
clc;
maskName = 'Mask.png';
rootDir = char(readlines('directory.txt'));
dist = 2; % increase to make curves steeper / increase angle, should be less then 45°
cellTypes = {'Healthy_RBCs','Rigid_RBCs'};
for minDistance = 20
    for idxType = 1:2
        cFolder = [rootDir '\' cellTypes{idxType}];
        filelist = dir(fullfile(cFolder, '**\*.mj2'));  
        % take the whole network as region of interest
        wallImg = imread([rootDir,'\',maskName]);
        if size(size(wallImg),2)>2
            wallImg = rgb2gray(wallImg);
        end
        wallBW = imbinarize(wallImg);
        ROIbounds = bwboundaries(wallBW); % one boundary for each disconnected branch

        roiIdxArray = 1:length(ROIbounds);
        if testROI ~= 0
            roiIdxArray = testROI;
        end
        for ROIidx = roiIdxArray
            bdy = ROIbounds{ROIidx};
            % peak file loop
            tic;
            for idxFile = 1:length(filelist)
                fileFolder = filelist(idxFile).folder;
                fileName = filelist(idxFile).name;
                filePath = [fileFolder '\' fileName];
                % output file
                outPath = [filePath(1:end-4) '_ROI_' num2str(ROIidx) '_traj.mat'];
                if isfile(outPath)
                    fprintf('info: exists --> %s\n',outPath);
                end
                if contains(fileFolder,cellTypes)>0 && (~isfile(outPath) || flag_overwrite)
                    % terminal output
                    fprintf('Working on --> %s\n',filePath);
                    %% load alignment
                    alignVec = [0 0];
                    try
                        load([filePath(1:end-4) '_alignVec.mat']);
                    catch
                        fprintf('loading of alignment for %s failed.\n',filePath(1:end-4))
                    end
                    %% load peaks
                    load([filePath(1:end-4) '_peaks.mat']);
                    %% optional: take out points
                    % if strcmp(fileNAME,'Experiment-4458.mj2')
                    %     allpoints = allpoints(330:end);
                    % end
                    %% convert to point cloud
                    FrameNum = length(allpoints);
                    idx0 = 1;
                    cloudPoints = [];
                    for idx = idx0:FrameNum
                        points = allpoints(idx).peak;
                        if ~isempty(points)
                            cloudPoints = [cloudPoints; points, idx * dist * ones(size(points,1),1)];
                        end
                    end
                    %% denoise
                    close all;
                    FLAGdenoise = 0;
                    ptCloud = pointCloud(cloudPoints);
                    if flag_plot
                        figure;
                        pcshow(ptCloud);
                        title('Original Data');
                    end
                    if FLAGdenoise
                        ptCloudB = pcdenoise(ptCloud,'NumNeighbors',10,'Threshold',.1);
                        if flag_plot
                            figure;
                            pcshow(ptCloudB);
                            title('Denoised Data');
                        end
                    else
                        ptCloudB = ptCloud;
                    end
                    ptCloud = ptCloudB;
                    %% ROI point selection
                    locations = ptCloud.Location;
                    % plasma alignment

                    X = locations(:,1)+alignVec(:,1);
                    Y = locations(:,2)+alignVec(:,2);
                    Z = locations(:,3);

                    % ROI conditions
                    COND = inpolygon(X,Y,bdy(:,1),bdy(:,2));
                    X = X(COND);
                    Y = Y(COND);
                    Z = Z(COND);
                    ptCloud = pointCloud([X,Y,Z]);
                    if flag_plot
                        figure;
                        pcshow(ptCloud);
                        title('ROI Data');
                    end
                    %% remove agglomerates
                    removeAgglo = 0;
                    if removeAgglo
                        % compute the local point-wise density
                        points = ptCloud.Location;
                        R = 5; % inclusive radius
                        pDens = zeros(1,size(points,1));
                        for idx=1:size(points,1)
                            pDist = sqrt( sum( (points-points(idx,:)).^2 ,2) );
                            nNeighbors = length( find(pDist<=R) );
                            pDens(idx) = nNeighbors;%/(4*pi*R.^3/3);
                        end
                        nAggloIdc = find(pDens < 30);
                        points = points(nAggloIdc,:);
                        ptCloud = pointCloud(points);
                    end
                    %% distance segmentation
                    [labels,numClusters] = pcsegdist(ptCloud,minDistance);
                    fprintf('# clusters: %d\n',numClusters);
                    %% extract clusters
                    if size(ptCloud.Location,1) > 1
                        location = ptCloud.Location;
                        [labels,sort_idc] = sort(labels);
                        location = location(sort_idc,:);
                        [~,idc] = unique(labels);
                        % segmentation into clusters
                        clear clu
                        for idx = 1:length(idc)-1
                            clu(idx).points = location(idc(idx):(idc(idx+1)-1),:);
                        end
                        if length(idc) == 1
                            clu.points = location;
                        end
                        % % filter by size
                        min_cluster_size = 3;
                        cluster_sizes = diff(idc);
                        clu = clu(cluster_sizes > min_cluster_size);
                        %% draw trajectories
                        trajpoints = [];
                        labels = [];
                        for idx = 1:length(clu)
                            cpoints = clu(idx).points;
                            trajpoints = [trajpoints; cpoints];
                            labels = [labels, idx*ones(1,size(cpoints,1))];
                        end
                        numClusters = length(clu);
                        if numClusters > 0
                            if flag_plot
                                figure;
                                pcshow(trajpoints,labels)
                                colormap(lines(numClusters))
                                title('Pre Filtered Cluster')
                            end
                            %% spline fit filter
                            clear filtclu; k = 1;
                            FLAGdispspline = 0;
                            for cluidx = 1:length(clu)
                                cpoints = clu(cluidx).points;
                                x = cpoints(:,1);
                                y = cpoints(:,2);
                                z = cpoints(:,3);
                                % data indices
                                didc = (1:length(x))';
                                % spline indices
                                Nspline = 1000;
                                sidc = linspace(1,length(x),Nspline);
                                % outlier filter
                                fpoints = zeros(size(cpoints));
                                fpoints(:,1) = hampel(x);
                                fpoints(:,2) = hampel(y);
                                fpoints(:,3) = hampel(z);
                                %knotarray = [min(xidc)-(10:10:20), linspace(min(xidc),max(xidc),4), max(xidc)+(10:10:20)];
                                knotnumber = max([2,round(length(x)/10)]);
                                knotarray = linspace(min(didc),max(didc),knotnumber);
                                spoints = zeros(Nspline,3);
                                for idx = 1:3
                                    slm = slmengine(didc,fpoints(:,idx),'plot','off','knots',knotarray);
                                    spoints(:,idx) = slmeval(sidc,slm);
                                end
                                if FLAGdispspline
                                    scatter3(cpoints(:,1),cpoints(:,2),cpoints(:,3),'ob');
                                    hold on
                                    scatter3(spoints(:,1),spoints(:,2),spoints(:,3),'.r');
                                    pause(0.5)
                                    hold off
                                end
                                %% spline distance filter
                                distthres = minDistance*0.7;
                                perc = 20;
                                minnorm = zeros(1,length(x));
                                for cidx = 1:length(x)
                                    cpoint = cpoints(cidx,:);
                                    norms = zeros(Nspline,1);
                                    for idx = 1:3
                                        norms = norms + (cpoint(idx)-spoints(:,idx)).^2;
                                    end
                                    norms = sqrt(norms);
                                    minnorm(cidx) = min(norms);
                                end
                                COND = minnorm > distthres;
                                if sum(COND)/length(COND) < perc/100
                                    % valid trajectory
                                    % save filtered trajectory
                                    savepoints = cpoints(~COND,:);
                                    minnorm = minnorm(~COND);
                                    % only take nearest if there are multiple points
                                    [~,uniqidc] = unique(savepoints(:,3));
                                    diffidc = diff(uniqidc);
                                    multptsidc = find(diffidc > 1);
                                    singleptsidc = find(diffidc == 1);
                                    chosenPoints = zeros(size(uniqidc,1),3);
                                    chosenPoints(singleptsidc,:) = savepoints(uniqidc(singleptsidc),:);
                                    chosenPoints(end,:) = savepoints(end,:);
                                    for idx = 1:length(multptsidc)
                                        interval = uniqidc(multptsidc(idx)):uniqidc(multptsidc(idx)+1)-1;
                                        chosen = find(minnorm(interval)==min(minnorm(interval)),1,'first');
                                        chosenPoints(multptsidc(idx),:) = savepoints(interval(chosen),:);
                                    end
                                    filtclu(k).points = chosenPoints;
                                    k = k+1;
                                end
                            end
                            %% filter points
                            trajpoints = [];
                            labels = [];
                            lengths = zeros(1,length(filtclu));
                            for idx = 1:length(filtclu)
                                cpoints = filtclu(idx).points;
                                trajpoints = [trajpoints; cpoints];
                                labels = [labels, idx*ones(1,size(cpoints,1))];
                                x = cpoints(:,1);
                                y = cpoints(:,2);
                                z = cpoints(:,3);
                                lengths(idx) = (max(x)-min(x))*(max(y)-min(y))*(max(z)-min(z));
                            end
                            size_filter = 0;
                            length_filter = 1;
                            if size_filter
                                % filter by size
                                min_cluster_size = 10;
                                [~,idc] = unique(labels);
                                cluster_sizes = diff(idc);
                                filtclu = filtclu(cluster_sizes > min_cluster_size);

                                % convert again
                                trajpoints = [];
                                labels = [];
                                for idx = 1:length(filtclu)
                                    cpoints = filtclu(idx).points;
                                    trajpoints = [trajpoints; cpoints];
                                    labels = [labels, idx*ones(1,size(cpoints,1))];
                                end
                            end

                            if length_filter
                                % filter by length
                                length_thres = 15000;
                                filtclu = filtclu(lengths > length_thres);
                                % convert again
                                trajpoints = [];
                                labels = [];
                                for idx = 1:length(filtclu)
                                    cpoints = filtclu(idx).points;
                                    trajpoints = [trajpoints; cpoints];
                                    labels = [labels, idx*ones(1,size(cpoints,1))];
                                end
                            end
                            numClusters = length(filtclu);
                            if ~isempty(ptCloud.Location) && numClusters > 0
                                % plot
                                if flag_plot
                                    figure;
                                    pcshow(trajpoints,labels)
                                    colormap(lines(numClusters))
                                    title('spline filter')
                                end
                                traj = filtclu;
                                % take out dist factor
                                for trajIdx = 1:length(traj)
                                    cPoints = traj(trajIdx).points;
                                    cPoints(:,3) = cPoints(:,3)/dist;
                                end
                                save(outPath,'traj');
                            end
                            if flag_test
                                return;
                            end
                        end
                    end
                end
                if flag_test
                    return;
                end
            end

        end
    end
    % merge Trajectories
    MergeNetworkTrajectories;
end
return
%% debugging: plot ROIs
%close all; figure;
hold on
for idxROI = 1:length(ROIbounds)
    pnts = ROIbounds{idxROI};
    plot(pnts(:,2),pnts(:,1),'-','color',[0.8,0.9,0.9])
end
hold off