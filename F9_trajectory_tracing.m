%--------------------------------------------------------------------------
% Script Name : F9_trajectory_tracing.m
% Authors     : Felix Maurer, Yazdan Rashidi
% Institution : Saarland University
% Email       : mail@felixmilanmaurer.com
% Date        : 2024
%
% Description :
%   This script traces cell trajectories in the 3D detected points (x,y,z),
%   using a predictive search from the flow estimations of coarse HRBC
%   tracking and plasma optical flow estimation.
%
% Usage :
%   - The script requires previous full network tracking and Lucan Kanade
%   flow estimation
%
% Dependencies :
%
% Reference :
%   This script is associated with the publication
%   Impact of Red Blood Cell Rigidity on in vivo Flow Dynamics and Lingering in Bifurcations
%   by Rashidi et al. 2025
% License :
%   MIT
%% source
addpath('src');
%% settings
% pause between plotting frames
framePause = inf; % [sec] set to inf for no plotting, good values > 0.5 s        
% for testing, activates intermediate steps plotting
flag_test = false;
% in case tracing already saved previously to skip tracing part
flag_trace = true;
%% file loop
clc;
cellTypes = {'Healthy_RBCs','Rigid_RBCs'};
% find data directories
maskName = 'Mask.png';
rootdir = char(readlines('directory.txt'));
if flag_trace
    rootDir = rootdir;
    filelist = dir(fullfile(rootDir, '**\',maskName));  % get list of files and folders in any subfolder
    folders = unique({filelist.folder});
    filelist = dir(fullfile(rootDir, '**\*peaks.mat'));  % get list of files and folders in any subfolder
    % find the files in the analyzed folders
    selIdc = zeros(1,length(filelist),'logical');
    for fileIdx = 1:length(filelist)
        fileFolder = filelist(fileIdx).folder;
        fileName = filelist(fileIdx).name;
        filePath = [fileFolder '\' fileName];
        if sum(strcmp(fileFolder(1:end-length('Healthy_RBCs')-1),folders))>0 || sum(strcmp(fileFolder(1:end-length('Rigid_RBCs')-1),folders))>0
            selIdc(fileIdx) = 1;
        end
    end

    % apply selection
    filelist = filelist(selIdc); 
    % optionally select specific videos
    selectSpecific = false;
    if selectSpecific
        % select more specific
        selIdc = contains({filelist.name},{'4854','4866','4880','4912'}) | contains({filelist.name},{'4853','4865','4879','4911'});
        selIdc = contains({filelist.name},{'4844','4858','4870','4904'}) | contains({filelist.name},{'4843','4857','4869','4903'});
        % apply selection
        filelist = filelist(selIdc);
    end
    if isempty(filelist)
        error('no files found.')
    end
    % gather files
    fprintf('loading data...\n')
    loadTimer = tic;
    clear cellPoints;
    for fileIdx = 1:length(filelist)
        fileFolder = filelist(fileIdx).folder;
        fileName = filelist(fileIdx).name;
        filePath = [fileFolder '\' fileName];

        % extract cell type
        slashIdc = strfind(fileFolder,'\');
        cellTypeStr = fileFolder(slashIdc(end)+1:end);
        if strcmp(cellTypeStr,'Healthy_RBCs')
            cellTypeNum = 1;
        end
        if strcmp(cellTypeStr,'Rigid_RBCs')
            cellTypeNum = 2;
        end

        % gather align vector
        alignVec = [0,0];
        try
            alignFiles = dir([filePath(1:end-length('peaks.mat')),'*alignVec.mat']);
            alignPath = [alignFiles(1).folder, '\', alignFiles(1).name];
            load(alignPath);
        catch
            fprintf('warning: unable to load alignment.\n')
        end
        if isnan(sum(alignVec))
            alignVec = [0 0];
        end

        % load peaks
        load(filePath);
        hrbcFrameRateFac = 1/2;
        rrbcFrameRateFac = 1;
        % skip frames from RRBC to equalize framerate
        skipFrames = true;
        if skipFrames
            if cellTypeNum == 2
                allpoints = allpoints([1:2:length(allpoints), 2:2:length(allpoints)]);
            end
            rrbcFrameRateFac = 1/2;
        end
        % save
        cellPoints(fileIdx).pnts = allpoints;
        cellPoints(fileIdx).cellType = cellTypeNum;
        cellPoints(fileIdx).alignVec = alignVec;
        cellPoints(fileIdx).filePath = filePath;
        cellPoints(fileIdx).fileFolder = fileFolder;
    end
    fprintf('loading done in %.2f seconds.\n',toc(loadTimer));
    fprintf('loaded files:\n');
    fprintf('--------------\n');
    for fileIdx = 1:length(cellPoints)
        fullPath = fullfile(cellPoints(fileIdx).fileFolder, cellPoints(fileIdx).filePath);
        fprintf('%d: %s\n', fileIdx, fullPath);
    end
    %% loop through bifurcations
    % take the folders of the peak files
    peakFolders = {cellPoints.fileFolder};
    cellTypeIdent = [cellPoints.cellType];
    % find data directories
    filelist = dir(fullfile(rootDir, '**\',maskName));  % get list of files and folders in any subfolder
    rootFolders = unique({filelist.folder});
    timer = tic;
    suffix = '_traj';
    for rootIdx = 1:length(rootFolders)
        rootdir = rootFolders{rootIdx};
        % load wall image
        wallImg = imread([rootdir,'\',maskName]);
        if length(size(wallImg)) > 2
            wallImg = wallImg(:,:,1);
        end
        wallBW = imbinarize(wallImg);
        wallBWfilt = bwmorph(wallBW,'thicken',3);
        % load the bifuractions
        load([rootdir,'\bifFlow.mat']); % HRBCs
        for cellTypeIdx = 1:2
            % find the corresponding peak files
            selIdc = strcmp([rootdir,'\',cellTypes{cellTypeIdx}],peakFolders);
            % struct to save everything
            clear clusters;
            kBifClu = ones(1,length(bifurcations));
            allAlign = [];
            for pntsIdx = find(selIdc)
                allAlign = [allAlign;cellPoints(pntsIdx).alignVec];
            end
            % take the peaks in that geometry
            for pntsIdx = find(selIdc)
                cPnts = cellPoints(pntsIdx);
                allPnts = cPnts.pnts;
                % adjust alignment
                alignVec = cPnts.alignVec;
                tic;
                pntNums = zeros(1,length(allPnts));
                for frameIdx = 1:length(allPnts)
                    framePnts = allPnts(frameIdx).peak;
                    pntNums(frameIdx) = size(framePnts,1);
                end
                borders = [0 cumsum(pntNums)];
                % make point cloud
                cloudPoints = zeros(sum(pntNums),3);
                inMask = zeros(sum(pntNums),1,'logical');
                for frameIdx = 1:length(allPnts)
                    framePnts = allPnts(frameIdx).peak;
                    if ~isempty(framePnts)
                        framePnts(:,1) = framePnts(:,1)+alignVec(1);
                        framePnts(:,2) = framePnts(:,2)+alignVec(2);
                        dist = 10;
                        try
                            inMask((borders(frameIdx)+1):borders(frameIdx+1),:) = diag(wallBWfilt(round(framePnts(:,1)),round(framePnts(:,2))));
                            cloudPoints((borders(frameIdx)+1):borders(frameIdx+1),:) = [framePnts, frameIdx * dist * ones(size(framePnts,1),1)];
                        end
                    end
                end
                toc;
                cloudPoints = cloudPoints(inMask,:); % filter by walls
                if flag_test
                    figure;
                    fprintf('plotting points in mask.\n')
                    scatter3(cloudPoints(:,1), cloudPoints(:,2), cloudPoints(:,3), 5, cloudPoints(:,3), 'filled');
                    title('Filtered Point Cloud');
                    xlabel('X');
                    ylabel('Y');
                    zlabel('Z');
                    axis equal;
                    grid on;
                    colorbar;
                    colormap('parula'); % Use a visually appealing colormap
                end
                %% go through all bifurcations
                for bifIdx = 1:length(bifurcations)
                    % find inlet vessel
                    cBif = bifurcations(bifIdx);
                    bifAngle = cBif.bifAvgAngle;
                    conVes = cBif.conVesBdy;
                    angleDif = zeros(1,length(conVes));
                    for vesIdx = 1:length(conVes)
                        vesStartPnts = flip(conVes(vesIdx).pnts(1:5,:),1);
                        dirVec = mean(diff(vesStartPnts));
                        vesAngle = -acos(dirVec(:,2)).*(1-2*(-dirVec(:,1)<0))/pi*180;
                        angleDif(vesIdx) = abs(bifAngle-vesAngle);
                        conVes(vesIdx).vesAngle = vesAngle;
                    end
                    % find inlet vessel
                    inVesIdx = find(angleDif == min(angleDif),1,'first');
                    selDaugh = ones(1,length(conVes),'logical');
                    selDaugh(inVesIdx) = 0;
                    bifurcations(bifIdx).inVesIdx = inVesIdx;
                    bifurcations(bifIdx).selDaugh = selDaugh;
                    originVesselIdx = inVesIdx;
                    try
                        vesPnts = conVes(originVesselIdx).pnts;
                    catch
                        vesPnts = [];
                    end
                    bifurcations(bifIdx).orVes = originVesselIdx;
                    if ~isempty(vesPnts)
                        originPnt = vesPnts(end,:);
                        bifurcationBW = zeros(size(wallBW),'logical');
                        flowMap = zeros(size(wallBW),'double');
                        angleMap = zeros(size(wallBW),'double');
                        for vesIdx = 1:length(conVes)
                            vesBW = conVes(vesIdx).vesBW;
                            %vesBW = bwmorph(vesBW,'thicken',3);
                            avgFlow = conVes(vesIdx).avgFlow;
                            avgAngle = conVes(vesIdx).avgAngle;
                            bifurcationBW = bifurcationBW | vesBW;
                            flowMap = flowMap + double(vesBW)*avgFlow;% * (1-2*inFlow(vesIdx));
                            angleMap = angleMap + double(vesBW)*avgAngle;
                        end
                        % add inner bifurcation data
                        bifBW = bifurcations(bifIdx).bifBW;
                        avgFlow = bifurcations(bifIdx).bifAvgFlow;
                        avgAngle = bifurcations(bifIdx).bifAvgAngle;
                        bifurcationBW = bifurcationBW | bifBW;
                        bifCtr = bifurcations(bifIdx).bifCtr;
                        % add flow
                        flowMap = flowMap + double(bifBW)*avgFlow;
                        angleMap = angleMap + double(bifBW)*avgAngle;
                        angleMap(~bifurcationBW) = nan;
                        % smooth flow map
                        smFlowMap = smoothdata(flowMap,1,'movmean',3);
                        smFlowMap = smoothdata(smFlowMap,2,'movmean',3);
                        % thicken bifurcation BW
                        %bifurcationBW = bwmorph(bifurcationBW,'thicken',8);
                        inMask = diag(bifurcationBW(round(cloudPoints(:,1)),round(cloudPoints(:,2))));
                        vesPoints = cloudPoints(inMask,:);
                        vesPointsTrafo = zeros(size(vesPoints));
                        for startPntIdx = 1:size(vesPoints,1)
                            pnt = vesPoints(startPntIdx,:);
                            distNorms = norm(originPnt-pnt(1:2));
                            pnt(3) = pnt(3)-distNorms/(3*smFlowMap(round(pnt(1)),round(pnt(2))));
                            vesPointsTrafo(startPntIdx,:) = pnt;
                        end

                        % align origin vessel with x-axis
                        vesIdx = originVesselIdx;
                        originAngle = ( (conVes(vesIdx).vesAngle)/180*pi)+pi/2+pi;%*(1-2*inFlow(vesIdx));
                        rotMat = [  cos(originAngle)    -sin(originAngle)   0;...
                            sin(originAngle)    cos(originAngle)    0;...
                            0                   0                   1];
                        globRotMat = rotMat;
                        distanceThres = 80; % maximum distance between points
                        for trajIter = 1%:2 % optionally add more iterations
                            if trajIter == 2
                                vesPoints = [];
                                for kClu = 1:length(clu)
                                    pnts = clu(kClu).pnts;
                                    vesPoints = [vesPoints; pnts];
                                    distanceThres = 30; % maximum distance between points
                                end
                            end
                            vesPointsO = vesPoints;
                            % align with x axis and bring to origin
                            vesPointsRot = (rotMat*vesPoints')';
                            minVesPointsRot = min(vesPointsRot);
                            vesPointsRot = vesPointsRot-minVesPointsRot;

                            allSelIdc = ones(1,size(vesPoints,1),'logical');
                            % main iteration
                            mainIterTimer = tic;
                            %% through trajectories loop
                            clear clu; kClu = 1;
                            while size(vesPointsRot,1) > 10
                                %% (1) trace single trajectory in one direction
                                fprintf('starting new trace number %d.\n',kClu);
                                vesPointsRot = vesPointsRot(allSelIdc,:);
                                allSelIdc = allSelIdc(allSelIdc);
                                conPnts = [];

                                iter = 1;
                                close all; fig = figure;
                                traceTraj = true;
                                while traceTraj
                                    if iter == 1
                                        %% (1) find start point
                                        startPntIdx = find(vesPointsRot(:,3)+vesPointsRot(:,1) == min(vesPointsRot(:,3)+vesPointsRot(:,1)),1,'first');
                                        startPntRot = vesPointsRot(startPntIdx,:);
                                        allSelIdc(startPntIdx) = 0;
                                    end
                                    if isempty(startPntRot)
                                        break;
                                    end
                                    %% (2) take points out
                                    vesPointsRot = vesPointsRot(allSelIdc,:);
                                    allSelIdc = allSelIdc(allSelIdc);
                                    %% (3) compute prediction of next point based on flow
                                    prediction_pnt = false;
                                    if prediction_pnt
                                        % transform starting point
                                        pntFlow = smFlowMap(round(startPnt(1)),round(startPnt(2)));
                                        pntAngle = angleMap(round(startPnt(1)),round(startPnt(2)));
                                        % direction vector
                                        rotMat = [  cos(pntAngle)    -sin(pntAngle)   0;...
                                            sin(pntAngle)    cos(pntAngle)    0;...
                                            0                   0                   1];
                                        dirVec = rotMat*[0; -1; 0];
                                    end
                                    % transform startpoint
                                    startPntTrafo = startPntRot+[1 0 1];
                                    % take into account previous direction
                                    prevPnts = conPnts(end+1-min([size(conPnts,1),3]):end,:);
                                    prevDir = mean(diff(prevPnts),1);
                                    if logical(~isnan(prevDir)) & logical(size(prevPnts,1) > 1)
                                        startPntTrafo = startPntRot + prevDir;
                                        prevDir = prevDir./norm(prevDir);
                                    else
                                        prevDir = [];
                                    end
                                    %% (4) find neighborhood to predicted point
                                    relVec = startPntTrafo - vesPointsRot;
                                    dists = sqrt(relVec(:,1).^2+relVec(:,2).^2+relVec(:,3).^2);
                                    nbPntLog = dists < distanceThres;
                                    nbPntLog = nbPntLog & (relVec(:,3) < 2*dist); % take out negative time
                                    if sum(nbPntLog) == 0
                                        break;
                                    end
                                    dists = dists(nbPntLog);
                                    nbPnts = vesPointsRot(nbPntLog,:);
                                    % relative vectors
                                    nbVec = (nbPnts-startPntRot);
                                    nbDists = sqrt(nbVec(:,1).^2+nbVec(:,2).^2+nbVec(:,3).^2);
                                    nbDir = nbVec./nbDists;
                                    % weighting by angle
                                    if ~isempty(prevDir)
                                        cosTheta = zeros(1,size(nbDir,1));
                                        for nbIdx = 1:size(nbDir,1)
                                            cosTheta(nbIdx) = dot(prevDir,nbDir(nbIdx,:));
                                        end
                                        nbAngleScore = (cosTheta/2+0.5).^(2);
                                    else
                                        nbAngleScore = ones(1,size(nbDir,1));
                                    end
                                    nbDistScore = 1./nbDists';
                                    %% (5) compute distance and angle weighted average
                                    avgVec = sum(nbVec.*repmat(nbDistScore',1,3).*repmat(nbAngleScore',1,3),1)./sum(nbDistScore.*nbAngleScore);
                                    totalScore = nbDistScore.*nbAngleScore;
                                    if (max(totalScore)<0.02)
                                        break
                                    end
                                    maxIdx = find(totalScore==max(totalScore),1,'first');
                                    maxVec = nbVec(maxIdx,:);
                                    %% (6) set new starting point
                                    conPnts = [conPnts;startPntRot];
                                    startPntRot = startPntRot + maxVec;
                                    %% (7) take out new point
                                    nbPntIdc = find(nbPntLog);
                                    allSelIdc(nbPntIdc(maxIdx)) = 0;
                                    %% plotting
                                    if ~isinf(framePause)
                                        plot3(vesPointsRot(:,1),vesPointsRot(:,2),vesPointsRot(:,3),'.b');hold on;
                                        plot3(startPntRot(:,1),startPntRot(:,2),startPntRot(:,3),'xr');
                                        plot3(nbPnts(:,1),nbPnts(:,2),nbPnts(:,3),'og');
                                        plot3(conPnts(:,1),conPnts(:,2),conPnts(:,3),'o-k');
                                        hold off;
                                        try
                                            xlim(startPntRot(:,1)+[-200,200]);
                                            ylim(startPntRot(:,2)+[-200,200]);
                                            zlim(startPntRot(:,3)+[-200,200]);
                                        end
                                        fig.Position(3:4) = 1000;
                                        fig.Position(2) = 10;
                                        v = [-5 -2 5];
                                        [caz,cel] = view(v);
                                        fprintf('totScore %f\n',max(totalScore));
                                        pause(framePause)
                                    else
                                        %close all;
                                    end

                                    %% iteration
                                    iter = iter+1;
                                end
                                oConPnts = conPnts;
                                fprintf('trajectory traced with length %d.\n',size(conPnts,1));
                                %% (2) trace in other direction (false linking remedy)
                                fprintf('starting backtrace.\n');
                                if ~isempty(oConPnts)
                                    % add trajectory again
                                    vesPointsRot = [vesPointsRot;conPnts];
                                    allSelIdc = [allSelIdc,ones(1,size(conPnts,1),'logical')];
                                    lastPrevDir = prevDir;
                                    conPnts = [];

                                    iter = 1;
                                    close all; fig = figure;
                                    traceTraj = true;
                                    while traceTraj
                                        if iter == 1
                                            %% (1) find start point
                                            startPntIdx = size(vesPointsRot,1);
                                            startPntRot = vesPointsRot(startPntIdx,:);
                                            allSelIdc(startPntIdx) = 0;
                                        end
                                        if isempty(startPntRot)
                                            break;
                                        end
                                        %% (2) take points out
                                        vesPointsRot = vesPointsRot(allSelIdc,:);
                                        allSelIdc = allSelIdc(allSelIdc);
                                        %% (3) compute prediction of next point based on flow
                                        prediction_pnt = false;
                                        if prediction_pnt
                                            % transform starting point
                                            pntFlow = smFlowMap(round(startPnt(1)),round(startPnt(2)));
                                            pntAngle = angleMap(round(startPnt(1)),round(startPnt(2)));
                                            % direction vector
                                            rotMat = [  cos(pntAngle)    -sin(pntAngle)   0;...
                                                sin(pntAngle)    cos(pntAngle)    0;...
                                                0                   0                   1];
                                            dirVec = rotMat*[0; -1; 0];
                                        end
                                        % transform start point
                                        startPntTrafo = startPntRot-[1 0 1];
                                        % take into account previous direction
                                        prevPnts = conPnts(end+1-min([size(conPnts,1),3]):end,:);
                                        if iter == 1
                                            prevDir = -lastPrevDir;
                                        else
                                            prevDir = mean(diff(prevPnts),1);
                                        end
                                        if logical(~isnan(prevDir)) & logical(size(prevPnts,1) > 1)
                                            startPntTrafo = startPntRot + prevDir;
                                            prevDir = prevDir./norm(prevDir);
                                        else
                                            prevDir = [];
                                        end
                                        %% (4) find neighborhood to predicted point
                                        relVec = startPntTrafo - vesPointsRot;
                                        dists = sqrt(relVec(:,1).^2+relVec(:,2).^2+relVec(:,3).^2);
                                        nbPntLog = dists < distanceThres;
                                        nbPntLog = nbPntLog & (relVec(:,3) > -2*dist); % take out negative time
                                        if sum(nbPntLog) == 0
                                            break;
                                        end
                                        dists = dists(nbPntLog);
                                        nbPnts = vesPointsRot(nbPntLog,:);
                                        % relative vectors
                                        nbVec = (nbPnts-startPntRot);
                                        nbDists = sqrt(nbVec(:,1).^2+nbVec(:,2).^2+nbVec(:,3).^2);
                                        nbDir = nbVec./nbDists;
                                        % weighting by angle
                                        if ~isempty(prevDir)
                                            cosTheta = zeros(1,size(nbDir,1));
                                            for nbIdx = 1:size(nbDir,1)
                                                cosTheta(nbIdx) = dot(prevDir,nbDir(nbIdx,:));
                                            end
                                            nbAngleScore = (cosTheta/2+0.5).^(2);
                                        else
                                            nbAngleScore = ones(1,size(nbDir,1));
                                        end
                                        nbDistScore = 1./nbDists';
                                        %% (5) compute distance and angle weighted average
                                        avgVec = sum(nbVec.*repmat(nbDistScore',1,3).*repmat(nbAngleScore',1,3),1)./sum(nbDistScore.*nbAngleScore);
                                        totalScore = nbDistScore.*nbAngleScore;
                                        if (max(totalScore)<0.02)
                                            break
                                        end
                                        maxIdx = find(totalScore==max(totalScore),1,'first');
                                        maxVec = nbVec(maxIdx,:);
                                        %% (6) set new starting point
                                        conPnts = [conPnts;startPntRot];
                                        startPntRot = startPntRot + maxVec;
                                        %% (7) take out new point
                                        nbPntIdc = find(nbPntLog);
                                        allSelIdc(nbPntIdc(maxIdx)) = 0;
                                        %% plotting
                                        if ~isinf(framePause)
                                            plot3(vesPointsRot(:,1),vesPointsRot(:,2),vesPointsRot(:,3),'.b');hold on;
                                            plot3(startPntRot(:,1),startPntRot(:,2),startPntRot(:,3),'xr');
                                            plot3(nbPnts(:,1),nbPnts(:,2),nbPnts(:,3),'og');
                                            plot3(oConPnts(:,1),oConPnts(:,2),oConPnts(:,3),'color',[0.9 0.9 0.9]);
                                            plot3(oConPnts(1,1),oConPnts(1,2),oConPnts(1,3),'s','color',[0.9 0.0 0.0]);
                                            plot3(conPnts(:,1),conPnts(:,2),conPnts(:,3),'o-k');
                                            hold off;
                                            try
                                                xlim(startPntRot(:,1)+[-200,200]);
                                                ylim(startPntRot(:,2)+[-200,200]);
                                                zlim(startPntRot(:,3)+[-200,200]);
                                            end
                                            fig.Position(3:4) = 1000;
                                            fig.Position(2) = 10;
                                            v = [-5 -2 5];
                                            [caz,cel] = view(v);
                                            fprintf('totScore %f\n',max(totalScore));
                                            pause(framePause)
                                        else
                                            %close all;
                                        end
                                        %% iteration
                                        iter = iter+1;
                                    end
                                end
                                fprintf('backtrace finished with trajectory of length %d.\n',size(conPnts,1));
                                % save trajectory
                                try
                                    cluPnts = (globRotMat'*(conPnts(2:end,:)+minVesPointsRot)')';
                                    [~,sortIdc] = sort(cluPnts(:,3));
                                    cluPnts = cluPnts(sortIdc,:);
                                    clu(kClu).pnts = cluPnts;
                                catch
                                    clu(kClu).pnts = [];
                                end

                                kClu = kClu + 1;

                                % status display
                                if ~mod(size(vesPointsRot,1)-1,round(length(allSelIdc)/10))
                                    fprintf('status: %.2f percent.\n',abs(1-size(vesPoints,1)/length(allSelIdc))*100);
                                    fprintf('remaining: %.2f seconds.\n',toc(mainIterTimer)/sum(allSelIdc) * sum(~allSelIdc));
                                    fprintf('%d clusters found.\n',length(clu));
                                end
                            end
                            if ~exist('clu','var')
                                clu = [];
                            end
                        end
                        fprintf('main iteration finished in %.2f seconds.\n',toc(mainIterTimer));
                        fprintf('%d trajectories found.\n',length(clu));
                        clc;
                        fprintf('-------------------------total time %.2f ---------------------\n',toc(timer));

                        %% link trajectories
                        linkingTraj = false;
                        if linkingTraj
                            % first moving average filter
                            clear filtclu; kClu = 1;
                            for cluIdx = 1:length(clu)
                                pnts = clu(cluIdx).pnts;
                                if ~isempty(pnts) && size(pnts,1) > 3
                                    pnts = movmean(pnts,3);
                                    filtclu(kClu).pnts = pnts;
                                    kClu = kClu + 1;
                                end
                            end
                            if exist("filtclu",'var')
                                clu = filtclu;
                            else
                                clu = [];
                            end
                            for cluIdx = 1:length(clu)
                                pnts = clu(cluIdx).pnts;
                                % extrapolate start
                                dirVec = mean(diff(pnts(1:3,:)));
                                clu(cluIdx).startExtr = pnts(1,:)-dirVec;
                                % extrapolate end
                                dirVec = mean(diff(pnts(end+1-(1:3),:)));
                                clu(cluIdx).endExtr = pnts(end,:)+dirVec;
                            end
                            allEnd = zeros(size(clu,1),3);
                            allStart = zeros(size(clu,1),3);
                            for cluIdx = 1:length(clu)
                                allEnd(cluIdx,:) = clu(cluIdx).endExtr;
                                allStart(cluIdx,:) = clu(cluIdx).startExtr;
                            end
                            selectIdt = ones(1,size(allStart,1),'logical');
                            remainIdt = ones(1,size(allStart,1),'logical');
                            for endIdx = 1:size(allEnd,1)
                                selectIdt(remainIdt==0) = 0;
                                selIdc = find(selectIdt);
                                relVec = allEnd(endIdx,:)-allStart(selectIdt,:);
                                dists = sqrt(relVec(:,1).^2+relVec(:,2).^2+relVec(:,3).^2);
                                minIdx = find(dists==min(dists),1,'first');
                                minDist = dists(minIdx);
                                if minDist < 20
                                    % link
                                    clu(selIdc(minIdx)).pnts = [clu(endIdx).pnts;clu(selIdc(minIdx)).pnts];
                                    selectIdt(selIdc(minIdx)) = 0;
                                    remainIdt(endIdx) = 0;
                                end
                            end
                            clu = clu(remainIdt);
                        end
                        %% measure features
                        for cluIdx = 1:length(clu)
                            pnts = clu(cluIdx).pnts;
                            % (1) how many neighbors in original cloud
                            neigh = zeros(1,size(pnts,1));
                            for startPntIdx = 1:size(pnts,1)
                                relVec = pnts(startPntIdx,:)-vesPointsO;
                                dists = sqrt(relVec(:,1).^2+relVec(:,2).^2+relVec(:,3).^2);
                                neigh(startPntIdx) = length(find(dists < distanceThres));
                            end
                            clu(cluIdx).neigh = neigh;
                        end
                        %% plot everything
                        plot_clusters = false;
                        if flag_test
                            plot_clusters = true;
                        end
                        if plot_clusters
                            close all
                            imshow(wallBW')
                            hold on
                            plot3(vesPointsO(:,1),vesPointsO(:,2),vesPointsO(:,3),'.','MarkerSize',5,'color',[0.7 0.7 0.7]);
                            allPntsZ = [];
                            for kClu = 1:length(clu)
                                pnts = clu(kClu).pnts;

                                if ~isempty(pnts) %&& size(pnts,1) > 20
                                    allPntsZ = [allPntsZ;pnts(1,3),pnts(end,3)];
                                    plot3(pnts(:,1),pnts(:,2),pnts(:,3),'-','LineWidth',4)
                                end
                            end
                            hold off
                            return
                        end
                        %% spline fit filter
                        spline_fit_filter = false;
                        if spline_fit_filter
                            clear filtclu; k = 1;
                            FLAGdispspline = 0;
                            for cluIdx = 1:length(clu)
                                cpoints = clu(cluIdx).pnts;
                                if ~isempty(cpoints) && size(cpoints,1) > 1
                                    x = cpoints(:,1);
                                    y = cpoints(:,2);
                                    z = cpoints(:,3);
                                    % diff
                                    dpoints = diff(cpoints);
                                    trajLength = sum(sqrt(dpoints(:,1).^2+dpoints(:,2).^2+dpoints(:,3).^2));
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
                                    knotnumber = max([5,round(trajLength/30)]);
                                    knotarray = linspace(min(didc),max(didc),knotnumber);
                                    spoints = zeros(Nspline,3);
                                    for idx = 1:3
                                        slm = slmengine(didc,fpoints(:,idx),'plot','off','knots',knotarray);
                                        spoints(:,idx) = slmeval(sidc,slm);
                                    end
                                    dpoints = diff(spoints);
                                    trajLength = sum(sqrt(dpoints(:,1).^2+dpoints(:,2).^2));
                                    filtclu(k).pnts = spoints;
                                    filtclu(k).length = trajLength;
                                    filtclu(k).avgV = abs(trajLength./(spoints(end,3)-spoints(1,3)));
                                    k = k+1;
                                end
                            end
                            plot_clusters = false;
                            if plot_clusters
                                close all
                                figure
                                hold on
                                plot3(vesPointsO(:,1),vesPointsO(:,2),vesPointsO(:,3),'.b','MarkerSize',5);
                                allPntsZ = [];
                                for kClu = 1:length(filtclu)
                                    pnts = filtclu(kClu).pnts;
                                    v = filtclu(kClu).avgV;
                                    if ~isempty(pnts)
                                        allPntsZ = [allPntsZ,v];
                                        plot3(pnts(:,1),pnts(:,2),pnts(:,3),'-','LineWidth',4)
                                    end
                                end
                                hold off
                                figure
                                ksdensity(v)
                            end
                        end
                        % take out distance again
                        for cluIdx = 1:length(clu)
                            pnts = clu(cluIdx).pnts;
                            if ~isempty(pnts)
                                pnts(:,3) = pnts(:,3)/dist;
                            end
                            clu(cluIdx).pnts = pnts;
                        end
                        %% saving
                        saveTimer = tic;
                        for idxClu = 1:length(clu)
                            if ~isempty(clu(idxClu).pnts)
                                clusters(bifIdx).clu(kBifClu(bifIdx)) = clu(idxClu);
                                kBifClu(bifIdx) = kBifClu(bifIdx) + 1;
                            end
                        end
                    else
                        saveTimer = tic;
                    end
                    fprintf('saving done in %.2f seconds.\n',toc(saveTimer));
                end
            end

            if cellTypeIdx == 1
                fprintf('merging all healthy RBC data.\n')
                % merge all clusters of this geometry
                for bifIdx = 1:length(bifurcations)
                    cClu = [];
                    if kBifClu(bifIdx) > 1
                        cClu = clusters(bifIdx).clu;
                        bifurcations(bifIdx).HRBCtraj = cClu;
                    else
                        bifurcations(bifIdx).HRBCtraj = [];
                    end
                    fprintf('%d clusters in bifurcation %d.\n',length(cClu),bifIdx);
                end
            else
                fprintf('merging all rigid RBC data.\n')
                % merge all clusters of this geometry
                for bifIdx = 1:length(bifurcations)
                    cClu = [];
                    if kBifClu(bifIdx) > 1
                        cClu = clusters(bifIdx).clu;
                        bifurcations(bifIdx).RRBCtraj = cClu;
                    else
                        bifurcations(bifIdx).RRBCtraj = [];
                    end
                    fprintf('%d clusters in bifurcation %d.\n',length(cClu),bifIdx);
                end
            end
        end
        %% determine bifurcation average velocity magnitude
        for bifIdx = 1:length(bifurcations)
            bifurcations(bifIdx).modeThres = 0;
            RRBCtraj = bifurcations(bifIdx).RRBCtraj;
            allRRBCvel = [];
            for trajIdx = 1:length(RRBCtraj)
                traj = RRBCtraj(trajIdx);
                if ~isempty(traj.pnts)
                    neigh = traj.neigh;
                    pnts =  traj.pnts;
                    velVec = diff(pnts(:,1:2));
                    if length(velVec(:)) > 3
                        velVec = velVec./repmat(diff(pnts(:,3)),1,2);
                        velMag = sqrt(velVec(:,1).^2+velVec(:,2).^2);
                        velMag = velMag(~isnan(velMag) & ~isinf(velMag));
                        avgMag = mean(velMag);
                        allRRBCvel = [allRRBCvel;avgMag];
                    end
                end
            end
            allRRBCvel = allRRBCvel(allRRBCvel > 0);

            % healthy
            HRBCtraj = bifurcations(bifIdx).HRBCtraj;
            allHRBCvel = [];
            for trajIdx = 1:length(HRBCtraj)
                traj = HRBCtraj(trajIdx);
                if ~isempty(traj.pnts)
                    neigh = traj.neigh;
                    pnts =  traj.pnts;
                    velVec = diff(pnts(:,1:2));
                    if length(velVec(:)) > 3
                        velVec = velVec./repmat(diff(pnts(:,3)),1,2);
                        velMag = sqrt(velVec(:,1).^2+velVec(:,2).^2);
                        velMag = velMag(~isnan(velMag) & ~isinf(velMag));
                        avgMag = mean(velMag);
                        allHRBCvel = [allHRBCvel;avgMag];
                    end
                end
            end
            allHRBCvel = allHRBCvel(allHRBCvel > 0);

            bifurcations(bifIdx).avgHRBCvel = median(allHRBCvel);
            bifurcations(bifIdx).avgRRBCvel = median(allRRBCvel);

        end
        %% saving
        fprintf('saving all data.\n')
        save([rootdir,'\bifTraj.mat'],"bifurcations");
    end
    toc(timer);
    % save workspace
    save([rootdir,'\TrackingWorkspace.mat']);
end
load([rootdir,'\TrackingWorkspace.mat']); % with skipping
%% plotting and further analysis of velocities in different branches
close all;
frameRate =  99.9255;
dist = 1; % always set to 1
deltaT = 1/frameRate;
microScale = 0.65; % micron/px: 20x is 0.65 µm/px, and 50x is 0.26 µm/px
velScale = microScale * dist/deltaT;
velScaleHRBC = velScale*hrbcFrameRateFac;
velScaleRRBC = velScale*rrbcFrameRateFac;
velScaleFlow = velScale*hrbcFrameRateFac;
wallImg = imread([rootdir,'\',maskName]);
fig = figure;
imshow(wallImg)
fig.Units = 'pixels';
fig.Position = [100 100 size(wallImg,2),size(wallImg,1)];
ax = gca;
ax.Position([1,2]) = 0;
ax.Position([3,4]) = 1;
drawArrow = @(x,y) quiver3( x(1),y(1),10000,x(2)-x(1),y(2)-y(1),0,0 ,'Color',[0 0 0],'LineWidth',1,'AutoScale','off');
hold on
vesCtrPnts = [];
allHRBCvel = [];
allRRBCvel = [];
allFlow = [];

for bifIdx = 1:length(bifurcations)
    cPnt = bifurcations(bifIdx).bifCtr;
    avgFlow = bifurcations(bifIdx).bifAvgFlow;
    bifAngle = bifurcations(bifIdx).bifAvgAngle;
    % add bifurcation
    cAngle = bifAngle;
    dirVec = [sin(cAngle/180*pi),cos(cAngle/180*pi)];
    dPnt = cPnt + 0;
    pos = get(gca, 'Position');
    xArrow = [cPnt(2),dPnt(2)];
    yArrow = size(wallImg,1)-[cPnt(1),dPnt(1)];
    HRBCVel = bifurcations(bifIdx).avgHRBCvel*velScaleHRBC;
    RRBCVel = bifurcations(bifIdx).avgRRBCvel*velScaleRRBC;
    HRBCTraj = bifurcations(bifIdx).HRBCtraj;
    RRBCTraj = bifurcations(bifIdx).RRBCtraj;
    % add vessels
    conVes = bifurcations(bifIdx).conVesBdy;
    bifurcations(bifIdx).velScaleHRBC = velScaleHRBC;
    bifurcations(bifIdx).velScaleRRBC = velScaleRRBC;
    
    inVesIdx = bifurcations(bifIdx).orVes;

    hold on
    for vesIdx = 1:length(conVes)
        % find boundaries
        vesBW = conVes(vesIdx).vesBW;
        vesAx = conVes(vesIdx).pnts;
        vesStart = vesAx(1,:);
        %realVesLength = max(conVes(vesIdx).length);
        medPnt = median(vesAx);
        vesLength = 2*norm(medPnt-vesStart);
        vesBounds = bwboundaries(vesBW);
        vesBounds = vesBounds{1};
        % average flow
        avgFlow = conVes(vesIdx).avgFlow;
        %% entire length of vessel velocities
        % HRBC velocity
        HRBCVel = [];
        HRBCDist = [];
        HRBCTime = [];
        clear HRBCVesTraj; kTraj = 1;
        for trajIdx = 1:length(HRBCTraj)
            pnts = HRBCTraj(trajIdx).pnts;
            % first filter unsensical trajectories
            trueTrajCond = [0;diff(pnts(:,3))]>=0;
            trueTrajCond = smooth([flip(trueTrajCond);trueTrajCond;flip(trueTrajCond)],3);
            trueTrajCond = trueTrajCond(size(pnts,1)+(1:size(pnts,1)))>0.7;trueTrajCond  = (1:size(pnts,1))>0;
            pnts = pnts(trueTrajCond,:);
            % check which points are inside
            inboundIdc = inpolygon(pnts(:,1),pnts(:,2),vesBounds(:,1),vesBounds(:,2));
            pnts = pnts(inboundIdc,:);
            distances = sqrt((pnts(:,1)-vesStart(1)).^2+(pnts(:,2)-vesStart(2)).^2) - norm(vesStart-medPnt);% distances = distances - min(distances);
            %proxIdc = 1:size(pnts,1);%distances > 1/2*norm(vesStart-medPnt) & distances < 3/2*norm(vesStart-medPnt); % to get the entire length, take all indices
            proxIdc = 1:size(pnts,1);
            proxIdc = distances > -50 & distances < 0;
            
            pnts = real(pnts(proxIdc,:));
            if size(pnts,1) > 2
                HRBCVesTraj(kTraj).pnts = pnts;kTraj = kTraj + 1;
                ds = norm(pnts(end,1:2)-pnts(1,1:2));
                dt = pnts(end,3)-pnts(1,3);
                diffPnts = diff(pnts);
                dSpace = sqrt(diffPnts(:,1).^2 + diffPnts(:,2).^2);
                dTime = abs(diff(pnts(:,3)));
                dSpaceMaxIdx = find(dSpace==max(dSpace),1,'first');

                trajPntDistances = zeros(size(pnts,1),size(pnts,1));
                for id1 = 1:size(pnts,1)
                    for id2 = 1:size(pnts,1)
                        trajPntDistances(id1,id2) = norm(pnts(id1,1:2)-pnts(id2,1:2));
                    end
                end
                [max1,max2] = find(trajPntDistances==max(trajPntDistances(:)));
                max1 = max1(1);
                max2 = max2(1);
                maxTrajDistance = max(trajPntDistances(:));
                maxTrajDistanceTime = abs(pnts(max2,3)-pnts(max1,3));
                ds = maxTrajDistance;
                dt = maxTrajDistanceTime;

                HRBCDist = [HRBCDist,ds];
                HRBCTime = [HRBCTime,dt];

                velMag = maxTrajDistance/maxTrajDistanceTime;
                avgMag = mean(velMag);
                if imag(complex(avgMag))>0
                    return
                end
                HRBCVel = [HRBCVel,avgMag];
            end
        end
        if ~exist('HRBCVesTraj','var')
            HRBCVesTraj = [];
        end
        cond = ~isnan(HRBCVel)&~isinf(HRBCVel);
        HRBCVel = HRBCVel(cond);
        HRBCNum = length(HRBCVel);
        HRBCVelArray = HRBCVel;
        HRBCVelErr = std(HRBCVel)/sqrt(length(HRBCVel))*velScaleHRBC;
        % save
        bifurcations(bifIdx).conVesBdy(vesIdx).HRBCVel = HRBCVel;
        bifurcations(bifIdx).conVesBdy(vesIdx).HRBCNum = HRBCNum;
        bifurcations(bifIdx).conVesBdy(vesIdx).HRBCVelArray = HRBCVelArray;
        bifurcations(bifIdx).conVesBdy(vesIdx).HRBCVelErr = HRBCVelErr;
        bifurcations(bifIdx).conVesBdy(vesIdx).HRBCVesTraj = HRBCVesTraj;

        % RRBC velocity
        RRBCVel = [];
        RRBCDist = [];
        RRBCTime = [];
        clear RRBCVesTraj; kTraj = 1;
        for trajIdx = 1:length(RRBCTraj)
            pnts = RRBCTraj(trajIdx).pnts;
            % first filter unsensical trajectories
            trueTrajCond = [0;diff(pnts(:,3))]>=0;
            trueTrajCond = smooth([flip(trueTrajCond);trueTrajCond;flip(trueTrajCond)],3);
            trueTrajCond = trueTrajCond(size(pnts,1)+(1:size(pnts,1)))>0.7;trueTrajCond  = (1:size(pnts,1))>0;
            pnts = pnts(trueTrajCond,:);
            % check which points are inside
            inboundIdc = inpolygon(pnts(:,1),pnts(:,2),vesBounds(:,1),vesBounds(:,2));
            pnts = pnts(inboundIdc,:);
            distances = sqrt((pnts(:,1)-vesStart(1)).^2+(pnts(:,2)-vesStart(2)).^2) - norm(vesStart-medPnt);% distances = distances - min(distances);
            proxIdc = distances > -50 & distances < 0;
            pnts = real(pnts(proxIdc,:));
            if size(pnts,1) > 2
                RRBCVesTraj(kTraj).pnts = pnts;kTraj = kTraj + 1;
                ds = norm(pnts(end,1:2)-pnts(1,1:2));
                dt = pnts(end,3)-pnts(1,3);
                diffPnts = diff(pnts);
                dSpace = sqrt(diffPnts(:,1).^2 + diffPnts(:,2).^2);
                dTime = abs(diff(pnts(:,3)));
                dSpaceMaxIdx = find(dSpace==max(dSpace),1,'first');
                % from bounding box
                trajPntDistances = zeros(size(pnts,1),size(pnts,1));
                for id1 = 1:size(pnts,1)
                    for id2 = 1:size(pnts,1)
                        trajPntDistances(id1,id2) = norm(pnts(id1,1:2)-pnts(id2,1:2));
                    end
                end
                [max1,max2] = find(trajPntDistances==max(trajPntDistances(:)));
                max1 = max1(1);
                max2 = max2(1);
                maxTrajDistance = max(trajPntDistances(:));
                maxTrajDistanceTime = abs(pnts(max2,3)-pnts(max1,3));
                ds = maxTrajDistance;
                dt = maxTrajDistanceTime;

                RRBCDist = [RRBCDist,ds];
                RRBCTime = [RRBCTime,dt];

                velMag = maxTrajDistance/maxTrajDistanceTime;
                avgMag = mean(velMag);
                RRBCVel = [RRBCVel,avgMag];
            end
        end
        if ~exist('RRBCVesTraj','var')
            RRBCVesTraj = [];
        end
        cond = ~isnan(RRBCVel)&~isinf(RRBCVel);
        RRBCVel = RRBCVel(cond);
        RRBCNum = length(RRBCVel);
        RRBCVelArray = RRBCVel;
        RRBCVelErr = std(RRBCVel)/sqrt(length(RRBCVel))*velScaleRRBC;
        % save
        bifurcations(bifIdx).conVesBdy(vesIdx).RRBCVel = RRBCVel;
        bifurcations(bifIdx).conVesBdy(vesIdx).RRBCNum = RRBCNum;
        bifurcations(bifIdx).conVesBdy(vesIdx).RRBCVelArray = RRBCVelArray;
        bifurcations(bifIdx).conVesBdy(vesIdx).RRBCVelErr = RRBCVelErr;
        % lengths and times
        bifurcations(bifIdx).conVesBdy(vesIdx).HRBCDist = HRBCDist;
        bifurcations(bifIdx).conVesBdy(vesIdx).HRBCTime = HRBCTime;
        bifurcations(bifIdx).conVesBdy(vesIdx).RRBCDist = RRBCDist;
        bifurcations(bifIdx).conVesBdy(vesIdx).RRBCTime = RRBCTime;
        bifurcations(bifIdx).conVesBdy(vesIdx).RRBCVesTraj = RRBCVesTraj;
        %% kernel density estimation
        try
        fprintf('starting kernel density analysis for bif %d ves %d.\n',bifIdx,vesIdx);
        % optionally estimate h value
        %VelVals = HRBCTime;
        %[~,~,bandwidth] = ksdensity(VelVals(~isnan(VelVals)&~isinf(VelVals)&VelVals>0));%,'support','positive','boundarycorrection','reflection');
        %h = max([1,round(bandwidth)]);
        h = 6;
        N = 101;
        maxTime = max([HRBCTime,RRBCTime]) + 5;
        minTime = 0;
        maxDist = max([HRBCDist,RRBCDist]) + 20;
        minDist = 0;
        x = linspace(minTime,maxTime,N);
        y = linspace(minDist,maxDist,N);
        xScale = median(diff(x));
        yScale = median(diff(y));

        % HRBC
        dataPnts = [HRBCTime'-minTime,HRBCDist'-minDist];
        dataPnts = dataPnts(HRBCDist>max([prctile(HRBCDist,50),prctile(RRBCDist,50)])-max(std(HRBCDist),std(RRBCDist)),:);
        

        [X,Y] = meshgrid(x,y);
        kde2d = zeros(size(X));
        kernel = @(r,h) exp(-r.^2/h.^2);
        for pntIdx = 1:size(dataPnts,1)
            xPeak = dataPnts(pntIdx,1);
            yPeak = dataPnts(pntIdx,2);
            for dx = -2*h:2*h
                for dy = -2*h:2*h
                    r = sqrt((dx)^2+(dy)^2);
                    try
                        kde2d(round(xPeak/xScale+dx),round((yPeak/yScale+dy))) = kde2d(round(xPeak/xScale+dx),round((yPeak/yScale+dy)))+kernel(r,h);
                    end
                end
            end
        end

        HRBCkde = kde2d';

        % RRBC
        dataPnts = [RRBCTime'-minTime,RRBCDist'-minDist];
        dataPnts = dataPnts(RRBCDist>max([prctile(HRBCDist,50),prctile(RRBCDist,50)])-max(std(HRBCDist),std(RRBCDist)),:);

        [X,Y] = meshgrid(x,y);
        kde2d = zeros(size(X));
        kernel = @(r,h) exp(-r.^2/h.^2);
        for pntIdx = 1:size(dataPnts,1)
            xPeak = dataPnts(pntIdx,1);
            yPeak = dataPnts(pntIdx,2);
            for dx = -2*h:2*h
                for dy = -2*h:2*h
                    r = sqrt((dx)^2+(dy)^2);
                    try
                        kde2d(round(xPeak/xScale+dx),round((yPeak/yScale+dy))) = kde2d(round(xPeak/xScale+dx),round((yPeak/yScale+dy)))+kernel(r,h);
                    end
                end
            end
        end

        RRBCkde = kde2d';
        %% distribution edges
        edgeThres = mean(0.5*mean(HRBCkde(HRBCkde>0)));
        bwKdeImg = HRBCkde>edgeThres;
        bwKdeImg(1,:) = 0;        bwKdeImg(end,:) = 0;        bwKdeImg(:,1) = 0;        bwKdeImg(:,end) = 0;
        HRBCedge = edge(bwKdeImg,'roberts');
        [yHRBCedge,xHRBCedge] = find(HRBCedge);
        xHRBCedge = x(xHRBCedge);
        yHRBCedge = y(yHRBCedge);

        edgeThres = mean(0.5*mean(RRBCkde(RRBCkde>0)));
        bwKdeImg = RRBCkde>edgeThres;
        bwKdeImg(1,:) = 0;        bwKdeImg(end,:) = 0;        bwKdeImg(:,1) = 0;        bwKdeImg(:,end) = 0;
        RRBCedge = edge(bwKdeImg,'roberts');
        [yRRBCedge,xRRBCedge] = find(RRBCedge);
        xRRBCedge = x(xRRBCedge);
        yRRBCedge = y(yRRBCedge);

        computeByPerc = false;
        maximaKde = false;
        if computeByPerc
            %% density percentile as function of length
            % maximum positions
            HRBCKdeMaxPos = zeros(1,N);
            RRBCKdeMaxPos = zeros(1,N);
            for kdeIdx = 1:N
                HRBCKdeMaxPos(kdeIdx) = find(HRBCkde(kdeIdx,:) == max(HRBCkde(kdeIdx,:)),1,'first');
                RRBCKdeMaxPos(kdeIdx) = find(RRBCkde(kdeIdx,:) == max(RRBCkde(kdeIdx,:)),1,'first');
            end
            % origin speed
            [yOrigin,xOrigin] = find(RRBCkde == max(RRBCkde(:)),1,'first');
            % corresponding percentiles
            HRBCpercentile = sum(HRBCkde(yOrigin,1:xOrigin))/sum(HRBCkde(yOrigin,:))*100; % in percent
            RRBCpercentile = sum(RRBCkde(yOrigin,1:xOrigin))/sum(RRBCkde(yOrigin,:))*100; % in percent
            % all percentile positions
            HRBCKdePrcPos = zeros(1,N);
            RRBCKdePrcPos = zeros(1,N);
            for kdeIdx = 1:N
                % HRBC
                HRBCcdf = cumsum(HRBCkde(kdeIdx,:))/sum(HRBCkde(kdeIdx,:))*100; % percent
                prctilePos = find(HRBCcdf > HRBCpercentile,1,'first');
                if isempty(prctilePos)
                    prctilePos = 0;
                end
                HRBCKdePrcPos(kdeIdx) = prctilePos;
                % RRBC
                RRBCcdf = cumsum(RRBCkde(kdeIdx,:))/sum(RRBCkde(kdeIdx,:))*100; % percent
                prctilePos = find(RRBCcdf > RRBCpercentile,1,'first');
                if isempty(prctilePos)
                    prctilePos = 0;
                end
                RRBCKdePrcPos(kdeIdx) = prctilePos;
            end
            HRBCKdePrcPos(HRBCKdePrcPos==0) = nan;
            RRBCKdePrcPos(HRBCKdePrcPos==0) = nan;

            flag_plotPercentiles = false;
            if flag_plotPercentiles
                close all;
                figure
                hold on
                plot(HRBCKdeMaxPos)
                plot(HRBCKdePrcPos)
                plot(RRBCKdeMaxPos)
                plot(RRBCKdePrcPos)
                xline(yOrigin)
                hold off

                legend('HRBCMax','HRBCPrc','RRBCMax','RRBCPrc');
                return
            end

            %% percentile position fit

            smParam = 0.001;

            xIdc = HRBCKdePrcPos;
            yIdc = 1:length(HRBCKdePrcPos);
            cond = ~isnan(xIdc)&xIdc>0;
            xIdc = xIdc(cond);
            yIdc = yIdc(cond);
            xSmFit = x(xIdc); ySmFit = y(yIdc);
            cond = xSmFit > 0;
            xSmFit = xSmFit(cond);
            ySmFit = ySmFit(cond);
            HRBCprcFit = fit(ySmFit',xSmFit','smoothingspline','smoothingparam',smParam);

            xIdc = RRBCKdePrcPos;
            yIdc = 1:length(RRBCKdePrcPos);
            cond = ~isnan(xIdc)&xIdc>0;
            xIdc = xIdc(cond);
            yIdc = yIdc(cond);
            xSmFit = x(xIdc); ySmFit = y(yIdc);
            cond = xSmFit > 0;
            xSmFit = xSmFit(cond);
            ySmFit = ySmFit(cond);
            RRBCprcFit = fit(ySmFit',xSmFit','smoothingspline','smoothingparam',smParam);

            % slope fits
            YInt = linspace(-1,1,100)*20;

            % RRBCs
            interValsY = (y(yOrigin)+YInt)';
            interValsX = (RRBCprcFit(interValsY));

            linfit = fit(interValsX,interValsY,'poly1');
            RRBCSlopeVel = linfit.p1;

            RRBCVel = RRBCSlopeVel*velScaleRRBC;
            bifurcations(bifIdx).conVesBdy(vesIdx).RRBCVel = RRBCVel;


            % HRBCs
            interValsY = (y(yOrigin)+YInt)';
            interValsX = (HRBCprcFit(interValsY));

            linfit = fit(interValsX,interValsY,'poly1');
            HRBCSlopeVel = linfit.p1;

            HRBCVel = HRBCSlopeVel*velScaleHRBC;
            bifurcations(bifIdx).conVesBdy(vesIdx).HRBCVel = HRBCVel;

            
        elseif maximaKde
            %% density maxima as function of length
            HRBCKdeMaxPos = zeros(1,N);
            RRBCKdeMaxPos = zeros(1,N);
            for kdeIdx = 1:N
                HRBCKdeMaxPos(kdeIdx) = find(HRBCkde(kdeIdx,:) == max(HRBCkde(kdeIdx,:)),1,'first');
                RRBCKdeMaxPos(kdeIdx) = find(RRBCkde(kdeIdx,:) == max(RRBCkde(kdeIdx,:)),1,'first');
            end
            %% maximum position fit

            smParam = 0.001;

            xMax = x(HRBCKdeMaxPos); yMax = y(1:length(RRBCKdeMaxPos));
            cond = xMax > 0;
            xMax = xMax(cond);
            yMax = yMax(cond);
            HRBCmaxFit = fit(yMax',xMax','smoothingspline','smoothingparam',smParam);

            xMax = x(RRBCKdeMaxPos); yMax = y(1:length(RRBCKdeMaxPos));
            cond = xMax > 0;
            xMax = xMax(cond);
            yMax = yMax(cond);
            RRBCmaxFit = fit(yMax',xMax','smoothingspline','smoothingparam',smParam);
            %% absolute maximum and slope at maximum
            YInt = linspace(-1,1,100)*5;

            % RRBCs
            [yMax,xMax] = find(RRBCkde == max(RRBCkde(:)),1,'first');
            RRBCMax = [x(xMax),y(yMax)];
            interValsY = (RRBCMax(2)+YInt)';
            interValsX = (RRBCmaxFit(interValsY));

            linfit = fit(interValsX,interValsY,'poly1');
            RRBCSlopeVel = linfit.p1;

            RRBCVel = RRBCSlopeVel*velScaleRRBC;
            bifurcations(bifIdx).conVesBdy(vesIdx).RRBCVel = RRBCVel;


            % HRBCs

            [yMax,xMax] = find(HRBCkde == max(HRBCkde(:)),1,'first');
            HRBCMax = [x(xMax),y(yMax)];

            interValsY = (HRBCMax(2)+YInt)';
            interValsX = (HRBCmaxFit(interValsY));

            linfit = fit(interValsX,interValsY,'poly1');
            HRBCSlopeVel = linfit.p1;

            HRBCVel = HRBCSlopeVel*velScaleHRBC;
            bifurcations(bifIdx).conVesBdy(vesIdx).HRBCVel = HRBCVel;



            % display
            fprintf('-> found velocities HRBC: %.2f, RRBC: %.2f\n',HRBCVel,RRBCVel);
        else
            [yHRBCMax,xHRBCMax] = find(HRBCkde == max(HRBCkde(:)),1,'first');
            [yRRBCMax,xRRBCMax] = find(RRBCkde == max(RRBCkde(:)),1,'first');
            
            HRBCMaxVels = HRBCDist(HRBCTime == round(x(xHRBCMax)))./round(x(xHRBCMax));
            HRBCMaxVels = sort(HRBCMaxVels);

            RRBCMaxVels = RRBCDist(RRBCTime == round(x(xRRBCMax)))./round(x(xRRBCMax));
            RRBCMaxVels = sort(RRBCMaxVels);

            HRBCVel = prctile([HRBCDist],70)./prctile([HRBCTime],30)*velScaleHRBC;
            RRBCVel = prctile([RRBCDist],70)./prctile([RRBCTime],30)*velScaleRRBC;

            bifurcations(bifIdx).conVesBdy(vesIdx).HRBCVel = HRBCVel;
            bifurcations(bifIdx).conVesBdy(vesIdx).RRBCVel = RRBCVel;

        end
        %% plotting
        KdePlots;    % also new velocity estimation from kde

        % display
        fprintf('-> found velocities HRBC: %.2f, RRBC: %.2f\n',HRBCVel,RRBCVel);

        bifurcations(bifIdx).conVesBdy(vesIdx).HRBCVel = HRBCVel;
        bifurcations(bifIdx).conVesBdy(vesIdx).RRBCVel = RRBCVel;
        
        flag_plot_kde_analysis = false;
        if flag_plot_kde_analysis
            close all;
            figure;
            tiledlayout(2,1);
            nexttile;
            hold on;
            surf(x,y,HRBCkde,'EdgeColor','none')
            plot3(HRBCTime,HRBCDist,10000*ones(1,length(HRBCTime)),'.r','MarkerSize',10)
            plot3(xHRBCedge,yHRBCedge,10000*ones(1,length(xHRBCedge)),'.g','MarkerSize',10)
            view(2)
            xlim([minTime,maxTime]);
            ylim([minDist,maxDist]);
            nexttile;
            hold on;
            surf(x,y,RRBCkde,'EdgeColor','none')
            plot3(RRBCTime,RRBCDist,10000*ones(1,length(RRBCTime)),'.r','MarkerSize',10)
            plot3(xRRBCedge,yRRBCedge,10000*ones(1,length(xRRBCedge)),'.g','MarkerSize',10)
            view(2)
            xlim([minTime,maxTime]);
            ylim([minDist,maxDist]);
            return
        end
        catch ME
            disp('An error occurred:');
            disp(ME.message);
            fprintf('warning: kernel density analysis failed for bif %d ves %d.\n',bifIdx,vesIdx);
            bifurcations(bifIdx).conVesBdy(vesIdx).RRBCVel = [nan];
            bifurcations(bifIdx).conVesBdy(vesIdx).HRBCVel = [nan];
        end
        %% beginning of vessel velocities
        % HRBC velocity
        HRBCVelBeg = [];
        for trajIdx = 1:length(HRBCTraj)
            pnts = HRBCTraj(trajIdx).pnts;
            % first filter unsensical trajectories
            trueTrajCond = [0;diff(pnts(:,3))]>=0;
            trueTrajCond = smooth([flip(trueTrajCond);trueTrajCond;flip(trueTrajCond)],3);
            trueTrajCond = trueTrajCond(size(pnts,1)+(1:size(pnts,1)))>0.7;trueTrajCond  = (1:size(pnts,1))>0;
            pnts = pnts(trueTrajCond,:);
            % check which points are inside
            inboundIdc = inpolygon(pnts(:,1),pnts(:,2),vesBounds(:,1),vesBounds(:,2));
            pnts = pnts(inboundIdc,:);

            distances = sqrt((pnts(:,1)-vesStart(1)).^2+(pnts(:,2)-vesStart(2)).^2) - norm(vesStart-medPnt);% distances = distances - min(distances);
            proxIdc = distances<0;
            pnts = real(pnts(proxIdc,:));
            velVec = diff(pnts(:,1:2));
            if length(velVec(:)) > 3
                velVec = velVec./repmat(diff(pnts(:,3)),1,2);
                velMag = sqrt(velVec(:,1).^2+velVec(:,2).^2);
                velMag = velMag(~isnan(velMag) & ~isinf(velMag));
                avgMag = mean(velMag);
                HRBCVelBeg = [HRBCVelBeg,avgMag];
            end
        end
        cond = ~isnan(HRBCVelBeg)&~isinf(HRBCVelBeg);
        HRBCVelBeg = HRBCVelBeg(cond);
        HRBCVelBeg = median(HRBCVelBeg)*velScaleHRBC;
        HRBCVelBegErr = std(HRBCVelBeg)/sqrt(length(HRBCVelBeg));
        % save
        bifurcations(bifIdx).conVesBdy(vesIdx).HRBCVelBeg = HRBCVelBeg;
        bifurcations(bifIdx).conVesBdy(vesIdx).HRBCVelBegErr = HRBCVelBegErr;
        % RRBC velocity
        RRBCVelBeg = [];
        for trajIdx = 1:length(RRBCTraj)
            pnts = RRBCTraj(trajIdx).pnts;
            % first filter unsensical trajectories
            trueTrajCond = [0;diff(pnts(:,3))]>=0;
            trueTrajCond = smooth([flip(trueTrajCond);trueTrajCond;flip(trueTrajCond)],3);
            trueTrajCond = trueTrajCond(size(pnts,1)+(1:size(pnts,1)))>0.7;trueTrajCond  = (1:size(pnts,1))>0;
            pnts = pnts(trueTrajCond,:);
            % check which points are inside
            inboundIdc = inpolygon(pnts(:,1),pnts(:,2),vesBounds(:,1),vesBounds(:,2));
            pnts = pnts(inboundIdc,:);

            distances = sqrt((pnts(:,1)-vesStart(1)).^2+(pnts(:,2)-vesStart(2)).^2) - norm(vesStart-medPnt);% distances = distances - min(distances);
            proxIdc = distances<-0.3*norm(vesStart-medPnt);
            pnts = real(pnts(proxIdc,:));
            velVec = diff(pnts(:,1:2));
            if length(velVec(:)) > 3
                velVec = velVec./repmat(diff(pnts(:,3)),1,2);
                velMag = sqrt(velVec(:,1).^2+velVec(:,2).^2);
                velMag = velMag(~isnan(velMag) & ~isinf(velMag));
                avgMag = mean(velMag);
                RRBCVelBeg = [RRBCVelBeg,avgMag];
            end
        end

        cond = ~isnan(RRBCVelBeg)&~isinf(RRBCVelBeg);
        RRBCVelBeg = RRBCVelBeg(cond);
        RRBCVelBeg = median(RRBCVelBeg)*velScaleRRBC;
        RRBCVelBegErr = std(RRBCVelBeg)/sqrt(length(RRBCVelBeg))*velScaleRRBC;
        % save
        bifurcations(bifIdx).conVesBdy(vesIdx).RRBCVelBeg = RRBCVelBeg;
        bifurcations(bifIdx).conVesBdy(vesIdx).RRBCVelBegErr = RRBCVelBegErr;
        %% end of vessel velocities
        % HRBC velocity
        HRBCVelEnd = [];
        for trajIdx = 1:length(HRBCTraj)
            pnts = HRBCTraj(trajIdx).pnts;
            % first filter unsensical trajectories
            trueTrajCond = [0;diff(pnts(:,3))]>=0;
            trueTrajCond = smooth([flip(trueTrajCond);trueTrajCond;flip(trueTrajCond)],3);
            trueTrajCond = trueTrajCond(size(pnts,1)+(1:size(pnts,1)))>0.7;trueTrajCond  = (1:size(pnts,1))>0;
            pnts = pnts(trueTrajCond,:);
            % check which points are inside
            inboundIdc = inpolygon(pnts(:,1),pnts(:,2),vesBounds(:,1),vesBounds(:,2));
            pnts = pnts(inboundIdc,:);

            distances = sqrt((pnts(:,1)-vesStart(1)).^2+(pnts(:,2)-vesStart(2)).^2) - norm(vesStart-medPnt);% distances = distances - min(distances);
            proxIdc = distances>0.3*norm(vesStart-medPnt);
            pnts = real(pnts(proxIdc,:));
            velVec = diff(pnts(:,1:2));
            if length(velVec(:)) > 3
                velVec = velVec./repmat(diff(pnts(:,3)),1,2);
                velMag = sqrt(velVec(:,1).^2+velVec(:,2).^2);
                velMag = velMag(~isnan(velMag) & ~isinf(velMag));
                avgMag = mean(velMag);
                HRBCVelEnd = [HRBCVelEnd,avgMag];
            end
        end
        cond = ~isnan(HRBCVelEnd)&~isinf(HRBCVelEnd);
        HRBCVelEnd = HRBCVelEnd(cond);
        HRBCVelEnd = median(HRBCVelEnd)*velScaleHRBC;
        HRBCVelEndErr = std(HRBCVelEnd)/sqrt(length(HRBCVelEnd));
        % save
        bifurcations(bifIdx).conVesBdy(vesIdx).HRBCVelEnd = HRBCVelEnd;
        bifurcations(bifIdx).conVesBdy(vesIdx).HRBCVelEndErr = HRBCVelEndErr;

        % RRBC velocity
        RRBCVelEnd = [];
        for trajIdx = 1:length(RRBCTraj)
            pnts = RRBCTraj(trajIdx).pnts;
            % first filter unsensical trajectories
            trueTrajCond = [0;diff(pnts(:,3))]>=0;
            trueTrajCond = smooth([flip(trueTrajCond);trueTrajCond;flip(trueTrajCond)],3);
            trueTrajCond = trueTrajCond(size(pnts,1)+(1:size(pnts,1)))>0.7;trueTrajCond  = (1:size(pnts,1))>0;
            pnts = pnts(trueTrajCond,:);
            % check which points are inside
            inboundIdc = inpolygon(pnts(:,1),pnts(:,2),vesBounds(:,1),vesBounds(:,2));
            pnts = pnts(inboundIdc,:);

            distances = sqrt((pnts(:,1)-vesStart(1)).^2+(pnts(:,2)-vesStart(2)).^2) - norm(vesStart-medPnt);% distances = distances - min(distances);
            proxIdc = distances>0.3*norm(vesStart-medPnt);
            pnts = real(pnts(proxIdc,:));
            velVec = diff(pnts(:,1:2));
            if length(velVec(:)) > 3
                velVec = velVec./repmat(diff(pnts(:,3)),1,2);
                velMag = sqrt(velVec(:,1).^2+velVec(:,2).^2);
                velMag = velMag(~isnan(velMag) & ~isinf(velMag));
                avgMag = mean(velMag);
                RRBCVelEnd = [RRBCVelEnd,avgMag];
            end
        end

        cond = ~isnan(RRBCVelEnd)&~isinf(RRBCVelEnd);
        RRBCVelEnd = RRBCVelEnd(cond);
        RRBCVelEnd = median(RRBCVelEnd)*velScaleRRBC;
        RRBCVelEndErr = std(RRBCVelEnd)/sqrt(length(RRBCVelEnd))*velScaleRRBC;
        % save
        bifurcations(bifIdx).conVesBdy(vesIdx).RRBCVelEnd = RRBCVelEnd;
        bifurcations(bifIdx).conVesBdy(vesIdx).RRBCVelEndErr = RRBCVelEndErr;
        %% plotting / annotation

        % annotate at vessel location
        vesPnts = conVes(vesIdx).pnts;
        cPnt = mean(vesPnts);
        minVal = 50;
        try
            minVal = min(abs(sum(cPnt-vesCtrPnts,2)));
        end
        if 0%minVal<5
            % skip
        else
            vesCtrPnts = [vesCtrPnts;cPnt];
            cAngle = conVes(vesIdx).avgAngle;
            dirVec = [sin(cAngle/180*pi),cos(cAngle/180*pi)];
            cPnt = mean(vesPnts(1:round(size(vesPnts,1)/2),:));
            dPnt = cPnt + 0;
            pos = get(gca, 'Position');
            xArrow = [cPnt(2),dPnt(2)];
            yArrow = size(wallImg,1)-[cPnt(1),dPnt(1)]+10;
            figure(fig);
            annotation('textbox',...
                [(xArrow(1)-min(xlim))/diff(xlim)*pos(3)+pos(1),(yArrow(1)-min(ylim))/diff(ylim)*pos(4)+pos(2),(xArrow(2)-xArrow(1))/diff(xlim)*pos(3),(yArrow(2)-yArrow(1))/diff(ylim)*pos(4)],'String',['' sprintf('%.2f',HRBCVel)],'color',[1 0.2 0.2]);
            annotation('textbox',...
                [(xArrow(1)-min(xlim))/diff(xlim)*pos(3)+pos(1),(yArrow(1)-min(ylim))/diff(ylim)*pos(4)+pos(2),(xArrow(2)-xArrow(1))/diff(xlim)*pos(3),(yArrow(2)-yArrow(1))/diff(ylim)*pos(4)],'String',['     ' sprintf('%.2f',RRBCVel)],'color',[0.2 0.2 1]);
        end
    end
end
hold off
saveVelocities = [allHRBCvel;allRRBCvel;allFlow];
% save bifurcations
save([rootdir,'\bifTraj.mat'],"bifurcations");
% save velocities
save([rootdir,'\velocities.mat'],'saveVelocities');
% print figure
print([rootdir,'\velocities.png'],'-dpng','-r300');