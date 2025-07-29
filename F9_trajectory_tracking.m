%--------------------------------------------------------------------------
% Script Name : F9_trajectory_tracing.m
% Authors     : Felix Maurer, Yazdan Rashidi
% Institution : Saarland University
% Email       : mail@felixmilanmaurer.com
% Date        : 2024
%
% Description :
%   This script traces cell trajectories in the 3D detected points (x,y,z),
%   using a predictive search from the flow estimations of coarse RBC
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
%% File Loop
flag_test = false;
flag_analyze = false;
cellTypes = {'Healthy_RBCs','Rigid_RBCs'};
% find data directories
maskName = 'Mask.png';
rootdir = char(readlines('directory.txt'));
if flag_analyze
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
        if sum(strcmp(fileFolder(1:end-4),folders))>0
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

    % gather files
    fprintf('loading data...\n')
    loadTimer = tic;
    clear cellPoints;
    for fileIdx = 1:length(filelist)
        fileFolder = filelist(fileIdx).folder;
        fileName = filelist(fileIdx).name;
        filePath = [fileFolder '\' fileName];

        % extract cell type
        cellTypeStr = fileFolder(end-2:end);
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
        rbcFrameRateFac = 1/2;
        wbcFrameRateFac = 1;
        % skip frames from WBC
        skipFrames = true;
        if skipFrames
            if cellTypeNum == 2
                allpoints = allpoints([1:2:length(allpoints), 2:2:length(allpoints)]);
            end
            wbcFrameRateFac = 1/2;
        end
        % save
        cellPoints(fileIdx).pnts = allpoints;
        cellPoints(fileIdx).cellType = cellTypeNum;
        cellPoints(fileIdx).alignVec = alignVec;
        cellPoints(fileIdx).filePath = filePath;
        cellPoints(fileIdx).fileFolder = fileFolder;
    end
    fprintf('loading done in %.2f seconds.\n',toc(loadTimer));
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
        load([rootdir,'\bifFlow.mat']); % RBCs
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
                %% go through all bifurcations
                framePause = inf; % seconds
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
                                    if logical(~isnan(prevDir)) && logical(size(prevPnts,1) > 1)
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
                                        xlim(startPntRot(:,1)+[-200,200]);
                                        ylim(startPntRot(:,2)+[-200,200]);
                                        zlim(startPntRot(:,3)+[-200,200]);
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
                                        if logical(~isnan(prevDir)) && logical(size(prevPnts,1) > 1)
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
                                            xlim(startPntRot(:,1)+[-200,200]);
                                            ylim(startPntRot(:,2)+[-200,200]);
                                            zlim(startPntRot(:,3)+[-200,200]);
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
                                    fprintf('status: %.2f percent.\n',(1-size(vesPoints,1)/length(allSelIdc))*100);
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
                            addpath('C:\Users\khadija\Documents\SLMtools\SLMtools');
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
                fprintf('merging all RBC data.\n')
                % merge all clusters of this geometry
                for bifIdx = 1:length(bifurcations)
                    cClu = [];
                    if kBifClu(bifIdx) > 1
                        cClu = clusters(bifIdx).clu;
                        bifurcations(bifIdx).RBCtraj = cClu;
                    else
                        bifurcations(bifIdx).RBCtraj = [];
                    end
                    fprintf('%d clusters in bifurcation %d.\n',length(cClu),bifIdx);
                end
            else
                fprintf('merging all WBC data.\n')
                % merge all clusters of this geometry
                for bifIdx = 1:length(bifurcations)
                    cClu = [];
                    if kBifClu(bifIdx) > 1
                        cClu = clusters(bifIdx).clu;
                        bifurcations(bifIdx).WBCtraj = cClu;
                    else
                        bifurcations(bifIdx).WBCtraj = [];
                    end
                    fprintf('%d clusters in bifurcation %d.\n',length(cClu),bifIdx);
                end
            end
        end
        %% determine bifurcation average velocity magnitude
        for bifIdx = 1:length(bifurcations)
            bifurcations(bifIdx).modeThres = 0;
            WBCtraj = bifurcations(bifIdx).WBCtraj;
            allWBCvel = [];
            for trajIdx = 1:length(WBCtraj)
                traj = WBCtraj(trajIdx);
                if ~isempty(traj.pnts)
                    neigh = traj.neigh;
                    pnts =  traj.pnts;
                    velVec = diff(pnts(:,1:2));
                    if length(velVec(:)) > 3
                        velVec = velVec./repmat(diff(pnts(:,3)),1,2);
                        velMag = sqrt(velVec(:,1).^2+velVec(:,2).^2);
                        velMag = velMag(~isnan(velMag) & ~isinf(velMag));
                        avgMag = mean(velMag);
                        allWBCvel = [allWBCvel;avgMag];
                    end
                end
            end
            allWBCvel = allWBCvel(allWBCvel > 0);

            % healthy
            RBCtraj = bifurcations(bifIdx).RBCtraj;
            allRBCvel = [];
            for trajIdx = 1:length(RBCtraj)
                traj = RBCtraj(trajIdx);
                if ~isempty(traj.pnts)
                    neigh = traj.neigh;
                    pnts =  traj.pnts;
                    velVec = diff(pnts(:,1:2));
                    if length(velVec(:)) > 3
                        velVec = velVec./repmat(diff(pnts(:,3)),1,2);
                        velMag = sqrt(velVec(:,1).^2+velVec(:,2).^2);
                        velMag = velMag(~isnan(velMag) & ~isinf(velMag));
                        avgMag = mean(velMag);
                        allRBCvel = [allRBCvel;avgMag];
                    end
                end
            end
            allRBCvel = allRBCvel(allRBCvel > 0);

            bifurcations(bifIdx).avgRBCvel = median(allRBCvel);
            bifurcations(bifIdx).avgWBCvel = median(allWBCvel);

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
velScaleRBC = velScale*rbcFrameRateFac;
velScaleWBC = velScale*wbcFrameRateFac;
velScaleFlow = velScale*rbcFrameRateFac;
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
allRBCvel = [];
allWBCvel = [];
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
    RBCVel = bifurcations(bifIdx).avgRBCvel*velScaleRBC;
    WBCVel = bifurcations(bifIdx).avgWBCvel*velScaleWBC;
    RBCTraj = bifurcations(bifIdx).RBCtraj;
    WBCTraj = bifurcations(bifIdx).WBCtraj;
    % add vessels
    conVes = bifurcations(bifIdx).conVesBdy;
    bifurcations(bifIdx).velScaleRBC = velScaleRBC;
    bifurcations(bifIdx).velScaleWBC = velScaleWBC;
    
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
        % RBC velocity
        RBCVel = [];
        RBCDist = [];
        RBCTime = [];
        clear RBCVesTraj; kTraj = 1;
        for trajIdx = 1:length(RBCTraj)
            pnts = RBCTraj(trajIdx).pnts;
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
                RBCVesTraj(kTraj).pnts = pnts;kTraj = kTraj + 1;
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

                RBCDist = [RBCDist,ds];
                RBCTime = [RBCTime,dt];

                velMag = maxTrajDistance/maxTrajDistanceTime;
                avgMag = mean(velMag);
                if imag(complex(avgMag))>0
                    return
                end
                RBCVel = [RBCVel,avgMag];
            end
        end
        if ~exist('RBCVesTraj','var')
            RBCVesTraj = [];
        end
        cond = ~isnan(RBCVel)&~isinf(RBCVel);
        RBCVel = RBCVel(cond);
        RBCNum = length(RBCVel);
        RBCVelArray = RBCVel;
        RBCVelErr = std(RBCVel)/sqrt(length(RBCVel))*velScaleRBC;
        % save
        bifurcations(bifIdx).conVesBdy(vesIdx).RBCVel = RBCVel;
        bifurcations(bifIdx).conVesBdy(vesIdx).RBCNum = RBCNum;
        bifurcations(bifIdx).conVesBdy(vesIdx).RBCVelArray = RBCVelArray;
        bifurcations(bifIdx).conVesBdy(vesIdx).RBCVelErr = RBCVelErr;
        bifurcations(bifIdx).conVesBdy(vesIdx).RBCVesTraj = RBCVesTraj;

        % WBC velocity
        WBCVel = [];
        WBCDist = [];
        WBCTime = [];
        clear WBCVesTraj; kTraj = 1;
        for trajIdx = 1:length(WBCTraj)
            pnts = WBCTraj(trajIdx).pnts;
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
                WBCVesTraj(kTraj).pnts = pnts;kTraj = kTraj + 1;
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

                WBCDist = [WBCDist,ds];
                WBCTime = [WBCTime,dt];

                velMag = maxTrajDistance/maxTrajDistanceTime;
                avgMag = mean(velMag);
                WBCVel = [WBCVel,avgMag];
            end
        end
        if ~exist('WBCVesTraj','var')
            WBCVesTraj = [];
        end
        cond = ~isnan(WBCVel)&~isinf(WBCVel);
        WBCVel = WBCVel(cond);
        WBCNum = length(WBCVel);
        WBCVelArray = WBCVel;
        WBCVelErr = std(WBCVel)/sqrt(length(WBCVel))*velScaleWBC;
        % save
        bifurcations(bifIdx).conVesBdy(vesIdx).WBCVel = WBCVel;
        bifurcations(bifIdx).conVesBdy(vesIdx).WBCNum = WBCNum;
        bifurcations(bifIdx).conVesBdy(vesIdx).WBCVelArray = WBCVelArray;
        bifurcations(bifIdx).conVesBdy(vesIdx).WBCVelErr = WBCVelErr;
        % lengths and times
        bifurcations(bifIdx).conVesBdy(vesIdx).RBCDist = RBCDist;
        bifurcations(bifIdx).conVesBdy(vesIdx).RBCTime = RBCTime;
        bifurcations(bifIdx).conVesBdy(vesIdx).WBCDist = WBCDist;
        bifurcations(bifIdx).conVesBdy(vesIdx).WBCTime = WBCTime;
        bifurcations(bifIdx).conVesBdy(vesIdx).WBCVesTraj = WBCVesTraj;
        %% kernel density estimation
        try
        fprintf('starting kernel density analysis for bif %d ves %d.\n',bifIdx,vesIdx);
        % optionally estimate h value
        %VelVals = RBCTime;
        %[~,~,bandwidth] = ksdensity(VelVals(~isnan(VelVals)&~isinf(VelVals)&VelVals>0));%,'support','positive','boundarycorrection','reflection');
        %h = max([1,round(bandwidth)]);
        h = 6;
        N = 101;
        maxTime = max([RBCTime,WBCTime]) + 5;
        minTime = 0;
        maxDist = max([RBCDist,WBCDist]) + 20;
        minDist = 0;
        x = linspace(minTime,maxTime,N);
        y = linspace(minDist,maxDist,N);
        xScale = median(diff(x));
        yScale = median(diff(y));

        % RBC
        dataPnts = [RBCTime'-minTime,RBCDist'-minDist];
        dataPnts = dataPnts(RBCDist>max([prctile(RBCDist,50),prctile(WBCDist,50)])-max(std(RBCDist),std(WBCDist)),:);
        

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

        RBCkde = kde2d';

        % WBC
        dataPnts = [WBCTime'-minTime,WBCDist'-minDist];
        dataPnts = dataPnts(WBCDist>max([prctile(RBCDist,50),prctile(WBCDist,50)])-max(std(RBCDist),std(WBCDist)),:);

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

        WBCkde = kde2d';
        %% distribution edges
        edgeThres = mean(0.5*mean(RBCkde(RBCkde>0)));
        bwKdeImg = RBCkde>edgeThres;
        bwKdeImg(1,:) = 0;        bwKdeImg(end,:) = 0;        bwKdeImg(:,1) = 0;        bwKdeImg(:,end) = 0;
        RBCedge = edge(bwKdeImg,'roberts');
        [yRBCedge,xRBCedge] = find(RBCedge);
        xRBCedge = x(xRBCedge);
        yRBCedge = y(yRBCedge);

        edgeThres = mean(0.5*mean(WBCkde(WBCkde>0)));
        bwKdeImg = WBCkde>edgeThres;
        bwKdeImg(1,:) = 0;        bwKdeImg(end,:) = 0;        bwKdeImg(:,1) = 0;        bwKdeImg(:,end) = 0;
        WBCedge = edge(bwKdeImg,'roberts');
        [yWBCedge,xWBCedge] = find(WBCedge);
        xWBCedge = x(xWBCedge);
        yWBCedge = y(yWBCedge);

        computeByPerc = false;
        maximaKde = false;
        if computeByPerc
            %% density percentile as function of length
            % maximum positions
            RBCKdeMaxPos = zeros(1,N);
            WBCKdeMaxPos = zeros(1,N);
            for kdeIdx = 1:N
                RBCKdeMaxPos(kdeIdx) = find(RBCkde(kdeIdx,:) == max(RBCkde(kdeIdx,:)),1,'first');
                WBCKdeMaxPos(kdeIdx) = find(WBCkde(kdeIdx,:) == max(WBCkde(kdeIdx,:)),1,'first');
            end
            % origin speed
            [yOrigin,xOrigin] = find(WBCkde == max(WBCkde(:)),1,'first');
            % corresponding percentiles
            RBCpercentile = sum(RBCkde(yOrigin,1:xOrigin))/sum(RBCkde(yOrigin,:))*100; % in percent
            WBCpercentile = sum(WBCkde(yOrigin,1:xOrigin))/sum(WBCkde(yOrigin,:))*100; % in percent
            % all percentile positions
            RBCKdePrcPos = zeros(1,N);
            WBCKdePrcPos = zeros(1,N);
            for kdeIdx = 1:N
                % RBC
                RBCcdf = cumsum(RBCkde(kdeIdx,:))/sum(RBCkde(kdeIdx,:))*100; % percent
                prctilePos = find(RBCcdf > RBCpercentile,1,'first');
                if isempty(prctilePos)
                    prctilePos = 0;
                end
                RBCKdePrcPos(kdeIdx) = prctilePos;
                % WBC
                WBCcdf = cumsum(WBCkde(kdeIdx,:))/sum(WBCkde(kdeIdx,:))*100; % percent
                prctilePos = find(WBCcdf > WBCpercentile,1,'first');
                if isempty(prctilePos)
                    prctilePos = 0;
                end
                WBCKdePrcPos(kdeIdx) = prctilePos;
            end
            RBCKdePrcPos(RBCKdePrcPos==0) = nan;
            WBCKdePrcPos(RBCKdePrcPos==0) = nan;

            flag_plotPercentiles = false;
            if flag_plotPercentiles
                close all;
                figure
                hold on
                plot(RBCKdeMaxPos)
                plot(RBCKdePrcPos)
                plot(WBCKdeMaxPos)
                plot(WBCKdePrcPos)
                xline(yOrigin)
                hold off

                legend('RBCMax','RBCPrc','WBCMax','WBCPrc');
                return
            end

            %% percentile position fit

            smParam = 0.001;

            xIdc = RBCKdePrcPos;
            yIdc = 1:length(RBCKdePrcPos);
            cond = ~isnan(xIdc)&xIdc>0;
            xIdc = xIdc(cond);
            yIdc = yIdc(cond);
            xSmFit = x(xIdc); ySmFit = y(yIdc);
            cond = xSmFit > 0;
            xSmFit = xSmFit(cond);
            ySmFit = ySmFit(cond);
            RBCprcFit = fit(ySmFit',xSmFit','smoothingspline','smoothingparam',smParam);

            xIdc = WBCKdePrcPos;
            yIdc = 1:length(WBCKdePrcPos);
            cond = ~isnan(xIdc)&xIdc>0;
            xIdc = xIdc(cond);
            yIdc = yIdc(cond);
            xSmFit = x(xIdc); ySmFit = y(yIdc);
            cond = xSmFit > 0;
            xSmFit = xSmFit(cond);
            ySmFit = ySmFit(cond);
            WBCprcFit = fit(ySmFit',xSmFit','smoothingspline','smoothingparam',smParam);

            % slope fits
            YInt = linspace(-1,1,100)*20;

            % WBCs
            interValsY = (y(yOrigin)+YInt)';
            interValsX = (WBCprcFit(interValsY));

            linfit = fit(interValsX,interValsY,'poly1');
            WBCSlopeVel = linfit.p1;

            WBCVel = WBCSlopeVel*velScaleWBC;
            bifurcations(bifIdx).conVesBdy(vesIdx).WBCVel = WBCVel;


            % RBCs
            interValsY = (y(yOrigin)+YInt)';
            interValsX = (RBCprcFit(interValsY));

            linfit = fit(interValsX,interValsY,'poly1');
            RBCSlopeVel = linfit.p1;

            RBCVel = RBCSlopeVel*velScaleRBC;
            bifurcations(bifIdx).conVesBdy(vesIdx).RBCVel = RBCVel;

            
        elseif maximaKde
            %% density maxima as function of length
            RBCKdeMaxPos = zeros(1,N);
            WBCKdeMaxPos = zeros(1,N);
            for kdeIdx = 1:N
                RBCKdeMaxPos(kdeIdx) = find(RBCkde(kdeIdx,:) == max(RBCkde(kdeIdx,:)),1,'first');
                WBCKdeMaxPos(kdeIdx) = find(WBCkde(kdeIdx,:) == max(WBCkde(kdeIdx,:)),1,'first');
            end
            %% maximum position fit

            smParam = 0.001;

            xMax = x(RBCKdeMaxPos); yMax = y(1:length(WBCKdeMaxPos));
            cond = xMax > 0;
            xMax = xMax(cond);
            yMax = yMax(cond);
            RBCmaxFit = fit(yMax',xMax','smoothingspline','smoothingparam',smParam);

            xMax = x(WBCKdeMaxPos); yMax = y(1:length(WBCKdeMaxPos));
            cond = xMax > 0;
            xMax = xMax(cond);
            yMax = yMax(cond);
            WBCmaxFit = fit(yMax',xMax','smoothingspline','smoothingparam',smParam);
            %% absolute maximum and slope at maximum
            YInt = linspace(-1,1,100)*5;

            % WBCs
            [yMax,xMax] = find(WBCkde == max(WBCkde(:)),1,'first');
            WBCMax = [x(xMax),y(yMax)];
            interValsY = (WBCMax(2)+YInt)';
            interValsX = (WBCmaxFit(interValsY));

            linfit = fit(interValsX,interValsY,'poly1');
            WBCSlopeVel = linfit.p1;

            WBCVel = WBCSlopeVel*velScaleWBC;
            bifurcations(bifIdx).conVesBdy(vesIdx).WBCVel = WBCVel;


            % RBCs

            [yMax,xMax] = find(RBCkde == max(RBCkde(:)),1,'first');
            RBCMax = [x(xMax),y(yMax)];

            interValsY = (RBCMax(2)+YInt)';
            interValsX = (RBCmaxFit(interValsY));

            linfit = fit(interValsX,interValsY,'poly1');
            RBCSlopeVel = linfit.p1;

            RBCVel = RBCSlopeVel*velScaleRBC;
            bifurcations(bifIdx).conVesBdy(vesIdx).RBCVel = RBCVel;



            % display
            fprintf('-> found velocities RBC: %.2f, WBC: %.2f\n',RBCVel,WBCVel);
        else
            [yRBCMax,xRBCMax] = find(RBCkde == max(RBCkde(:)),1,'first');
            [yWBCMax,xWBCMax] = find(WBCkde == max(WBCkde(:)),1,'first');
            
            RBCMaxVels = RBCDist(RBCTime == round(x(xRBCMax)))./round(x(xRBCMax));
            RBCMaxVels = sort(RBCMaxVels);

            WBCMaxVels = WBCDist(WBCTime == round(x(xWBCMax)))./round(x(xWBCMax));
            WBCMaxVels = sort(WBCMaxVels);

            RBCVel = prctile([RBCDist],70)./prctile([RBCTime],30)*velScaleRBC;
            WBCVel = prctile([WBCDist],70)./prctile([WBCTime],30)*velScaleWBC;

            bifurcations(bifIdx).conVesBdy(vesIdx).RBCVel = RBCVel;
            bifurcations(bifIdx).conVesBdy(vesIdx).WBCVel = WBCVel;

        end
        %% plotting
        KdePlots;    % also new velocity estimation from kde

        % display
        fprintf('-> found velocities RBC: %.2f, WBC: %.2f\n',RBCVel,WBCVel);

        bifurcations(bifIdx).conVesBdy(vesIdx).RBCVel = RBCVel;
        bifurcations(bifIdx).conVesBdy(vesIdx).WBCVel = WBCVel;
        
        flag_plot_kde_analysis = false;
        if flag_plot_kde_analysis
            close all;
            figure;
            tiledlayout(2,1);
            nexttile;
            hold on;
            surf(x,y,RBCkde,'EdgeColor','none')
            plot3(RBCTime,RBCDist,10000*ones(1,length(RBCTime)),'.r','MarkerSize',10)
            plot3(xRBCedge,yRBCedge,10000*ones(1,length(xRBCedge)),'.g','MarkerSize',10)
            view(2)
            xlim([minTime,maxTime]);
            ylim([minDist,maxDist]);
            nexttile;
            hold on;
            surf(x,y,WBCkde,'EdgeColor','none')
            plot3(WBCTime,WBCDist,10000*ones(1,length(WBCTime)),'.r','MarkerSize',10)
            plot3(xWBCedge,yWBCedge,10000*ones(1,length(xWBCedge)),'.g','MarkerSize',10)
            view(2)
            xlim([minTime,maxTime]);
            ylim([minDist,maxDist]);
            return
        end
        catch
            fprintf('warning: kernel density analysis failed for bif %d ves %d.\n',bifIdx,vesIdx);
            bifurcations(bifIdx).conVesBdy(vesIdx).WBCVel = [nan];
            bifurcations(bifIdx).conVesBdy(vesIdx).RBCVel = [nan];
        end
        %% beginning of vessel velocities
        % RBC velocity
        RBCVelBeg = [];
        for trajIdx = 1:length(RBCTraj)
            pnts = RBCTraj(trajIdx).pnts;
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
                RBCVelBeg = [RBCVelBeg,avgMag];
            end
        end
        cond = ~isnan(RBCVelBeg)&~isinf(RBCVelBeg);
        RBCVelBeg = RBCVelBeg(cond);
        RBCVelBeg = median(RBCVelBeg)*velScaleRBC;
        RBCVelBegErr = std(RBCVelBeg)/sqrt(length(RBCVelBeg));
        % save
        bifurcations(bifIdx).conVesBdy(vesIdx).RBCVelBeg = RBCVelBeg;
        bifurcations(bifIdx).conVesBdy(vesIdx).RBCVelBegErr = RBCVelBegErr;
        % WBC velocity
        WBCVelBeg = [];
        for trajIdx = 1:length(WBCTraj)
            pnts = WBCTraj(trajIdx).pnts;
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
                WBCVelBeg = [WBCVelBeg,avgMag];
            end
        end

        cond = ~isnan(WBCVelBeg)&~isinf(WBCVelBeg);
        WBCVelBeg = WBCVelBeg(cond);
        WBCVelBeg = median(WBCVelBeg)*velScaleWBC;
        WBCVelBegErr = std(WBCVelBeg)/sqrt(length(WBCVelBeg))*velScaleWBC;
        % save
        bifurcations(bifIdx).conVesBdy(vesIdx).WBCVelBeg = WBCVelBeg;
        bifurcations(bifIdx).conVesBdy(vesIdx).WBCVelBegErr = WBCVelBegErr;
        %% end of vessel velocities
        % RBC velocity
        RBCVelEnd = [];
        for trajIdx = 1:length(RBCTraj)
            pnts = RBCTraj(trajIdx).pnts;
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
                RBCVelEnd = [RBCVelEnd,avgMag];
            end
        end
        cond = ~isnan(RBCVelEnd)&~isinf(RBCVelEnd);
        RBCVelEnd = RBCVelEnd(cond);
        RBCVelEnd = median(RBCVelEnd)*velScaleRBC;
        RBCVelEndErr = std(RBCVelEnd)/sqrt(length(RBCVelEnd));
        % save
        bifurcations(bifIdx).conVesBdy(vesIdx).RBCVelEnd = RBCVelEnd;
        bifurcations(bifIdx).conVesBdy(vesIdx).RBCVelEndErr = RBCVelEndErr;

        % WBC velocity
        WBCVelEnd = [];
        for trajIdx = 1:length(WBCTraj)
            pnts = WBCTraj(trajIdx).pnts;
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
                WBCVelEnd = [WBCVelEnd,avgMag];
            end
        end

        cond = ~isnan(WBCVelEnd)&~isinf(WBCVelEnd);
        WBCVelEnd = WBCVelEnd(cond);
        WBCVelEnd = median(WBCVelEnd)*velScaleWBC;
        WBCVelEndErr = std(WBCVelEnd)/sqrt(length(WBCVelEnd))*velScaleWBC;
        % save
        bifurcations(bifIdx).conVesBdy(vesIdx).WBCVelEnd = WBCVelEnd;
        bifurcations(bifIdx).conVesBdy(vesIdx).WBCVelEndErr = WBCVelEndErr;
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
                [(xArrow(1)-min(xlim))/diff(xlim)*pos(3)+pos(1),(yArrow(1)-min(ylim))/diff(ylim)*pos(4)+pos(2),(xArrow(2)-xArrow(1))/diff(xlim)*pos(3),(yArrow(2)-yArrow(1))/diff(ylim)*pos(4)],'String',['' sprintf('%.2f',RBCVel)],'color',[1 0.2 0.2]);
            annotation('textbox',...
                [(xArrow(1)-min(xlim))/diff(xlim)*pos(3)+pos(1),(yArrow(1)-min(ylim))/diff(ylim)*pos(4)+pos(2),(xArrow(2)-xArrow(1))/diff(xlim)*pos(3),(yArrow(2)-yArrow(1))/diff(ylim)*pos(4)],'String',['     ' sprintf('%.2f',WBCVel)],'color',[0.2 0.2 1]);
        end
    end
end
hold off
saveVelocities = [allRBCvel;allWBCvel;allFlow];
% save bifurcations
save([rootdir,'\bifTraj.mat'],"bifurcations");
% save velocities
save([rootdir,'\velocities.mat'],'saveVelocities');
% print figure
print([rootdir,'\velocities.png'],'-dpng','-r300');