%--------------------------------------------------------------------------
% Script Name : F10_lingering_analysis.m
% Authors     : Felix Maurer, Yazdan Rashidi
% Institution : Saarland University
% Email       : mail@felixmilanmaurer.com
% Date        : 2024
%
% Description :
%   This script combines measurements from all geometries to analyze
%   lingering and lateral migration, and speed statistics.
%
% Usage :
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
flag_plotting = false;
%% source
addpath('src');
%% file loop
% analyze or skip and load previous analysis
flag_analyze = true;
if flag_analyze
    cellTypes = {'Healthy_RBCs','Rigid_RBCs'};
    % find data directories
    maskName = 'Mask.png';
    parentDir = char(readlines('directory.txt'));
    allDirList = dir(fullfile([parentDir,'\**\Mask.png']));

    % mother velocities
    allHRBCVelMoth = [];
    allRRBCVelMoth = [];
    allHRBCVelMothEnd = [];
    allRRBCVelMothEnd = [];
    allHRBCVelMothBeg = [];
    allRRBCVelMothBeg = [];
    ErrallHRBCVelMoth = [];
    ErrallRRBCVelMoth = [];
    ErrallHRBCVelMothEnd = [];
    ErrallRRBCVelMothEnd = [];
    ErrallHRBCVelMothBeg = [];
    ErrallRRBCVelMothBeg = [];
    allHRBCNumMoth = [];
    allRRBCNumMoth = [];

    % daughter velocities
    allHRBCVelDaugh = [];
    allRRBCVelDaugh = [];
    allHRBCVelDaughEnd = [];
    allRRBCVelDaughEnd = [];
    allHRBCVelDaughBeg = [];
    allRRBCVelDaughBeg = [];
    ErrallHRBCVelDaugh = [];
    ErrallRRBCVelDaugh = [];
    ErrallHRBCVelDaughEnd = [];
    ErrallRRBCVelDaughEnd = [];
    ErrallHRBCVelDaughBeg = [];
    ErrallRRBCVelDaughBeg = [];
    allHRBCNumDaugh = [];
    allRRBCNumDaugh = [];

    % bifurcation velocities
    allHRBCVelBif = [];
    allRRBCVelBif = [];
    ErrallHRBCVelBif = [];
    ErrallRRBCVelBif = [];
    allHRBCNumBif = [];
    allRRBCNumBif = [];

    % lingering
    allHRBCLin = [];
    allRRBCLin = [];
    ErrallHRBCLin = [];
    ErrallRRBCLin = [];

    % geometry
    allMothLens = []; kEntries = 1;

    bifNum = 0;
    for dirIdx = 1:length(allDirList)
        rootDir = allDirList(dirIdx).folder;
        filelist = dir(fullfile([rootDir,'\**\*bifTraj.mat']));
        fileFolder = filelist(1).folder;
        fileName = filelist(1).name;
        filePath = [fileFolder,'\',fileName];
        %%
        load(filePath);
        for bifIdx = 1:length(bifurcations)
            if bifIdx > 1
                bifBW = bifBW+bifurcations(bifIdx).bifBW;
            else
                bifBW = bifurcations(bifIdx).bifBW;
            end
        end
        % get all trajectory clusters
        clear allClusters;
        clusters_from_network_tracking = true;
        if clusters_from_network_tracking
            frameRates = [50,100];
            cellTypeIdc = [1,2];
            for cellTypeIdx = cellTypeIdc
                filelist = dir(fullfile([rootDir,'\**\' cellTypes{cellTypeIdx} '\ROI*merge.mat']));
                kClu = 1;
                for fileIdx = 1:length(filelist)
                    fileFolder = filelist(fileIdx).folder;
                    fileName = filelist(fileIdx).name;
                    filePath = [fileFolder,'\',fileName];
                    load(filePath)
                    for cluIdx = 1:length(clu)
                        pnts = clu(cluIdx).points;
                        % optional frame skipping to equalize framerate
                        % if cellTypeIdx == 2 % skipping
                        %     pnts = pnts(1:2:end,:);
                        % end
                        allClusters(cellTypeIdx).clu(kClu).pnts = pnts; kClu = kClu + 1;
                        % plot for debugging
                        %plot(pnts(:,2),pnts(:,1),'.-','MarkerSize',20);
                    end
                end
            end
        else % take clusters from bifurcation struct
            frameRates = [50,50]; % skipping for RRBC means 50
            kClu1 = 1;
            kClu2 = 1;
            for bifIdx = 1:length(bifurcations)
                % HRBCs
                HRBCtraj = bifurcations(bifIdx).HRBCtraj;
                for trajIdx = 1:length(HRBCtraj)
                    pnts = HRBCtraj(trajIdx).pnts;
                    allClusters(1).clu(kClu1).pnts = pnts; kClu1 = kClu1 + 1;
                end
                % RRBCs
                RRBCtraj = bifurcations(bifIdx).RRBCtraj;
                for trajIdx = 1:length(RRBCtraj)
                    pnts = RRBCtraj(trajIdx).pnts;
                    allClusters(2).clu(kClu2).pnts = pnts; kClu2 = kClu2 + 1;
                end
            end
        end
        bifNum = bifNum+length(bifurcations);
        totVesBW = bifurcations(1).bifBW;
        for bifIdx = 1:length(bifurcations)
            totVesBW = totVesBW | bifurcations(bifIdx).bifBW;
            cVes = bifurcations(bifIdx).conVesBdy;
            for vesIdx = 1:length(cVes)
                totVesBW = totVesBW | cVes(vesIdx).vesBW;
            end
        end
        imshow(totVesBW); hold on;
        for bifIdx = 1:length(bifurcations)
            try
                bifPnts = bifurcations(bifIdx).bifInnerPnts;
                bifBW = bifurcations(bifIdx).bifBW;
                bifBound = boundary(bifPnts(:,1),bifPnts(:,2));
                bifBoundPnts = bifPnts(bifBound,:);
                cBif = bifurcations(bifIdx);
                bifAngle = cBif.bifAvgAngle;
                cVes = cBif.conVesBdy;
                angleDif = zeros(1,length(cVes));
                for vesIdx = 1:length(cVes)
                    %hold on
                    vesStartPnts = flip(cVes(vesIdx).pnts(1:5,:),1);
                    dirVec = mean(diff(vesStartPnts));
                    vesAngle = -acos(dirVec(:,2)).*(1-2*(-dirVec(:,1)<0))/pi*180;
                    angleDif(vesIdx) = abs(bifAngle-vesAngle);
                end
                inVesIdx = find(angleDif == min(angleDif),1,'first');
                inVesBW = cVes(inVesIdx).vesBW;
                [x,y] = find(inVesBW);
                inVesBound = boundary(x,y);
                                
                % mother velocity
                HRBCVelMoth = cVes(inVesIdx).HRBCVel;
                RRBCVelMoth = cVes(inVesIdx).RRBCVel;
                
                % mother end velocity
                HRBCVelMothEnd= cVes(inVesIdx).HRBCVelEnd;
                RRBCVelMothEnd = cVes(inVesIdx).RRBCVelEnd;

                % mother beginning velocity
                HRBCVelMothBeg= cVes(inVesIdx).HRBCVelBeg;
                RRBCVelMothBeg = cVes(inVesIdx).RRBCVelBeg;

                % mother numbers
                HRBCNumMoth = cVes(inVesIdx).HRBCNum;
                RRBCNumMoth = cVes(inVesIdx).RRBCNum;

                % compute velocity error mother
                HRBCVelMothErr = std(cVes(inVesIdx).HRBCVelArray)*bifurcations(1).velScaleHRBC;
                RRBCVelMothErr = std(cVes(inVesIdx).RRBCVelArray)*bifurcations(1).velScaleRRBC;
                
                % find velocity in daughter branches
                selDaugh = ones(1,length(cVes),'logical');
                selDaugh(inVesIdx) = 0;

                %% finding the apex
                totBifBW = bifBW;
                for vesIdx = 1:length(cVes)
                    totBifBW = totBifBW + cVes(vesIdx).vesBW;
                end
                wallBW = edge(totBifBW);
                wallBW = bwskel(wallBW);

                bifCtr = cBif.bifCtr;
                [x,y] = find(wallBW);
                distancesToBifCtr = sqrt((x-bifCtr(:,1)).^2+(y-bifCtr(:,2)).^2);
                selIdc = distancesToBifCtr < 30;
                x = x(selIdc);
                y = y(selIdc);
                wallBW = 0*wallBW;
                for idx = 1:length(x)
                    wallBW(x(idx),y(idx)) = 1;
                end
                %imshow(wallBW)
                wallBW = logical(wallBW);

                fprintf('--> analyzing bifurcation %d\n',bifIdx);
                bdy = bwboundaries(wallBW);
                pnts1 = bdy{1};
                pnts2 = bdy{2};
                pnts3 = bdy{3};
                idc = zeros(size(pnts1,1)*size(pnts2,1)*size(pnts3,1),3);
                circum = zeros(1,size(idc,1));
                kTri = 1;
                timer = tic;
                for id1 = 1:size(pnts1,1)
                    for id2 = 1:size(pnts2,1)
                        for id3 = 1:size(pnts3,1)
                            % find triangle with lowest circumference
                            dist1 = norm(pnts1(id1,:)-pnts2(id2,:));
                            dist2 = norm(pnts1(id1,:)-pnts3(id3,:));
                            dist3 = norm(pnts2(id2,:)-pnts3(id3,:));
                            idc(kTri,:) = [id1,id2,id3];
                            circum(kTri) = dist1 + dist2 + dist3;
                            kTri = kTri+1;
                        end
                    end
                    % status display
                    if ~mod(id1-1,size(pnts1,1)/10)
                        time = toc(timer);
                        fprintf('status:%.2f percent.\n',id1/size(pnts1,1)*100);
                        fprintf('remaining:%.2f seconds.\n',time/id1*(size(pnts1,1)-id1));
                    end
                end
                % find minimum
                minIdx = find(circum==min(circum),1,'first');
                bifPnts = [pnts1(idc(minIdx,1),:); pnts2(idc(minIdx,2),:); pnts3(idc(minIdx,3),:)];
                bifurcations(bifIdx).bifPnts = bifPnts;
                bifurcations(bifIdx).bifPntsIdc = idc(minIdx,:);
                

                % find flow apex
                mothPnts = cVes(inVesIdx).pnts;
                mothPnt = mothPnts(1,:);
                distances = sqrt(((mothPnt(1)-bifPnts(:,1)).^2+(mothPnt(2)-bifPnts(:,2)).^2));
                % apex search
                mothDir = -mean(diff(mothPnts(1:4,:)));
                expApex = mothPnt+mothDir/norm(mothDir)*max(distances);
                apexDist = sqrt(((expApex(1)-bifPnts(:,1)).^2+(expApex(2)-bifPnts(:,2)).^2));
                apexIdx = find(apexDist==min(apexDist),1,'first');
                apex = bifPnts(apexIdx,:);
                
                noApexIdc = logical([1 1 1]);
                noApexIdc(apexIdx) = 0;
                noApexPnts = bifPnts(noApexIdc,:);
                ApexShiftVec = noApexPnts-apex;
                ApexShiftVec = mean(ApexShiftVec);
                ApexShiftVec = ApexShiftVec./norm(ApexShiftVec);
                shiftedApex = apex + ApexShiftVec*4.6; % pixel shift
                
                % go through clusters
                clear bifFrameNum;
                clear bifTrajLen;
                for cellTypeIdx = 1:2
                    cClu = allClusters(cellTypeIdx).clu;
                    for cluIdx = 1:length(cClu)
                        pnts = cClu(cluIdx).pnts;
                        selIdc = inpolygon(pnts(:,1),pnts(:,2),bifBoundPnts(:,1),bifBoundPnts(:,2));
                        bifFrameNum(cellTypeIdx).clu(cluIdx).frames = sum(selIdc);
                        % compute distance
                        if sum(selIdc) <2
                            bifTrajLen(cellTypeIdx).clu(cluIdx).distance = 0;
                        else
                            dPnts = diff(pnts(selIdc,:));
                            trajLength = sum(sqrt(dPnts(:,1).^2+dPnts(:,2).^2));
                            bifTrajLen(cellTypeIdx).clu(cluIdx).distance = trajLength;
                        end
                    end
                end
                % redefine apex by trajectories
                clear trajApex;
                for cellTypeIdx = 1:2
                    allPnts = [];
                    cClu = allClusters(cellTypeIdx).clu;
                    for cluIdx = 1:length(cClu)
                        pnts = cClu(cluIdx).pnts;
                        bifFrames = bifFrameNum(cellTypeIdx).clu(cluIdx).frames;
                        if bifFrames > 2
                            allPnts = [allPnts;pnts];
                        end
                    end
                    if ~isempty(allPnts)
                        allApexDistances = sqrt((apex(1)-allPnts(:,1)).^2+(apex(2)-allPnts(:,2)).^2);
                        minIdx = find(allApexDistances == min(allApexDistances),1,'first');
                        trajApex(cellTypeIdx).apex = allPnts(minIdx,1:2);
                    else
                        trajApex(cellTypeIdx).apex = apex;
                    end
                end
                % measure apex distance
                clear distApex;
                for cellTypeIdx = 1:2
                    cClu = allClusters(cellTypeIdx).clu;
                    cApex = shiftedApex;
                    for cluIdx = 1:length(cClu)
                        pnts = cClu(cluIdx).pnts;
                        apexDistances = sqrt((cApex(1)-pnts(:,1)).^2+(cApex(2)-pnts(:,2)).^2);
                        distApex(cellTypeIdx).clu(cluIdx).minApexDist = min(apexDistances);
                    end
                end

                LenRefPixels = norm(apex-mothPnt);
                
                plot(mothPnts(:,2),mothPnts(:,1),'-g')
                plot(mothPnt(2),mothPnt(1),'.r');
                plot(bifPnts(:,2),bifPnts(:,1),'.b');
                
                
                microScale = 0.65;
                LinRefLength = LenRefPixels*microScale; % micro meter
                % advection time (beginning is always closest to bifurcation center)
                AdvectionTimeHealthy = LinRefLength/HRBCVelMothBeg; % um/s and um
                AdvectionTimeRigid = LinRefLength/RRBCVelMothBeg; % um/s and um
                % calculate
                HRBCResidenceFrames = [bifFrameNum(1).clu.frames];
                RRBCResidenceFrames = [bifFrameNum(2).clu.frames];
                HRBCApexDistancesPixel = [distApex(1).clu.minApexDist];
                RRBCApexDistancesPixel = [distApex(2).clu.minApexDist];
                HRBCBifTrajLen = [bifTrajLen(1).clu.distance];
                RRBCBifTrajLen = [bifTrajLen(2).clu.distance];
                % conditioning
                condHRBC = HRBCResidenceFrames>2;
                condRRBC = RRBCResidenceFrames>2;
                HRBCResidenceFrames = HRBCResidenceFrames(condHRBC);
                RRBCResidenceFrames = RRBCResidenceFrames(condRRBC);
                HRBCApexDistancesPixel = HRBCApexDistancesPixel(condHRBC);
                RRBCApexDistancesPixel = RRBCApexDistancesPixel(condRRBC);
                HRBCBifTrajLen = HRBCBifTrajLen(condHRBC);
                RRBCBifTrajLen = RRBCBifTrajLen(condRRBC);
                % lingering
                ApexDistanceHealthy = HRBCApexDistancesPixel*microScale;
                ApexDistanceRigid = RRBCApexDistancesPixel*microScale;
                ResidenceTimeHealthy = HRBCResidenceFrames/frameRates(1);
                ResidenceTimeRigid = RRBCResidenceFrames/frameRates(2);
                LingeringHealthy = ResidenceTimeHealthy./AdvectionTimeHealthy;
                LingeringRigid = ResidenceTimeRigid./AdvectionTimeHealthy;
                AdvectionLengthHealthy = HRBCBifTrajLen*microScale;
                AdvectionLengthRigid = RRBCBifTrajLen*microScale;

                % saving
                bifurcations(bifIdx).AdvectionTimeHealthy = AdvectionTimeHealthy;
                bifurcations(bifIdx).AdvectionTimeRigid = AdvectionTimeRigid;
                bifurcations(bifIdx).HRBCVelMoth = HRBCVelMoth;
                bifurcations(bifIdx).RRBCVelMoth = RRBCVelMoth;
                bifurcations(bifIdx).HRBCVelMothEnd = HRBCVelMothEnd;
                bifurcations(bifIdx).RRBCVelMothEnd = RRBCVelMothEnd;
                bifurcations(bifIdx).HRBCVelMothBeg = HRBCVelMothBeg;
                bifurcations(bifIdx).RRBCVelMothBeg = RRBCVelMothBeg;
                bifurcations(bifIdx).LinRefLenght = LinRefLength;
                bifurcations(bifIdx).ResidenceTimeHealthy = ResidenceTimeHealthy;
                bifurcations(bifIdx).ResidenceTimeRigid = ResidenceTimeRigid;
                bifurcations(bifIdx).LingeringHealthy = LingeringHealthy;
                bifurcations(bifIdx).LingeringRigid = LingeringRigid;
                bifurcations(bifIdx).ApexDistanceHealthy = ApexDistanceHealthy;
                bifurcations(bifIdx).ApexDistanceRigid = ApexDistanceRigid;
                bifurcations(bifIdx).AdvectionLengthHealthy = AdvectionLengthHealthy;
                bifurcations(bifIdx).AdvectionLengthRigid = AdvectionLengthRigid;
                %% calculate bifurcation velocity
                % healthy
                cond = ~isnan(AdvectionLengthHealthy+ResidenceTimeHealthy) & ~isinf(AdvectionLengthHealthy+ResidenceTimeHealthy);
                AdvectionSpeed = AdvectionLengthHealthy(cond)./ResidenceTimeHealthy(cond);
                allHRBCVelBif = [allHRBCVelBif,real(mean(AdvectionSpeed(~isinf(AdvectionSpeed)&~isnan(AdvectionSpeed))))];
                allHRBCNumBif = [allHRBCNumBif,length(AdvectionSpeed(~isinf(AdvectionSpeed)&~isnan(AdvectionSpeed)))];
                ErrallHRBCVelBif = [ErrallHRBCVelBif,std(AdvectionSpeed(~isinf(AdvectionSpeed)&~isnan(AdvectionSpeed)))];%/sqrt(length(AdvectionSpeed(~isinf(AdvectionSpeed)&~isnan(AdvectionSpeed))))];
                % rigid
                cond = ~isnan(AdvectionLengthRigid+ResidenceTimeRigid) & ~isinf(AdvectionLengthRigid+ResidenceTimeRigid);
                AdvectionSpeed = AdvectionLengthRigid(cond)./ResidenceTimeRigid(cond);
                allRRBCVelBif = [allRRBCVelBif,real(mean(AdvectionSpeed(~isinf(AdvectionSpeed)&~isnan(AdvectionSpeed))))];
                allRRBCNumBif = [allRRBCNumBif,length(AdvectionSpeed(~isinf(AdvectionSpeed)&~isnan(AdvectionSpeed)))];
                ErrallRRBCVelBif = [ErrallRRBCVelBif,std(AdvectionSpeed(~isinf(AdvectionSpeed)&~isnan(AdvectionSpeed)))];%/sqrt(length(AdvectionSpeed(~isinf(AdvectionSpeed)&~isnan(AdvectionSpeed))))];
                %% lingering
                % healthy
                allHRBCLin = [allHRBCLin,mean(LingeringHealthy(~isnan(LingeringHealthy)&~isinf(LingeringHealthy)))];
                ErrallHRBCLin = [ErrallHRBCLin,std(LingeringHealthy(~isnan(LingeringHealthy)&~isinf(LingeringHealthy)))/sqrt(length(LingeringHealthy(~isnan(LingeringHealthy)&~isinf(LingeringHealthy))))];
                % rigid
                allRRBCLin = [allRRBCLin,mean(LingeringRigid(~isnan(LingeringRigid)&~isinf(LingeringRigid)))];
                ErrallRRBCLin = [ErrallRRBCLin,std(LingeringRigid(~isnan(LingeringRigid)&~isinf(LingeringRigid)))/sqrt(length(LingeringRigid(~isnan(LingeringRigid)&~isinf(LingeringRigid))))];
                %% save others
                % minimum number of cells in vessel
                minNumCells = 5;
                % mother velocities
                if isempty(HRBCVelMoth) || HRBCNumMoth < minNumCells
                    HRBCVelMoth = nan;
                end
                if isempty(RRBCVelMoth) || RRBCNumMoth < minNumCells
                    RRBCVelMoth = nan;
                end

                % daughter velocities
                HRBCVelDaugh = [cVes(selDaugh).HRBCVel];
                RRBCVelDaugh = [cVes(selDaugh).RRBCVel];

                %daughter numbers
                HRBCNumDaugh = [cVes(selDaugh).HRBCNum];
                RRBCNumDaugh = [cVes(selDaugh).RRBCNum];
                % daughter velocities
                HRBCVelDaugh(HRBCNumDaugh<minNumCells) = nan;
                RRBCVelDaugh(RRBCNumDaugh<minNumCells) = nan;


                allHRBCVelMoth = [allHRBCVelMoth,HRBCVelMoth];
                allRRBCVelMoth = [allRRBCVelMoth,RRBCVelMoth];

                allHRBCVelMothEnd = [allHRBCVelMothEnd,HRBCVelMothEnd];
                allRRBCVelMothEnd = [allRRBCVelMothEnd,RRBCVelMothEnd];

                allHRBCVelMothBeg = [allHRBCVelMothBeg,HRBCVelMothBeg];
                allRRBCVelMothBeg = [allRRBCVelMothBeg,RRBCVelMothBeg];

                ErrallHRBCVelMoth = [ErrallHRBCVelMoth, HRBCVelMothErr];
                ErrallRRBCVelMoth = [ErrallRRBCVelMoth, RRBCVelMothErr];

                % mother number 
                allHRBCNumMoth = [allHRBCNumMoth,HRBCNumMoth];
                allRRBCNumMoth = [allRRBCNumMoth,RRBCNumMoth];

                allHRBCVelDaugh = [allHRBCVelDaugh,HRBCVelDaugh];
                allRRBCVelDaugh = [allRRBCVelDaugh,RRBCVelDaugh];

                allHRBCVelDaughEnd = [allHRBCVelDaughEnd,cVes(selDaugh).HRBCVelEnd];
                allRRBCVelDaughEnd = [allRRBCVelDaughEnd,cVes(selDaugh).RRBCVelEnd];

                allHRBCVelDaughBeg = [allHRBCVelDaughBeg,cVes(selDaugh).HRBCVelBeg];
                allRRBCVelDaughBeg = [allRRBCVelDaughBeg,cVes(selDaugh).RRBCVelBeg];

                ErrHRBCVelDaugh = [];
                ErrRRBCVelDaugh = [];
                selDaughIdc = find(selDaugh);
                for daughIdx = 1:length(selDaughIdc)
                    errHRBC = std([cVes(selDaugh(daughIdx)).HRBCVelArray])*bifurcations(1).velScaleHRBC;
                    if isempty(errHRBC)
                        errHRBC = nan;
                    end
                    errRRBC = std([cVes(selDaugh(daughIdx)).RRBCVelArray])*bifurcations(1).velScaleRRBC;
                    if isempty(errRRBC)
                        errRRBC = nan;
                    end
                    ErrHRBCVelDaugh = [ErrHRBCVelDaugh,errHRBC];
                    ErrRRBCVelDaugh = [ErrRRBCVelDaugh,errRRBC];
                end
                ErrallHRBCVelDaugh = [ErrallHRBCVelDaugh,ErrHRBCVelDaugh];
                ErrallRRBCVelDaugh = [ErrallRRBCVelDaugh,ErrRRBCVelDaugh];
                
                % daughter number 
                allHRBCNumDaugh = [allHRBCNumDaugh,cVes(selDaugh).HRBCNum];
                allRRBCNumDaugh = [allRRBCNumDaugh,cVes(selDaugh).RRBCNum];

                % geometry
                MothLen = max(max(cVes(inVesIdx).length));
                allMothLens = [allMothLens,MothLen];
                %%
                if flag_plotting % optional plotting
                    sum(LingeringHealthy > 1)/length(LingeringHealthy)
                    sum(LingeringRigid > 1)/length(LingeringRigid)
                    if (length(LingeringHealthy) > 5 && length(LingeringRigid) > 5)
                        if (sum(isnan(mean(LingeringHealthy)+mean(LingeringRigid)))==0)
                            %figure
                            try % choose which figures to show
                                if 0 % figure of sort (CDF)
                                    figure
                                    hold on
                                    plot((1:length(LingeringHealthy))/length(LingeringHealthy),sort(LingeringHealthy),'.-r')

                                    %ksdensity(LingeringHealthy)

                                    plot((1:length(LingeringRigid))/length(LingeringRigid),sort(LingeringRigid),'.-b')
                                    %ksdensity(LingeringRigid)
                                    hold off
                                    legend(sprintf('H: %.4f',mean(LingeringHealthy)),sprintf('R: %.4f',mean(LingeringRigid)))
                                    pause(0.01);
                                    title(filePath)
                                end
                                if 0 % apex distance
                                    figure
                                    hold on
                                    plot(ApexDistanceHealthy,LingeringHealthy,'.r');
                                    plot(ApexDistanceRigid,LingeringRigid,'.b');
                                    legend(sprintf('H: %.4f',mean(LingeringHealthy)),sprintf('R: %.4f',mean(LingeringRigid)))
                                    title(filePath)
                                    hold off
                                end
                                if 0 % advection length
                                    figure
                                    hold on
                                    plot(AdvectionLengthHealthy,LingeringHealthy,'.r');
                                    plot(AdvectionLengthRigid,LingeringRigid,'.b');
                                    legend(sprintf('H: %.4f',mean(LingeringHealthy)),sprintf('R: %.4f',mean(LingeringRigid)))
                                    title(filePath)
                                    hold off
                                end
                                pause(0.01)
                            end
                        end
                    end
                end
            catch ME
                disp(ME.message);
            end
        end % end of bifurcation loop
        hold off;
        saveNameMotherDet = filePath;
        saveNameMotherDet = strrep(saveNameMotherDet,'\','-');
        saveNameMotherDet = strrep(saveNameMotherDet,':','-');
        saveNameMotherDet = saveNameMotherDet(1:end-4);
        print(sprintf('img-%s-BIF%d-mother_detection.png',saveNameMotherDet,bifIdx),'-dpng','-r300');
        close all;

        save([rootDir,'\bifLin.mat'],"bifurcations");
    end
    % save velocity data
    velData.allHRBCVelMoth = allHRBCVelMoth;
    velData.allRRBCVelMoth = allRRBCVelMoth;

    velData.allHRBCVelDaugh = allHRBCVelDaugh;
    velData.allRRBCVelDaugh = allRRBCVelDaugh;

    velData.allHRBCVelBif = allHRBCVelBif;
    velData.allRRBCVelBif = allRRBCVelBif;

    velData.allHRBCLin = allHRBCLin;
    velData.allRRBCLin = allRRBCLin;

    save('velData.mat','velData');
    save('LingeringAnalysisResultsData.mat');
end
%% continue with statistical meta-analysis
load('LingeringAnalysisResultsData.mat')
load("velData.mat");
clc, close all;
%% renormalize
renormalize = false;
if renormalize
    % mother
    allRRBCVelMothBeg = allRRBCVelMothBeg./allRRBCVelMoth;
    allRRBCVelMothEnd = allRRBCVelMothEnd./allRRBCVelMoth;
    allHRBCVelMothBeg = allHRBCVelMothBeg./allHRBCVelMoth;
    allHRBCVelMothEnd = allHRBCVelMothEnd./allHRBCVelMoth;
    % bifurcation
    allRRBCVelBif = allRRBCVelBif./allRRBCVelMoth;
    allHRBCVelBif = allHRBCVelBif./allHRBCVelMoth;
    % daughter
    repRRBCData = repmat(allRRBCVelMoth,2,1);repRRBCData = repRRBCData(:)';
    repHRBCData = repmat(allHRBCVelMoth,2,1);repHRBCData = repHRBCData(:)';
    allRRBCVelDaughBeg = allRRBCVelDaughBeg./repRRBCData;
    allRRBCVelDaughEnd = allRRBCVelDaughEnd./repRRBCData;
    allHRBCVelDaughBeg = allHRBCVelDaughBeg./repHRBCData;
    allHRBCVelDaughEnd = allHRBCVelDaughEnd./repHRBCData;
    allHRBCVelDaugh = allHRBCVelDaugh./repHRBCData;
    allRRBCVelDaugh = allRRBCVelDaugh./repRRBCData;
end
%%
% color definition
linewidth = 1;
lambda = linspace(0,1,100);
R = double(lambda > 0.5).*(lambda-0.5)/0.5;
G = 0*R;
B = flip(R);
smR = fit(linspace(-1,1,100)',R','smoothingspline','smoothingparam',1);
smG = fit(linspace(-1,1,100)',G','smoothingspline','smoothingparam',1);
smB = fit(linspace(-1,1,100)',B','smoothingspline','smoothingparam',1);
%% mother speeds
close all
xData = allRRBCVelMoth; % daughter beginning is closer to bifurcation
yData = allHRBCVelMoth; % daughter end is further away from bifurcation

xErr = ErrallRRBCVelMoth;
yErr = ErrallHRBCVelMoth;

xErr(isnan(xErr)) = median(xErr,'omitnan');
yErr(isnan(yErr)) = median(yErr,'omitnan');

cond = ~isinf(xData)&~isinf(yData)&~isnan(xData)&~isnan(yData);

xData = xData(cond);
yData = yData(cond); if isempty(xData) || isempty(yData), error('x- or y-data empty'), end, if length(xData)<3 || length(yData)<3, error('not enough data points'), end
xErr = xErr(cond);
yErr = yErr(cond);


axLims = [min([xData,yData]) max([xData,yData])];
axLims = [0 1000];

fig = figure();
hold on
x = axLims;


% plot linear model
mdl = fitlm(xData,yData,'intercept',true,'RobustOpts','on');
lm_slope = mdl.Coefficients.Estimate(2);
lm_intcp = mdl.Coefficients.Estimate(1);
lm_var = mdl.Coefficients.SE(2);
lm_intvar = mdl.Coefficients.SE(1);
xLine = linspace(-1000,1000,2);
lineMid = (lm_slope)*xLine + lm_intcp;
lineTop = (lm_slope+lm_var)*xLine + lm_intcp+lm_intvar;
lineBot = (lm_slope-lm_var)*xLine + lm_intcp-lm_intvar;
fill([xLine,flip(xLine)],[lineTop, fliplr(lineBot)],[0.85 0.85 0.85],'EdgeColor','none')
lineTop = (lm_slope+lm_var)*xLine + lm_intcp-lm_intvar;
lineBot = (lm_slope-lm_var)*xLine + lm_intcp+lm_intvar;
fill([xLine,flip(xLine)],[lineTop, fliplr(lineBot)],[0.85 0.85 0.85],'EdgeColor','none')
plot(x,x,'--','LineWidth',1.4,'color',[0.2 0.2 0.2])

linewidth = 1;
for pntIdx = 1:length(xData)
    xVal = xData(pntIdx);
    yVal = yData(pntIdx);
    scale = abs(xData-yData);dist = -(xVal-yVal)/sqrt(2)/100;
    dist(isnan(dist)) = 0; pointColor = [smR(dist),smG(dist),smB(dist)];
    pointColor(pointColor>1)=1;pointColor(isnan(pointColor))=0;
    pointColor(pointColor<0)=0;
    errorbar(xData(pntIdx),yData(pntIdx),yErr(pntIdx),yErr(pntIdx),xErr(pntIdx),xErr(pntIdx),'s','color',pointColor,'MarkerSize',6,'MarkerFaceColor','none','MarkerEdgeColor','none','LineWidth',linewidth)
    s = scatter(xData(pntIdx),yData(pntIdx),20,'color',pointColor,'Marker','o','MarkerFaceColor',pointColor,'MarkerFaceAlpha',0.2,'MarkerEdgeColor',pointColor);
end

xlh = xlabel(['\fontname{cmss10}v_{\fontname{mwa_cmmi10}R}^{\fontname{mwa_cmmi10}M}\fontname{cmss10} [' char(181) 'm/s]']);%,'Interpreter','tex')
ylh = ylabel(['\fontname{cmss10}v_{\fontname{mwa_cmmi10}H}^{\fontname{mwa_cmmi10}M}\fontname{cmss10} [' char(181) 'm/s]']);%,'Interpreter','tex')

signRankTestLine;

ylh.Position(1) = ylh.Position(1) + 180;
xlh.Position(2) = xlh.Position(2) + 100;


xticks([0,500,1000])
yticks([0,500,1000])
xticklabels({0,'',1000})
yticklabels({0,'',1000})
ax.TickLength = [0.03 0.03];
title('$v$(moth)','interpreter','latex')
exportgraphics(fig,'OUT-mother_vels.png','Resolution',600);
saveas(fig,'OUT-mother_vels.svg');

%% daughter vels

xData = allRRBCVelDaugh; % daughter beginning is closer to bifurcation
yData = allHRBCVelDaugh; % daughter end is further away from bifurcation

xErr = ErrallRRBCVelDaugh;
yErr = ErrallHRBCVelDaugh;

xErr(isnan(xErr)) = median(xErr,'omitnan');
yErr(isnan(yErr)) = median(yErr,'omitnan');

cond = ~isinf(xData)&~isinf(yData)&~isnan(xData)&~isnan(yData);

xData = xData(cond);
yData = yData(cond); if isempty(xData) || isempty(yData), error('x- or y-data empty'), end, if length(xData)<3 || length(yData)<3, error('not enough data points'), end
xErr = xErr(cond);
yErr = yErr(cond);

axLims = [0 max([xData,yData])];
axLims = [0 1000];

fig = figure();
hold on
x = axLims;


% plot linear model
mdl = fitlm(xData,yData,'intercept',true,'RobustOpts','on');
lm_slope = mdl.Coefficients.Estimate(2);
lm_intcp = mdl.Coefficients.Estimate(1);
lm_var = mdl.Coefficients.SE(2);
lm_intvar = mdl.Coefficients.SE(1);
xLine = linspace(-1000,1000,2);
lineMid = (lm_slope)*xLine + lm_intcp;
lineTop = (lm_slope+lm_var)*xLine + lm_intcp+lm_intvar;
lineBot = (lm_slope-lm_var)*xLine + lm_intcp-lm_intvar;
fill([xLine,flip(xLine)],[lineTop, fliplr(lineBot)],[0.85 0.85 0.85],'EdgeColor','none')
lineTop = (lm_slope+lm_var)*xLine + lm_intcp-lm_intvar;
lineBot = (lm_slope-lm_var)*xLine + lm_intcp+lm_intvar;
fill([xLine,flip(xLine)],[lineTop, fliplr(lineBot)],[0.85 0.85 0.85],'EdgeColor','none')

plot(x,x,'--','LineWidth',1.4,'color',[0.2 0.2 0.2])

linewidth = 1;
for pntIdx = 1:length(xData)
    xVal = xData(pntIdx);
    yVal = yData(pntIdx);
    scale = abs(xData-yData);dist = -(xVal-yVal)/sqrt(2)/100;
    dist(isnan(dist)) = 0; pointColor = [smR(dist),smG(dist),smB(dist)];
    pointColor(pointColor>1)=1;pointColor(isnan(pointColor))=0;
    pointColor(pointColor<0)=0;
    errorbar(xData(pntIdx),yData(pntIdx),yErr(pntIdx),yErr(pntIdx),xErr(pntIdx),xErr(pntIdx),'s','color',pointColor,'MarkerSize',6,'MarkerFaceColor','none','MarkerEdgeColor','none','LineWidth',linewidth)
    s = scatter(xData(pntIdx),yData(pntIdx),20,'color',pointColor,'Marker','o','MarkerFaceColor',pointColor,'MarkerFaceAlpha',0.2,'MarkerEdgeColor',pointColor);
end

xlh = xlabel(['\fontname{cmss10}v_{\fontname{mwa_cmmi10}R}^{\fontname{mwa_cmmi10}D}\fontname{cmss10} [' char(181) 'm/s]']);%,'Interpreter','tex')
ylh = ylabel(['\fontname{cmss10}v_{\fontname{mwa_cmmi10}H}^{\fontname{mwa_cmmi10}D}\fontname{cmss10} [' char(181) 'm/s]']);%,'Interpreter','tex')

signRankTestLine;

ylh.Position(1) = ylh.Position(1) + 180;
xlh.Position(2) = xlh.Position(2) + 100;

xticks([0,500,1000])
yticks([0,500,1000])
xticklabels({0,'',1000})
yticklabels({0,'',1000})
ax.TickLength = [0.03 0.03];

title('$v$(daugh)','interpreter','latex')

exportgraphics(fig,'OUT-daughter_vels.png','Resolution',600);
saveas(fig,'OUT-daughter_vels.svg');
%% bifurcation vels

xData = allRRBCVelBif; % daughter beginning is closer to bifurcation
yData = allHRBCVelBif; % daughter end is further away from bifurcation

xErr = ErrallRRBCVelBif;
yErr = ErrallHRBCVelBif;

xErr(isnan(xErr)) = median(xErr,'omitnan');
yErr(isnan(yErr)) = median(yErr,'omitnan');

cond = ~isinf(xData)&~isinf(yData)&~isnan(xData)&~isnan(yData);

xData = xData(cond);
yData = yData(cond); if isempty(xData) || isempty(yData), error('x- or y-data empty'), end, if length(xData)<3 || length(yData)<3, error('not enough data points'), end
xErr = xErr(cond);
yErr = yErr(cond);


axLims = [0 max([xData,yData])];
axLims = [0 1000];

fig = figure();
hold on
x = axLims;

% plot linear model
mdl = fitlm(xData,yData,'intercept',true,'RobustOpts','on');
lm_slope = mdl.Coefficients.Estimate(2);
lm_intcp = mdl.Coefficients.Estimate(1);
lm_var = mdl.Coefficients.SE(2);
lm_intvar = mdl.Coefficients.SE(1);
xLine = linspace(-1000,1000,2);
lineMid = (lm_slope)*xLine + lm_intcp;
lineTop = (lm_slope+lm_var)*xLine + lm_intcp+lm_intvar;
lineBot = (lm_slope-lm_var)*xLine + lm_intcp-lm_intvar;
fill([xLine,flip(xLine)],[lineTop, fliplr(lineBot)],[0.85 0.85 0.85],'EdgeColor','none')
lineTop = (lm_slope+lm_var)*xLine + lm_intcp-lm_intvar;
lineBot = (lm_slope-lm_var)*xLine + lm_intcp+lm_intvar;
fill([xLine,flip(xLine)],[lineTop, fliplr(lineBot)],[0.85 0.85 0.85],'EdgeColor','none')


plot(x,x,'--','LineWidth',1.4,'color',[0.2 0.2 0.2])

linewidth = 1;
for pntIdx = 1:length(xData)
    xVal = xData(pntIdx);
    yVal = yData(pntIdx);
    scale = abs(xData-yData);dist = -(xVal-yVal)/sqrt(2)/100;
    dist = real(dist);
    dist(isnan(dist)) = 0; pointColor = [smR(dist),smG(dist),smB(dist)];
    pointColor(pointColor>1)=1;pointColor(isnan(pointColor))=0;
    pointColor(pointColor<0)=0;
    errorbar(xData(pntIdx),yData(pntIdx),yErr(pntIdx),yErr(pntIdx),xErr(pntIdx),xErr(pntIdx),'s','color',pointColor,'MarkerSize',6,'MarkerFaceColor','none','MarkerEdgeColor','none','LineWidth',linewidth)
    s = scatter(xData(pntIdx),yData(pntIdx),20,'color',pointColor,'Marker','o','MarkerFaceColor',pointColor,'MarkerFaceAlpha',0.2,'MarkerEdgeColor',pointColor);
end

xlh = xlabel(['\fontname{cmss10}v_{\fontname{mwa_cmmi10}R}^{\fontname{mwa_cmmi10}B}\fontname{cmss10} [' char(181) 'm/s]']);%,'Interpreter','tex')
ylh = ylabel(['\fontname{cmss10}v_{\fontname{mwa_cmmi10}H}^{\fontname{mwa_cmmi10}B}\fontname{cmss10} [' char(181) 'm/s]']);%,'Interpreter','tex')

signRankTestLine;

title('$v$(bif)','interpreter','latex')

ylh.Position(1) = ylh.Position(1) + 180;
xlh.Position(2) = xlh.Position(2) + 100;

xticks([0,500,1000])
yticks([0,500,1000])
xticklabels({0,'',1000})
yticklabels({0,'',1000})
ax.TickLength = [0.03 0.03];

exportgraphics(fig,'OUT-bifurcation_vels.png','Resolution',600);
saveas(fig,'OUT-bifurcation_vels.svg');
%% difference mother and bifurcation

xData = (allRRBCVelBif-allRRBCVelMoth); % daughter beginning is closer to bifurcation
yData = (allHRBCVelBif-allHRBCVelMoth); % daughter end is further away from bifurcation

xErr = sqrt(ErrallRRBCVelBif.^2+ErrallRRBCVelMoth.^2);
yErr = sqrt(ErrallHRBCVelBif.^2+ErrallHRBCVelMoth.^2);

xErr(isnan(xErr)) = median(xErr,'omitnan');
yErr(isnan(yErr)) = median(yErr,'omitnan');

cond = ~isnan(xData)&~isnan(yData)&~isinf(xData)&~isinf(yData);

xData = xData(cond);
yData = yData(cond); if isempty(xData) || isempty(yData), error('x- or y-data empty'), end, if length(xData)<3 || length(yData)<3, error('not enough data points'), end
xErr = xErr(cond);
yErr = yErr(cond);


axLims = [min([xData,yData]) max([xData,yData])];
axLims = [-750,750];

fig = figure();
hold on
x = axLims;

% plot linear model
weights = 1./(xErr+yErr); weights = weights/mean(weights);  
mdl = fitlm(xData,yData,'intercept',true,'RobustOpts','on','weights',weights);
lm_slope = mdl.Coefficients.Estimate(2);
lm_intcp = mdl.Coefficients.Estimate(1);
lm_var = mdl.Coefficients.SE(2);
lm_intvar = mdl.Coefficients.SE(1);
xLine = linspace(-1000,1000,2);
lineMid = (lm_slope)*xLine + lm_intcp;
lineTop = (lm_slope+lm_var)*xLine + lm_intcp+lm_intvar;
lineBot = (lm_slope-lm_var)*xLine + lm_intcp-lm_intvar;
fill([xLine,flip(xLine)],[lineTop, fliplr(lineBot)],[0.85 0.85 0.85],'EdgeColor','none')
lineTop = (lm_slope+lm_var)*xLine + lm_intcp-lm_intvar;
lineBot = (lm_slope-lm_var)*xLine + lm_intcp+lm_intvar;
fill([xLine,flip(xLine)],[lineTop, fliplr(lineBot)],[0.85 0.85 0.85],'EdgeColor','none')


plot(x,x,'--','LineWidth',1.4,'color',[0.2 0.2 0.2])
plot(x,0*x,'-','LineWidth',1.4,'color',[0.2 0.2 0.2])
plot(0*x,x,'-','LineWidth',1.4,'color',[0.2 0.2 0.2])

linewidth = 1;
for pntIdx = 1:length(xData)
    xVal = xData(pntIdx);
    yVal = yData(pntIdx);
    scale = abs(xData-yData);dist = -(xVal-yVal)/sqrt(2)/100;
    dist = real(dist);
    dist(isnan(dist)) = 0; pointColor = [smR(dist),smG(dist),smB(dist)];
    pointColor(pointColor>1)=1;pointColor(isnan(pointColor))=0;
    pointColor(pointColor<0)=0;
    errorbar(xData(pntIdx),yData(pntIdx),yErr(pntIdx),yErr(pntIdx),xErr(pntIdx),xErr(pntIdx),'s','color',pointColor,'MarkerSize',6,'MarkerFaceColor','none','MarkerEdgeColor','none','LineWidth',linewidth)
    s = scatter(xData(pntIdx),yData(pntIdx),20,'color',pointColor,'Marker','o','MarkerFaceColor',pointColor,'MarkerFaceAlpha',0.2,'MarkerEdgeColor',pointColor);
end

xlh = xlabel(['\Delta\fontname{cmss10}v_{\fontname{mwa_cmmi10}R}\fontname{cmss10} [' char(181) 'm/s]']);%,'Interpreter','tex')
ylh = ylabel(['\Delta\fontname{cmss10}v_{\fontname{mwa_cmmi10}H}\fontname{cmss10} [' char(181) 'm/s]']);%,'Interpreter','tex')

signRankTestLine;

ylh.Position(1) = ylh.Position(1) + 180;
xlh.Position(2) = xlh.Position(2) + 100;


xticks([-750 0 750])
yticks([-750 0 750])
xticklabels({-750,'',750})
yticklabels({-750,'',750})
ax.TickLength = [0.03 0.03];

title('$\Delta v$(bif-moth)','interpreter','latex')
exportgraphics(fig,'OUT-diff_moth_bif.png','Resolution',600);
saveas(fig,'OUT-diff_moth_bif.svg');
%% difference bifurcation and daughter

repRRBCData = repmat(allRRBCVelBif,2,1);
repHRBCData = repmat(allHRBCVelBif,2,1);

xData = allRRBCVelDaugh-repRRBCData(:)'; % daughter beginning is closer to bifurcation
yData = allHRBCVelDaugh-repHRBCData(:)'; % daughter end is further away from bifurcation


repErrRRBCData = repmat(ErrallRRBCVelBif',2,1)';
repErrHRBCData = repmat(ErrallHRBCVelBif',2,1)';

xErr = sqrt(repErrRRBCData.^2+ErrallRRBCVelDaugh.^2);
yErr = sqrt(repErrHRBCData.^2+ErrallHRBCVelDaugh.^2);

xErr(isnan(xErr)) = median(xErr,'omitnan');
yErr(isnan(yErr)) = median(yErr,'omitnan');

cond = ~isnan(xData)&~isnan(yData)&~isinf(xData)&~isinf(yData);

xData = xData(cond);
yData = yData(cond); if isempty(xData) || isempty(yData), error('x- or y-data empty'), end, if length(xData)<3 || length(yData)<3, error('not enough data points'), end
xErr = xErr(cond);
yErr = yErr(cond);

axLims = [min([xData,yData]) max([xData,yData])];
axLims = [-750,750];


fig = figure();
hold on
x = axLims;

% plot linear model
weights = 1./(xErr+yErr); weights = weights/mean(weights);  
mdl = fitlm(xData,yData,'intercept',true,'RobustOpts','on','weights',weights);
lm_slope = mdl.Coefficients.Estimate(2);
lm_intcp = mdl.Coefficients.Estimate(1);
lm_var = mdl.Coefficients.SE(2);
lm_intvar = mdl.Coefficients.SE(1);
xLine = linspace(-1000,1000,2);
lineMid = (lm_slope)*xLine + lm_intcp;
lineTop = (lm_slope+lm_var)*xLine + lm_intcp+lm_intvar;
lineBot = (lm_slope-lm_var)*xLine + lm_intcp-lm_intvar;
fill([xLine,flip(xLine)],[lineTop, fliplr(lineBot)],[0.85 0.85 0.85],'EdgeColor','none')
lineTop = (lm_slope+lm_var)*xLine + lm_intcp-lm_intvar;
lineBot = (lm_slope-lm_var)*xLine + lm_intcp+lm_intvar;
fill([xLine,flip(xLine)],[lineTop, fliplr(lineBot)],[0.85 0.85 0.85],'EdgeColor','none')


plot(x,x,'--','LineWidth',1.4,'color',[0.2 0.2 0.2])
plot(x,0*x,'-','LineWidth',1.4,'color',[0.2 0.2 0.2])
plot(0*x,x,'-','LineWidth',1.4,'color',[0.2 0.2 0.2])

linewidth = 1;
for pntIdx = 1:length(xData)
    xVal = xData(pntIdx);
    yVal = yData(pntIdx);
    scale = abs(xData-yData);dist = -(xVal-yVal)/sqrt(2)/100;%)-(xVal-yVal)/sqrt(2)/mean(abs(xData-yData));%)-10*(xVal-yVal)/sqrt(2)/mean(abs(xData-yData));%)-100*(xVal-yVal)/sqrt(2)/mean(abs(xData-xData));%)-100*(xVal-yVal)/sqrt(2)/mean(abs(xVal-yVal));%)-10*(xVal-yVal)/sqrt(2)/mean(abs(xVal-yVal));%)-(xVal-yVal)/sqrt(2)/mean(abs(xVal-yVal));%)-(xVal-yVal)/sqrt(2)/(mean(abs(xVal-yVal));%)mean(abs(xVal-yVal));%mean(abs(xVal-yVal));%dist/50*250;
    dist = real(dist);
    dist(isnan(dist)) = 0; pointColor = [smR(dist),smG(dist),smB(dist)];
    pointColor(pointColor>1)=1;pointColor(isnan(pointColor))=0;
    pointColor(pointColor<0)=0;
    errorbar(xData(pntIdx),yData(pntIdx),yErr(pntIdx),yErr(pntIdx),xErr(pntIdx),xErr(pntIdx),'s','color',pointColor,'MarkerSize',6,'MarkerFaceColor','none','MarkerEdgeColor','none','LineWidth',linewidth)
    s = scatter(xData(pntIdx),yData(pntIdx),20,'color',pointColor,'Marker','o','MarkerFaceColor',pointColor,'MarkerFaceAlpha',0.2,'MarkerEdgeColor',pointColor);
end
xlh = xlabel(['\Delta\fontname{cmss10}v_{\fontname{mwa_cmmi10}R}\fontname{cmss10} [' char(181) 'm/s]']);%,'Interpreter','tex')
ylh = ylabel(['\Delta\fontname{cmss10}v_{\fontname{mwa_cmmi10}H}\fontname{cmss10} [' char(181) 'm/s]']);%,'Interpreter','tex')

signRankTestLine;

ylh.Position(1) = ylh.Position(1) + 180;
xlh.Position(2) = xlh.Position(2) + 100;


xticks([-750 0 750])
yticks([-750 0 750])
xticklabels({-750,'',750})
yticklabels({-750,'',750})

ax.TickLength = [0.03 0.03];
title('$\Delta v$(daugh-bif)','interpreter','latex')

exportgraphics(fig,'OUT-diff_bif_daugh.png','Resolution',600);
saveas(fig,'OUT-diff_bif_daugh.svg');
%% difference daughter and mother

repRRBCData = repmat(allRRBCVelMoth,2,1);
repHRBCData = repmat(allHRBCVelMoth,2,1);

xData = -repRRBCData(:)'+allRRBCVelDaugh; % daughter beginning is closer to bifurcation
yData = -repHRBCData(:)'+allHRBCVelDaugh; % daughter end is further away from bifurcation

ErrallRRBCVelDaugh(isnan(ErrallRRBCVelDaugh)) = median(ErrallRRBCVelDaugh,'omitnan');
ErrallHRBCVelDaugh(isnan(ErrallHRBCVelDaugh)) = median(ErrallHRBCVelDaugh,'omitnan');

repErrRRBCData = repmat(ErrallRRBCVelMoth',2,1)';
repErrHRBCData = repmat(ErrallHRBCVelMoth',2,1)';
repErrRRBCData(isnan(repErrRRBCData)) = median(repErrRRBCData,'omitnan');
repErrHRBCData(isnan(repErrHRBCData)) = median(repErrHRBCData,'omitnan');
xErr = sqrt(repErrRRBCData.^2+ErrallRRBCVelDaugh.^2);
yErr = sqrt(repErrHRBCData.^2+ErrallHRBCVelDaugh.^2);

cond = ~isnan(xData)&~isnan(yData)&~isinf(xData)&~isinf(yData);

xErr = xErr(cond);
yErr = yErr(cond);

xData = xData(cond);
yData = yData(cond); if isempty(xData) || isempty(yData), error('x- or y-data empty'), end, if length(xData)<3 || length(yData)<3, error('not enough data points'), end

axLims = [min([xData,yData]) max([xData,yData])];
axLims = [-750,750];


fig = figure();
hold on
x = axLims;


% plot linear model
weights = 1./(xErr+yErr); weights = weights/mean(weights);  
mdl = fitlm(xData,yData,'intercept',true,'RobustOpts','on','weights',weights);
lm_slope = mdl.Coefficients.Estimate(2);
lm_intcp = mdl.Coefficients.Estimate(1);
lm_var = mdl.Coefficients.SE(2);
lm_intvar = mdl.Coefficients.SE(1);
xLine = linspace(-1000,1000,2);
lineMid = (lm_slope)*xLine + lm_intcp;
lineTop = (lm_slope+lm_var)*xLine + lm_intcp+lm_intvar;
lineBot = (lm_slope-lm_var)*xLine + lm_intcp-lm_intvar;
fill([xLine,flip(xLine)],[lineTop, fliplr(lineBot)],[0.85 0.85 0.85],'EdgeColor','none')
lineTop = (lm_slope+lm_var)*xLine + lm_intcp-lm_intvar;
lineBot = (lm_slope-lm_var)*xLine + lm_intcp+lm_intvar;
fill([xLine,flip(xLine)],[lineTop, fliplr(lineBot)],[0.85 0.85 0.85],'EdgeColor','none')


plot(x,x,'--','LineWidth',1.4,'color',[0.2 0.2 0.2])
plot(x,0*x,'-','LineWidth',1.4,'color',[0.2 0.2 0.2])
plot(0*x,x,'-','LineWidth',1.4,'color',[0.2 0.2 0.2])

linewidth = 1;
for pntIdx = 1:length(xData)
    xVal = xData(pntIdx);
    yVal = yData(pntIdx);
    scale = abs(xData-yData);dist = -(xVal-yVal)/sqrt(2)/100;
    dist = real(dist);
    dist(isnan(dist)) = 0; pointColor = [smR(dist),smG(dist),smB(dist)];
    pointColor(pointColor>1)=1;pointColor(isnan(pointColor))=0;
    pointColor(pointColor<0)=0;
    errorbar(xData(pntIdx),yData(pntIdx),yErr(pntIdx),yErr(pntIdx),xErr(pntIdx),xErr(pntIdx),'s','color',pointColor,'MarkerSize',6,'MarkerFaceColor','none','MarkerEdgeColor','none','LineWidth',linewidth)
    s = scatter(xData(pntIdx),yData(pntIdx),20,'color',pointColor,'Marker','o','MarkerFaceColor',pointColor,'MarkerFaceAlpha',0.2,'MarkerEdgeColor',pointColor);
end

xlh = xlabel(['\Delta\fontname{cmss10}v_{\fontname{mwa_cmmi10}R}\fontname{cmss10} [' char(181) 'm/s]']);%,'Interpreter','tex')
ylh = ylabel(['\Delta\fontname{cmss10}v_{\fontname{mwa_cmmi10}H}\fontname{cmss10} [' char(181) 'm/s]']);%,'Interpreter','tex')

signRankTestLine;

ylh.Position(1) = ylh.Position(1) + 180;
xlh.Position(2) = xlh.Position(2) + 100;


xticks([-750 0 750])
yticks([-750 0 750])
xticklabels({-750,'',750})
yticklabels({-750,'',750})

ax.TickLength = [0.03 0.03];
title('$\Delta v$(daugh-moth)','interpreter','latex')

exportgraphics(fig,'OUT-diff_moth_daugh.png','Resolution',600);
saveas(fig,'OUT-diff_moth_daugh.svg');

%% residence time

% identity plot
xData = allRRBCLin;
yData = allHRBCLin;


cond = ~isnan(xData)&~isnan(yData)&~isinf(xData)&~isinf(yData);

xErr = ErrallRRBCLin;
yErr = ErrallHRBCLin;

xErr(isnan(xErr)) = median(xErr,'omitnan');
yErr(isnan(yErr)) = median(yErr,'omitnan');

xErr = xErr(cond);
yErr = yErr(cond);

xData = xData(cond);
yData = yData(cond); if isempty(xData) || isempty(yData), error('x- or y-data empty'), end, if length(xData)<3 || length(yData)<3, error('not enough data points'), end

axLims = [min([xData,yData]) max([xData,yData])];
axLims = [0,5];


fig = figure();
hold on
x = axLims;


% plot linear model
weights = 1./(xErr+yErr); weights = weights/mean(weights);  
mdl = fitlm(xData,yData,'intercept',true,'RobustOpts','on','weights',weights);
lm_slope = mdl.Coefficients.Estimate(2);
lm_intcp = mdl.Coefficients.Estimate(1);
lm_var = mdl.Coefficients.SE(2);
lm_intvar = mdl.Coefficients.SE(1);
xLine = linspace(-1000,1000,2);
lineMid = (lm_slope)*xLine + lm_intcp;
lineTop = (lm_slope+lm_var)*xLine + lm_intcp+lm_intvar;
lineBot = (lm_slope-lm_var)*xLine + lm_intcp-lm_intvar;
fill([xLine,flip(xLine)],[lineTop, fliplr(lineBot)],[0.85 0.85 0.85],'EdgeColor','none')
lineTop = (lm_slope+lm_var)*xLine + lm_intcp-lm_intvar;
lineBot = (lm_slope-lm_var)*xLine + lm_intcp+lm_intvar;
fill([xLine,flip(xLine)],[lineTop, fliplr(lineBot)],[0.85 0.85 0.85],'EdgeColor','none')


plot(x,x,'--','LineWidth',1.4,'color',[0.2 0.2 0.2])
plot(x,0*x,'-','LineWidth',1.4,'color',[0.2 0.2 0.2])
plot(0*x,x,'-','LineWidth',1.4,'color',[0.2 0.2 0.2])

linewidth = 1;
for pntIdx = 1:length(xData)
    xVal = xData(pntIdx);
    yVal = yData(pntIdx);
    scale = abs(xData-yData);dist = -(xVal-yVal)/sqrt(2)/mean(abs(scale(~isnan(scale))));%)-(xVal-yVal)/sqrt(2)/mean(abs(xData-yData));%)-10*(xVal-yVal)/sqrt(2)/mean(abs(xData-yData));%)-100*(xVal-yVal)/sqrt(2)/mean(abs(xData-xData));%)-100*(xVal-yVal)/sqrt(2)/mean(abs(xVal-yVal));%)-10*(xVal-yVal)/sqrt(2)/mean(abs(xVal-yVal));%)-(xVal-yVal)/sqrt(2)/mean(abs(xVal-yVal));%)-(xVal-yVal)/sqrt(2)/(mean(abs(xVal-yVal));%)mean(abs(xVal-yVal));%mean(abs(xVal-yVal));%-(xVal-yVal)/sqrt(2);
    scale = abs(xData-yData);dist = -(xVal-yVal)/sqrt(2)/mean(abs(scale(~isnan(scale))));%)-(xVal-yVal)/sqrt(2)/mean(abs(xData-yData));%)-10*(xVal-yVal)/sqrt(2)/mean(abs(xData-yData));%)-100*(xVal-yVal)/sqrt(2)/mean(abs(xData-xData));%)-100*(xVal-yVal)/sqrt(2)/mean(abs(xVal-yVal));%)-10*(xVal-yVal)/sqrt(2)/mean(abs(xVal-yVal));%)-(xVal-yVal)/sqrt(2)/mean(abs(xVal-yVal));%)-(xVal-yVal)/sqrt(2)/(mean(abs(xVal-yVal));%)mean(abs(xVal-yVal));%mean(abs(xVal-yVal));%dist/50*250;
    dist = real(dist);
    dist(isnan(dist)) = 0; pointColor = [smR(dist),smG(dist),smB(dist)];
    pointColor(pointColor>1)=1;pointColor(isnan(pointColor))=0;
    pointColor(pointColor<0)=0;
    errorbar(xData(pntIdx),yData(pntIdx),yErr(pntIdx),yErr(pntIdx),xErr(pntIdx),xErr(pntIdx),'s','color',pointColor,'MarkerSize',6,'MarkerFaceColor','none','MarkerEdgeColor','none','LineWidth',linewidth)
    s = scatter(xData(pntIdx),yData(pntIdx),20,'color',pointColor,'Marker','o','MarkerFaceColor',pointColor,'MarkerFaceAlpha',0.2,'MarkerEdgeColor',pointColor);
end

xlh = xlabel(['\fontname{cmss10}P_\lambda^{\fontname{mwa_cmmi10}R}']);%,'Interpreter','tex')
ylh = ylabel(['\fontname{cmss10}P_\lambda^{\fontname{mwa_cmmi10}H}']);%,'Interpreter','tex')

signRankTestLine;

ylh.Position(1) = ylh.Position(1) + .6;
xlh.Position(2) = xlh.Position(2) + .5;


xticks([0 2.5 5])
yticks([0 2.5 5])
xticklabels({0,'',5})
yticklabels({0,'',5})

ax.TickLength = [0.02 0.02];
title('$\Delta v$(daugh-moth)','interpreter','latex')

exportgraphics(fig,'OUT-lingering.png','Resolution',600);


saveas(fig,'OUT-lingering.svg');
lingering_distances = (xData-yData)/sqrt(2);