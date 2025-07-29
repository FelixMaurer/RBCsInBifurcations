%--------------------------------------------------------------------------
% Script Name : F3_geometry_extraction.m
% Author      : Felix Maurer
% Institution : Saarland University
% Email       : mail@felixmilanmaurer.com
% Date        : 2024
%
% Description :
%   This script extracts geometrical features from the manually drawn
%   vessel mask. These features are the bifurcation centers, bifurcation
%   borders, marking the beginning or end of attached vessels, the vessel
%   center lines, the vessel walls. It includes a transformation for each
%   vessel into a local coordinate system where the vessel axis is a
%   straight line (x).
%
% Usage :
%   - required input is a hand drawn mask, where inside vessel areas are
%   white and outside vessel areas are black. Open ends of vessels in the
%   mask should be sharply cutoff, perpendicular to the vessel axis. Files
%   should be called 'Mask.png' by default
%
% Dependencies :
%   - Shape Language Modelling by John D'Errico
%   mathworks.com/matlabcentral/fileexchange/24443-slm-shape-language-modeling
% Reference :
%   This script is associated with the publication
%   Impact of Red Blood Cell Rigidity on in vivo Flow Dynamics and Lingering in Bifurcations
%   by Rashidi et al. 2025
% License :
%   MIT
%% source
addpath('src');
%% settings
% these flags are for figures
flag_outputVesselDetection = true;
flag_outputVesselImages = true;
flag_outputBifDetection = true;
flag_outputBifTrafo = true;
%% file locations
addpath('src/SLMtools'); % path to shape language modelling
addpath('src');
rootDir = char(readlines('directory.txt'));
%% find mask files
filelist = dir(fullfile(rootDir, '**\Mask.png'));  % get list of files and folders in any subfolder
filelist = filelist(~[filelist.isdir]);  % remove folders from list
%% analyze geometry
smParam = 0.02;
FLAGplotROIbounds = true;
for fileIdx = 1:length(filelist)
    % file location and declaration
    fileFolder = filelist(fileIdx).folder;
    fileName = filelist(fileIdx).name;
    filePath = [fileFolder '\' fileName];
    outDir = [fileFolder,'\geometry'];
    condition = true; % optional filter conditions
    if condition
        mkdir(outDir);
        % terminal output
        fprintf('-> working on %s\n',filePath);
        wallsImg = (imread(filePath));
        if size(size(wallsImg),2) > 2
            % only take one channel in case of RGB
            wallsImg = wallsImg(:,:,1);
        end
        % take the dimensions
        imgH = size(wallsImg,1);
        imgW = size(wallsImg,2);
        %% image measurements
        % binarize image
        vesselBW = imbinarize(wallsImg);
        rest = true; % optional filter conditions
        if rest
            %% segmentation
            for segIdx = 1:2 % two successive segmentations to improve skeleton
                if segIdx == 1
                    membraneBW = vesselBW;
                    % find skeleton
                    % optional thickening to avoid noisy edges
                    %skelBW = bwskel(bwmorph(membraneBW,'thicken',2));
                    skelBW = bwskel(membraneBW);
                end
                % trace the lines of the skeleton
                % the central boundary ctrBdy contains each vessels axis
                ctrBdy = traceSkel(skelBW);
                % the central boundaries will be also stored in the binary 
                % image ctrBW
                ctrBW = zeros(imgH,imgW,'logical');
                clear smCtrBdy; % smoothed central boundary
                nPnts = zeros(1,length(ctrBdy));
                for ctrIdx = 1:length(ctrBdy)
                    nPnts(ctrIdx) = size(ctrBdy(ctrIdx).pnts,1);
                end
                % loop through all central boundaries
                selIdc = 1:length(ctrBdy);
                selIdc = selIdc(nPnts>20);
                ctrBdy = ctrBdy(selIdc);
                for ctrIdx = 1:length(ctrBdy)
                    ctrPnts = ctrBdy(ctrIdx).pnts;
                    x = ctrPnts(:,1);
                    y = ctrPnts(:,2);
                    extra = 0;
                    pntIdc = (-extra+1:length(x)+extra)';
                    % fit a spline to the points of the line
                    xfit = fit(pntIdc,[x(end-extra+1:end);x;x(1:extra)],'smoothingspline','smoothingparam',smParam);
                    yfit = fit(pntIdc,[y(end-extra+1:end);y;y(1:extra)],'smoothingspline','smoothingparam',smParam);
                    % new smoothed line points
                    xsm = xfit(1:length(x));
                    ysm = yfit(1:length(x));
                    % at each point of the boundary find the normal vector
                    dX = diff(xsm);
                    dY = diff(ysm);
                    % change orientation based on curvature
                    % by choice the normal orientation is radially outwards
                    ddX = diff(dX);
                    ddY = diff(dY);
                    ddX = [ddX; ddX(1)];
                    ddY = [ddY; ddY(1)];
                    curv = mean((dX.*ddY-ddX.*dY)./sqrt(dX.*dX+dY.*dY).^3);
                    if curv >=0
                        normals = [-dY,dX];
                    else
                        normals = [dY,-dX];
                    end
                    normals = normals./sqrt(normals(:,1).^2+normals(:,2).^2);
                    % store all smoothed lines and normals
                    smCtrBdy(ctrIdx).pnts = [xsm,ysm];
                    smCtrBdy(ctrIdx).normals = normals;
                    % write the smoothed boundaries into the binary image
                    for pntIdx = 1:length(xsm)-1
                        ctrBW(round(xsm(pntIdx)),round(ysm(pntIdx))) = 1;
                    end
                end
                %% trace central lines with length parametrization
                % parameterize the lines by length
                clear smParamBdy;
                for ctrIdx =1:length(smCtrBdy)
                    ctrPnts = smCtrBdy(ctrIdx).pnts;
                    x = ctrPnts(:,1);
                    y = ctrPnts(:,2);
                    % compute increments
                    dx = diff(x);
                    dy = diff(y);
                    ds = sqrt(dx.^2+dy.^2);
                    s = cumsum(ds); % path length integral
                    s = [0; s];
                    extra = 0;
                    sfit = [s(end+1-extra:end);s;s(1:extra)];
                    % find coordinates as functions of path length x(s),y(s) 
                    xfit = fit(sfit,[x(end+1-extra:end);x;x(1:extra)],'smoothingspline','smoothingparam',smParam);
                    yfit = fit(sfit,[y(end+1-extra:end);y;y(1:extra)],'smoothingspline','smoothingparam',smParam);

                    dSample = 1;% desired fixed length increment
                    NSamples = max(s)/dSample;
                    fprintf('%d length samples.\n',NSamples);

                    sampleS = linspace(min(s),max(s),NSamples);
                    sampleS = [sampleS,sampleS(1)];
                    xsm = xfit(sampleS);
                    ysm = yfit(sampleS);

                    % at each point of the boundary find the normal vector
                    dX = diff(xsm);
                    dY = diff(ysm);
                    % change orientation based on curvature
                    ddX = diff(dX);
                    ddY = diff(dY);
                    ddX = [ddX; ddX(1)];
                    ddY = [ddY; ddY(1)];
                    curv = mean((dX.*ddY-ddX.*dY)./sqrt(dX.*dX+dY.*dY).^3);
                    if curv >=0
                        normals = [-dY,dX];
                    else
                        normals = [dY,-dX];
                    end
                    normals = normals./sqrt(normals(:,1).^2+normals(:,2).^2);

                    smParamBdy(ctrIdx).pnts = [xsm(1:end-1),ysm(1:end-1)];
                    smParamBdy(ctrIdx).normals = normals;
                    smParamBdy(ctrIdx).length = s;
                end

                close all;
                plot_lines = false; % plot for debugging
                if plot_lines
                    fac = 5;
                    fig = figure;
                    surf(double(vesselBW'),'EdgeColor','none')
                    view(2);
                    colormap('gray')

                    fig.Units = 'centimeters';
                    fig.Position(2) = 3;
                    fig.Position(3) = 10;
                    fig.Position(4) = fig.Position(3)*imgW/imgH;
                    ax = gca;
                    ax.Position(1:2) = 0;
                    ax.Position(3:4) = 1;
                    hold on
                    for memIdx = 1:length(smParamBdy)
                        pnts = smParamBdy(memIdx).pnts;
                        plot3(pnts(:,1),pnts(:,2),ones(1,length(pnts(:,2)))*1000,'color',[0.9 0.9 0.8],'LineWidth',1.4)
                    end
                    colors = lines(length(smParamBdy));
                    for ctrIdx =1:length(smParamBdy)
                        ctrPnts = smParamBdy(ctrIdx).pnts;
                        ctrNorms = smParamBdy(ctrIdx).normals;
                        for pntIdx = 1:size(ctrPnts,1)
                            plot3(ctrPnts(pntIdx,1),ctrPnts(pntIdx,2),1000,'color',colors(ctrIdx,:),'LineStyle','none','Marker','.');
                            plot3(ctrPnts(pntIdx,1)+[0 ctrNorms(pntIdx,1)]*fac,ctrPnts(pntIdx,2)+[0 ctrNorms(pntIdx,2)]*fac,1000*ones(1,2),'-r')
                            plot3(ctrPnts(pntIdx,1)-[0 ctrNorms(pntIdx,1)]*fac,ctrPnts(pntIdx,2)-[0 ctrNorms(pntIdx,2)]*fac,1000*ones(1,2),'-g')
                        end
                    end
                    hold off
                    axis equal
                    toc
                    xlim([1,imgH])
                    ylim([1,imgW])
                end
                if segIdx == 1 
                    %% find corresponding vessels for each bifurcation
                    clc;
                    % find bifurcations
                    skelBW = bwskel(vesselBW);
                    % find points of neighbors 3 meaning bifurcation center
                    numberNeighboringPixels = 3;
                    lut = makelut(@(x)sum(x(:))>=(numberNeighboringPixels+1),3);
                    lutBW = bwlookup(skelBW,lut) & skelBW;
                    stats = regionprops(lutBW,'Centroid');
                    bifPos = [stats.Centroid];
                    % extract bifurcation centers
                    bifCtr = [bifPos(1:2:end)', bifPos(2:2:end)'];
                    allConVesselIdc = [];
                    allConDirection = [];
                    for bifIdx = 1:size(bifCtr,1)
                        cBifCtr = bifCtr(bifIdx,[2 1]);
                        conVesselIdc = [];
                        conDirection = [];
                        for ctrIdx = 1:length(smParamBdy)
                            cPnts = smParamBdy(ctrIdx).pnts;
                            relVec = (cPnts-cBifCtr);
                            distances = sqrt(relVec(:,1).^2+relVec(:,2).^2);
                            if min(distances) < 5
                                conVesselIdc = [conVesselIdc ctrIdx];
                                minIdx = find(distances==min(distances),1);
                                if minIdx < length(distances)/2 % in beginning: direction positive
                                    conDirection = [conDirection 1];
                                else % in the end
                                    conDirection = [conDirection -1];
                                end
                            end
                        end
                        allConVesselIdc = [allConVesselIdc,conVesselIdc];
                        allConDirection = [allConDirection,conDirection];
                    end
                    [allConVesselIdc,sortIdc] = sort(allConVesselIdc);
                    allConDirection = allConDirection(sortIdc);
                    endIdent = logical(diff(allConVesselIdc));
                    endIdent = [endIdent,1];
                    [~,uniqIdc] = unique(allConVesselIdc);
                    endIdent = endIdent(uniqIdc);
                    allConDirection = allConDirection(uniqIdc);
                    %% extend skeleton ends
                    % loop through ends
                    for ctrIdx = find(endIdent)
                        pnts = smParamBdy(ctrIdx).pnts;
                        if allConDirection(ctrIdx) < 0
                            pnts = flip(pnts,1);
                        end
                        dirVec = mean(diff(pnts(end-5:end,:)));
                        dirVec = dirVec/norm(dirVec);
                        kLine = 0;
                        imgPnt = round(pnts(end,:) + dirVec*kLine);
                        runFlag = true;
                        while inbounds(skelBW,imgPnt(:,1),imgPnt(:,2)) && runFlag
                            if membraneBW(imgPnt(:,1),imgPnt(:,2)) == 0
                                runFlag = false;
                            else
                                skelBW(imgPnt(:,1),imgPnt(:,2)) = 1;
                                kLine = kLine + 1;
                                imgPnt = round(pnts(end,:) + dirVec*kLine);
                            end
                        end
                    end
                end
            end
            %% sampling the vessel walls from central axis
            clear profiles
            for ctrIdx = 1:length(smParamBdy)
                ctrPnts = smParamBdy(ctrIdx).pnts;
                ctrNorms = smParamBdy(ctrIdx).normals;
                ctrLength = smParamBdy(ctrIdx).length;
                % sample image profiles
                fac = 6; dFac = 2;
                lineLengthReached = false;
                facFound = false;
                oNumZeros = inf;
                while ~lineLengthReached
                    allSamples = [];
                    for pntIdx = 1:size(ctrPnts,1)
                        linePnts = [ctrPnts(pntIdx,2)+linspace(-ctrNorms(pntIdx,2)*fac,ctrNorms(pntIdx,2)*fac,1+2*ceil(fac));ctrPnts(pntIdx,1)+linspace(-ctrNorms(pntIdx,1)*fac,ctrNorms(pntIdx,1)*fac,1+2*ceil(fac))]';
                        linePnts = round(linePnts); % pixel grid
                        relVec = (linePnts(end,:)-linePnts(1,:));
                        % sample from an image
                        % unify direction
                        samples = zeros(1,size(linePnts,1));
                        for linePntIdx = 1:size(linePnts,1)
                            if inbounds(membraneBW,linePnts(linePntIdx,2),linePnts(linePntIdx,1))
                                samples(linePntIdx) = membraneBW(linePnts(linePntIdx,2),linePnts(linePntIdx,1));
                            end
                        end
                        allSamples = [allSamples,samples'];
                    end
                    % estimate width of vessel
                    unwrap = logical(allSamples);
                    topDiff = (diff(unwrap,1));
                    botDiff = (diff(imcomplement(unwrap),1));
                    topDiff(round(size(unwrap,1)/2):end,:) = 0;
                    botDiff(round(1:size(unwrap,1)/2),:) = 0;
                    topDiff(topDiff<0) = 0;
                    botDiff(botDiff<0) = 0;
                    topDiff = logical(topDiff);
                    botDiff = logical(botDiff);
                    [yTop,xTop] = find(topDiff);
                    [yBot,xBot] = find(botDiff);
                    widths = zeros(1,size(unwrap,2));

                    for idx = 1:length(xTop)
                        yTopVal = yTop(idx);
                        diffArray = abs(xBot-xTop(idx));
                        minIdx = find(diffArray == min(diffArray));
                        yBotVal = yBot(minIdx);
                        minIdx = minIdx(yBotVal==min(yBotVal));
                        if diffArray(minIdx) == 0
                            widths(xTop(idx)) = yBot(minIdx)-yTop(idx);
                        end
                    end

                    numZeros = sum(widths==0);

                    if facFound
                        profiles(ctrIdx).samples = allSamples;
                        lineLengthReached = true;
                    else
                        if numZeros > oNumZeros || numZeros == 0 || sum(allSamples(1,:)+allSamples(end,:)) == 0 || (numZeros == oNumZeros && numZeros < length(widths)/2)
                            facFound = true;
                        else
                            oNumZeros = numZeros;
                            fac = fac + dFac;
                        end
                    end
                end
            end
            fprintf('line length found with factor %.f.\n',fac);
            %% find corresponding vessels for each bifurcation
            % again with final layout after all segmentation steps
            clc;
            % find bifurcations
            skelBW = bwskel(vesselBW);
            % find points of neighbors 3 meaning bifurcation center
            numberNeighboringPixels = 3;
            lut = makelut(@(x)sum(x(:))>=(numberNeighboringPixels+1),3);
            lutBW = bwlookup(skelBW,lut) & skelBW;
            stats = regionprops(lutBW,'Centroid');
            bifPos = [stats.Centroid];
            % extract bifurcation centers
            bifCtr = [bifPos(1:2:end)', bifPos(2:2:end)'];

            clear bifurcations;

            for bifIdx = 1:size(bifCtr,1)
                cBifCtr = bifCtr(bifIdx,[2 1]);
                conVesselIdc = [];
                conDirection = [];
                for ctrIdx = 1:length(smParamBdy)
                    cPnts = smParamBdy(ctrIdx).pnts;
                    relVec = (cPnts-cBifCtr);
                    distances = sqrt(relVec(:,1).^2+relVec(:,2).^2);
                    if min(distances) < 5
                        conVesselIdc = [conVesselIdc ctrIdx];
                        minIdx = find(distances==min(distances),1);
                        if minIdx < length(distances)/2 % in beginning: direction positive
                            conDirection = [conDirection 1];
                        else % in the end
                            conDirection = [conDirection -1];
                        end
                    end
                end
                bifurcations(bifIdx).bifCtr = cBifCtr;
                for vesIdx = 1:length(conVesselIdc)
                    pnts = smParamBdy(conVesselIdc(vesIdx)).pnts;
                    normals = smParamBdy(conVesselIdc(vesIdx)).normals;
                    unwrap = profiles(conVesselIdc(vesIdx)).samples;
                    vesLength = smParamBdy(conVesselIdc(vesIdx)).length;
                    if conDirection(vesIdx) < 0
                        pnts = flip(pnts,1);
                        normals = flip(normals,1);
                        unwrap = flip(unwrap,2);
                        %unwrap = flip(unwrap,1);
                    else
                        %normals = -normals;
                        %unwrap = flip(unwrap,1);
                    end
                    unwrap = bwareaopen(unwrap,100);
                    % remove disconnected
                    unwrap = logical(unwrap - imcomplement(imfill(imcomplement(unwrap),round([size(unwrap,1)/2 size(unwrap,2)/2]))));
                    % remove zeros from border
                    startIdx = find(sum(unwrap,2)>0,1,'first');
                    endIdx = find(sum(unwrap,2)>0,1,'last');
                    roi = size(unwrap,1)/2+[-1,1]*max(abs([startIdx,endIdx]-size(unwrap,1)/2))+[1 0];
                    unwrap = unwrap(roi(1):roi(2),:);
                    % save everything
                    bifurcations(bifIdx).conVesBdy(vesIdx).pnts = pnts;
                    bifurcations(bifIdx).conVesBdy(vesIdx).normals = normals;
                    bifurcations(bifIdx).conVesBdy(vesIdx).unwrap = unwrap;
                    bifurcations(bifIdx).conVesBdy(vesIdx).length = vesLength;
                end
            end
            %% sort bifurcation vessels clock-wise
            for bifIdx = 1:length(bifurcations)
                bifImgPnts = [];
                conVesBdy = bifurcations(bifIdx).conVesBdy;
                orient = zeros(1,length(conVesBdy));
                for vesIdx = 1:length(conVesBdy)
                    pnts = conVesBdy(vesIdx).pnts;
                    dirVec = mean(diff(pnts(1:5,:)));
                    dirVec = dirVec./norm(dirVec);
                    if vesIdx == 1
                        refVec = dirVec';
                    end
                    refNorm = [refVec(2);-refVec(1)];
                    measVec = [dirVec.*refVec,dirVec.*refNorm];
                    % measure angle
                    orient(vesIdx) = (measVec(1)+1)*(1-2*double(measVec(2)<0));
                end
                [~,sortIdc] = sort(orient);
                bifurcations(bifIdx).conVesBdy = conVesBdy(flip(sortIdc));
            end
            %% correct all normal orientations
            for bifIdx = 1:length(bifurcations)
                conVesBdy = bifurcations(bifIdx).conVesBdy;
                for vesIdx = 1:length(conVesBdy)
                    pnts = conVesBdy(vesIdx).pnts;
                    norms = conVesBdy(vesIdx).normals;
                    unwrap = conVesBdy(vesIdx).unwrap;
                    dirVec = (pnts(3,:)-pnts(1,:));
                    dirVec = dirVec/norm(dirVec);
                    if [dirVec(2) -dirVec(1)]*norms(3,:)' > 0
                        % change orientation
                        norms = -norms;
                        unwrap = flip(unwrap,1);
                    end
                    bifurcations(bifIdx).conVesBdy(vesIdx).normals = norms;
                    bifurcations(bifIdx).conVesBdy(vesIdx).unwrap = unwrap;
                end
            end
            %% take out ends
            take_out_ends = true;
            figure
            hold on
            if take_out_ends
                for bifIdx = 1:length(bifurcations)
                    conVes = bifurcations(bifIdx).conVesBdy;
                    for vesIdx = 1:length(conVes)
                        unwrap = conVes(vesIdx).unwrap;
                        pnts = conVes(vesIdx).pnts;
                        normals = conVes(vesIdx).normals;
                        vesLength = conVes(vesIdx).length;
                        %imshow(unwrap)
                        topDiff = (diff(unwrap,1));
                        botDiff = (diff(imcomplement(unwrap),1));
                        topDiff(round(size(unwrap,1)/2):end,:) = 0;
                        botDiff(round(1:size(unwrap,1)/2),:) = 0;
                        topDiff(topDiff<0) = 0;
                        botDiff(botDiff<0) = 0;
                        topDiff = logical(topDiff);
                        botDiff = logical(botDiff);
                        [yTop,xTop] = find(topDiff);
                        [yBot,xBot] = find(botDiff);
                        widths = zeros(1,size(unwrap,2));
                        symmetry = zeros(1,size(unwrap,2));
                        for idx = 1:length(xTop)
                            yTopVal = yTop(idx);
                            diffArray = abs(xBot-xTop(idx));
                            minIdx = find(diffArray == min(diffArray));
                            yBotVal = yBot(minIdx);
                            minIdx = minIdx(yBotVal==min(yBotVal));
                            if diffArray(minIdx) == 0
                                widths(xTop(idx)) = yBot(minIdx)-yTop(idx);
                                symmetry(xTop(idx)) = (yBot(minIdx)+yTop(idx))-size(unwrap,1);
                            end
                        end
                        widths(widths==0) = size(unwrap,1);
                        proStartIdx = round(length(widths)/2);
                        profile = movmean(widths(proStartIdx:end),3);
                        smFit = fit((1:length(profile))',profile','smoothingspline','smoothingparam',smParam);
                        smProfile = smFit(1:length(profile));
                        [peaks,peakPos] = findpeaks(smProfile,'MinPeakProminence',3);
                        endIdx = length(widths); % in case there is no corner mistake take last
                        % no peak case
                        idc = (-20:-1)+length(widths)+1;
                        linfitEnd = fit(idc',widths(idc)','poly1');
                        idc = (round(-length(widths)*0.2):round(length(widths)*0.2))+round(length(widths)*0.4);

                        linfitMid = fit(idc',widths(idc)','poly1');

                        dirEnd = [1 linfitEnd.p1];
                        dirMid = [1 linfitMid.p1];
                        dirEnd = dirEnd/norm(dirEnd);
                        dirMid = dirMid/norm(dirMid);

                        xInt = length(widths);
                        if dirEnd*dirMid' < 0.95 && linfitEnd.p1 < 0
                            % corner mistake
                            % find intersect
                            xInt = round((linfitEnd.p2-linfitMid.p2)/((linfitMid.p1-linfitEnd.p1)));
                        end

                        if ~isempty(peakPos)
                            % there is a corner mistake
                            endIdx = proStartIdx+max(peakPos)-10; % subtract to get away from turn point
                        end
                        % check symmetry
                        symEnd = find(abs(symmetry(proStartIdx:end))>300,1,'first');
                        if isempty(symEnd)
                            symEnd = inf;
                        end
                        endIdx = min(endIdx,symEnd);
                        endIdx = min(endIdx,xInt);
                        % save widths, zero width means unknown
                        widths = widths(1:endIdx);
                        unwrap = unwrap(:,1:endIdx);
                        pnts = pnts(1:endIdx,:);
                        normals = normals(1:endIdx,:);
                        bifurcations(bifIdx).conVesBdy(vesIdx).widths = widths;
                        bifurcations(bifIdx).conVesBdy(vesIdx).pnts = pnts;
                        bifurcations(bifIdx).conVesBdy(vesIdx).normals = normals;
                        bifurcations(bifIdx).conVesBdy(vesIdx).unwrap = unwrap;
                    end
                end
            end
            hold off
            %% find bifurcation area
            for bifIdx = 1:length(bifurcations)
                cVesselBdy = bifurcations(bifIdx).conVesBdy;
                for vesIdx = 1:length(cVesselBdy)
                    close all
                    unwarpOr = cVesselBdy(vesIdx).unwrap;
                    cUnWrap = unwarpOr;
                    cUnWrap(:,1) = 1;
                    cUnWrap(:,end) = 1;
                    cUnWrap = cUnWrap(:,1:round(size(cUnWrap,2)*3/4));
                    edgeImg = bwskel(imfill(bwskel(logical(cUnWrap - imerode(cUnWrap,strel('disk',1,4))+edge(cUnWrap))),'holes'));
                    topHalf = edgeImg(1:round(size(cUnWrap,1)/2),:);
                    botHalf = edgeImg(round(size(cUnWrap,1)/2):end,:);

                    clear LinePnts;
                    clear halfs;
                    halfs(1).img = topHalf;
                    halfs(2).img = botHalf;
                    for halfIdx = 1:2
                        skelBWTrace = halfs(halfIdx).img;
                        [y,x] = find(skelBWTrace);
                        % find points of neighbors 1 meaning open end
                        numberNeighboringPixels = 1;
                        lut = makelut(@(x)sum(x(:))==(numberNeighboringPixels+1),3);
                        lutBW = bwlookup(skelBWTrace,lut) & skelBWTrace;
                        stats = regionprops(lutBW,'Centroid');
                        endPos = [stats.Centroid];
                        endX = endPos(1:2:end);
                        endY = endPos(2:2:end);
                        endPnt = [endX',endY'];
                        startPnt = endPnt(find(endPnt(:,1)==min(endPnt(:,1)),1),:);

                        thisEndPnt = startPnt;
                        cPnt = thisEndPnt;
                        linePoints = [];
                        skelDistPnts = [x,y];
                        distVec = cPnt-skelDistPnts;
                        dists = sqrt(distVec(:,1).^2+distVec(:,2).^2);
                        pntIdx = find(dists==min(dists),1,'first');
                        % take out this point
                        cond = ones(1,size(skelDistPnts,1),'logical');
                        cond(pntIdx) = 0;
                        skelDistPnts = skelDistPnts(cond,:);
                        % save
                        linePoints = [linePoints;cPnt];
                        k_nn = 0;
                        minDist = 1;
                        while minDist < 2
                            k_nn = k_nn+1;
                            distVec = cPnt-skelDistPnts;
                            dists = sqrt(distVec(:,1).^2+distVec(:,2).^2);
                            minDist = min(dists);
                            pntIdx = find(dists==min(dists),1,'first');
                            % determine new point
                            cPnt = skelDistPnts(pntIdx,:);
                            % take out this point
                            cond = ones(1,size(skelDistPnts,1),'logical');
                            cond(pntIdx) = 0;
                            skelDistPnts = skelDistPnts(cond,:);
                            % save
                            linePoints = [linePoints;cPnt];
                        end
                        %plot(linePoints(:,1),linePoints(:,2),'.')
                        linePoints(:,2) = linePoints(:,2) + (halfIdx==2)*(round(size(cUnWrap,1)/2)-1);
                        LinePnts(halfIdx).pnts = flip(linePoints,2);
                    end

                    xTop = LinePnts(1).pnts(:,2);
                    yTop = LinePnts(1).pnts(:,1);
                    xBot = LinePnts(2).pnts(:,2);
                    yBot = LinePnts(2).pnts(:,1);

                    xTop = flip(xTop,1);
                    yTop = flip(yTop,1);
                    xBot = flip(xBot,1);
                    yBot = flip(yBot,1);

                    maxIdx = min(size(xTop,1),size(xBot,1));
                    xTop = xTop(1:maxIdx);
                    yTop = yTop(1:maxIdx);
                    xBot = xBot(1:maxIdx);
                    yBot = yBot(1:maxIdx);

                    conVecs = diff([xTop,yTop]);
                    % normalize
                    conVecs = conVecs./sqrt(conVecs(:,1).^2+conVecs(:,2).^2);
                    conVecsTop = conVecs;

                    conVecs = diff([xBot,yBot]);
                    % normalize
                    conVecs = conVecs./sqrt(conVecs(:,1).^2+conVecs(:,2).^2);
                    conVecsBot = conVecs;

                    orth = zeros(1,size(conVecsTop,1));
                    for conVecIdx = 1:size(conVecsTop,1)
                        orth(conVecIdx) = conVecsTop(conVecIdx,:)*conVecsBot(conVecIdx,:)';
                    end

                    bifEndIdx = find(orth>0.95,1,'last');
                    bifEndPos = [xTop(bifEndIdx),yTop(bifEndIdx);xBot(bifEndIdx),yBot(bifEndIdx)];
                    outImg = cVesselBdy(vesIdx).unwrap;
                    % save in struct
                    bifurcations(bifIdx).conVesBdy(vesIdx).bifEndPos = bifEndPos;
                    bifurcations(bifIdx).conVesBdy(vesIdx).bifEndIdx = bifEndIdx;
                    bifurcations(bifIdx).conVesBdy(vesIdx).topPnts = [xTop,yTop];%[xBot,size(cUnWrap,1)-yBot];
                    bifurcations(bifIdx).conVesBdy(vesIdx).botPnts = [xBot,yBot];%[xTop,size(cUnWrap,1)-yTop];
                end
            end
            % filter vessels
            for bifIdx = 1:length(bifurcations)
                conVes = bifurcations(bifIdx).conVesBdy;
                selIdc = ones(1,length(conVes),'logical');
                for vesIdx = 1:length(conVes)
                    if isempty(conVes(vesIdx).bifEndPos)
                        selIdc(vesIdx) = 0;
                    end
                end
                bifurcations(bifIdx).conVesBdy = conVes(selIdc);
            end
            close all

            %% find bifurcation image points
            for bifIdx = 1:length(bifurcations)
                bifImgPnts = [];
                conVesBdy = bifurcations(bifIdx).conVesBdy;
                for vesIdx = 1:length(conVesBdy)
                    bifEndPos = conVesBdy(vesIdx).bifEndPos;
                    bifEndIdx = conVesBdy(vesIdx).bifEndIdx;
                    pnts = conVesBdy(vesIdx).pnts;
                    bifEndPos = min(bifEndPos,size(pnts,1));
                    norms = conVesBdy(vesIdx).normals;
                    unwrap = conVesBdy(vesIdx).unwrap;
                    startPnt = pnts(bifEndPos(:,1),:)';
                    startNorm = norms(bifEndPos(:,1),:)';
                    bifPnt1 = startPnt(:,1)+startNorm(:,1)*(bifEndPos(1,2)-round(size(unwrap,1)/2));
                    bifPnt2 = startPnt(:,2)+startNorm(:,2)*(bifEndPos(2,2)-round(size(unwrap,1)/2));
                    bifImgPnts = [bifImgPnts,bifPnt1,bifPnt2];
                end
                bifurcations(bifIdx).bifImgPnts = bifImgPnts;
            end
            selIdc = ones(1,length(bifurcations),'logical');
            for bifIdx = 1:length(bifurcations)
                bifImgPnts = bifurcations(bifIdx).bifImgPnts;
                if isempty(bifImgPnts)
                    selIdc(bifIdx) = 0;
                end
            end
            bifurcations = bifurcations(selIdc);
            %% correct bifurcation point overlap
            fprintf('correcting overlap...');
            allBifEnds = [];
            for bifIdx = 1:length(bifurcations)
                conVes = bifurcations(bifIdx).conVesBdy;
                for vesIdx = 1:length(conVes)
                    bifEndPos = conVes.bifEndPos;
                    allBifEnds  = [allBifEnds;bifEndPos];
                end
                bifImgPnts = bifurcations(bifIdx).bifImgPnts';
                % boundary
                bdyIdc = boundary(bifImgPnts,1);
                if length(bdyIdc) > max(bdyIdc)+1
                    bdyIdc = boundary(bifImgPnts,0);
                end
                try
                    bdyIdc(end) = length(bdyIdc);
                    [~,sorti] = sort(bdyIdc);
                    if sorti(2) == 2
                        pntList = [bdyIdc(2:2:length(bdyIdc)),...
                            bdyIdc(3:2:length(bdyIdc))];
                    else
                        pntList = [bdyIdc(3:2:length(bdyIdc)),...
                            bdyIdc([4:2:length(bdyIdc) 2])];
                    end
                    signList = find(pntList(:,1)-pntList(:,2)>0)';
                    pntList(pntList==length(bdyIdc)) = 1;
                    % go through and flip
                    if ~isempty(signList)
                        for listIdx = signList
                            pnt1 = bifImgPnts(pntList(listIdx,1),:);
                            pnt2 = bifImgPnts(pntList(listIdx,2),:);
                            offDist = norm(pnt2-pnt1)/2;

                            % first vessel of swap pair
                            bifPntIdx = 1;
                            vesIdx = listIdx;
                            oldBifEnd1 = conVes(vesIdx).bifEndPos(bifPntIdx,:);
                            newBifEnd1X = round(oldBifEnd1(1)+offDist);
                            % top
                            topPnts = conVes(vesIdx).topPnts;
                            newXIdx = find(topPnts(:,1)==newBifEnd1X);
                            newXIdx = newXIdx(find(topPnts(newXIdx,2) == max(topPnts(newXIdx,2)),1,'first'));
                            newBifEnd1Y = topPnts(newXIdx,2);
                            % bottom
                            botPnts = conVes(vesIdx).botPnts;
                            newXIdx = find(botPnts(:,1)==newBifEnd1X);
                            newXIdx = newXIdx(find(botPnts(newXIdx,2) == min(botPnts(newXIdx,2)),1,'first'));
                            newBifEnd2Y = botPnts(newXIdx,2);
                            newBifEnd2X = newBifEnd1X;
                            % save
                            conVes(vesIdx).bifEndPos(1,:) = [newBifEnd1X,newBifEnd1Y];
                            conVes(vesIdx).bifEndPos(2,:) = [newBifEnd2X,newBifEnd2Y];

                            % second vessel of swap pair
                            bifPntIdx = 2;
                            vesIdx = listIdx;
                            oldBifEnd1 = conVes(vesIdx).bifEndPos(bifPntIdx,:);
                            newBifEnd1X = round(oldBifEnd1(1)+offDist);
                            % top
                            topPnts = conVes(vesIdx).topPnts;
                            newXIdx = find(topPnts(:,1)==newBifEnd1X);
                            newXIdx = newXIdx(find(topPnts(newXIdx,2) == max(topPnts(newXIdx,2)),1,'first'));
                            newBifEnd1Y = topPnts(newXIdx,2);
                            % bottom
                            botPnts = conVes(vesIdx).botPnts;
                            newXIdx = find(botPnts(:,1)==newBifEnd1X);
                            newXIdx = newXIdx(find(botPnts(newXIdx,2) == min(botPnts(newXIdx,2)),1,'first'));
                            newBifEnd2Y = botPnts(newXIdx,2);
                            newBifEnd2X = newBifEnd1X;
                            % save
                            conVes(vesIdx).bifEndPos(1,:) = [newBifEnd1X,newBifEnd1Y];
                            conVes(vesIdx).bifEndPos(2,:) = [newBifEnd2X,newBifEnd2Y];
                        end
                    end
                    bifurcations(bifIdx).conVesBdy = conVes;
                end
                %% show status
                StatusBar(bifIdx, length(bifurcations), 'bifurcation');
            end
            %% plot unwrapped and bifurcation ends
            plot_unwrapped = false;
            if flag_outputVesselImages || flag_outputVesselDetection
                plot_unwrapped = true;
            end
            if plot_unwrapped
                fig = figure;
                surf(double(flip(membraneBW,1)),'EdgeColor','none')
                view(2);
                colormap('gray')
                xlim([1,size(membraneBW,2)])
                ylim([1,size(membraneBW,1)])
                fig.Units = 'centimeters';
                fig.Position(2) = 3;
                fig.Position(3) = 18;
                fig.Position(4) = fig.Position(3)*imgH/imgW;
                ax = gca;
                ax.Position(1:2) = 0;
                ax.Position(3:4) = 1;
                for bifIdx = 1:length(bifurcations)
                    conVes = bifurcations(bifIdx).conVesBdy;
                    % colors for each vessel
                    colors = lines(length(conVes));
                    % plot vessels
                    hold on
                    for vesIdx = 1:length(conVes)
                        ctrPnts = conVes(vesIdx).pnts;
                        unwrap = conVes(vesIdx).unwrap;
                        plot3(ctrPnts(:,2),size(membraneBW,1)+1-ctrPnts(:,1),ones(1,size(ctrPnts,1))*1000,'color',colors(vesIdx,:),'LineWidth',1.6)
                        ax = gca;
                        pos = ax.Position;
                        yPosition = size(membraneBW,1)-ctrPnts(round(length(ctrPnts)/2),1);
                        xPosition = ctrPnts(round(length(ctrPnts)/2),2);
                        posx = (xPosition-min(xlim))/(max(xlim)-min(xlim))*pos(3)+pos(1);
                        posy = (yPosition-min(ylim))/(max(ylim)-min(ylim))*pos(4)+pos(2);
                        annotation('textbox','Position',[posx,posy,0,0],'String',sprintf('b%02dv%02d',bifIdx,vesIdx),'EdgeColor','y','color',colors(vesIdx,:),'FontSize',9)
                    end
                end
                hold off
                outPath = [outDir,'\','vesselDetection.png'];
                if flag_outputVesselDetection
                    print(gcf,outPath,'-dpng');
                end
                %% plot bifurcation ends
                fprintf('vessel coordinate transformation...');
                for bifIdx = 1:length(bifurcations)
                    conVes = bifurcations(bifIdx).conVesBdy;
                    for vesIdx = 1:length(conVes)
                        bifEndPos = conVes(vesIdx).bifEndPos;
                        close all;
                        memUnWrap = double(conVes(vesIdx).unwrap);
                        padding = 20;
                        memUnWrap = padarray(memUnWrap,padding,0);
                        fig = figure;
                        surf((memUnWrap),'EdgeColor','none')
                        view(2);
                        colormap('gray')
                        clear lengths;
                        hold on
                        colors = lines(length(conVes));
                        plot3([1 size(memUnWrap,2)],[1,1],[1000,1000],'color',colors(vesIdx,:),'LineWidth',10);
                        plot3([1 size(memUnWrap,2)],[padding+1,padding+1],[1000,1000],'color',[0.8 0.8 0.8],'LineWidth',1,'LineStyle','--');
                        plot3([1 size(memUnWrap,2)],size(memUnWrap,1)-[padding,padding]+1,[1000,1000],'color',[0.8 0.8 0.8],'LineWidth',1,'LineStyle','--');
                        plot3(bifEndPos(:,1),padding+bifEndPos(:,2)+[0;1],[1000,1000],'color',[1 0 0],'Marker','x','MarkerSize',13);
                        hold off
                        xlim([1,size(memUnWrap,2)])
                        ylim([1,size(memUnWrap,1)])
                        fig.Units = 'centimeters';
                        fig.Position(2) = 3;
                        fig.Position(3) = 10;
                        fig.Position(4) = fig.Position(3)*size(memUnWrap,1)/size(memUnWrap,2);
                        ax = gca;
                        ax.Position(1:2) = 0;
                        ax.Position(3:4) = 1;

                        pos = ax.Position;

                        xPosition = size(memUnWrap,2)/2;
                        yPosition = 0;

                        posx = (xPosition-min(xlim))/(max(xlim)-min(xlim))*pos(3)+pos(2);
                        posy = (yPosition-min(ylim))/(max(ylim)-min(ylim))*pos(4)+pos(1);
                        annotation('textbox','Position',[posx,posy,1,1],'String',sprintf('b%02d v%02d',bifIdx,vesIdx),'EdgeColor','none','color',colors(vesIdx,:),'FontSize',12)
                        outPath = [outDir,'\',sprintf('bif%03d_ves%03d.png',bifIdx,vesIdx)];
                        if flag_outputVesselImages
                            print(gcf,outPath,'-dpng');
                        end
                    end
                    %% show status
                    StatusBar(bifIdx, length(bifurcations), 'bifurcation');
                end
            end
            %% go through bifurcations and draw boundaries
            fprintf('drawing boundaries...');
            close all
            imshow(membraneBW)
            hold on
            for bifIdx = 1:length(bifurcations)
                bifImgPnts = [];
                conVesBdy = bifurcations(bifIdx).conVesBdy;
                if length(conVesBdy) <3
                    lineColor = [1 0.6 0];
                else
                    lineColor = [0 0.8 0.3];
                end
                for vesIdx = 1:length(conVesBdy)
                    bifEndPos = conVesBdy(vesIdx).bifEndPos;
                    bifEndIdx = conVesBdy(vesIdx).bifEndIdx;
                    unwrap= conVesBdy(vesIdx).unwrap;
                    pnts = conVesBdy(vesIdx).pnts;
                    bifEndPos = min(bifEndPos,size(pnts,1));
                    norms = conVesBdy(vesIdx).normals;
                    startPnt = pnts(bifEndPos(:,1),:)';
                    startNorm = norms(bifEndPos(:,1),:)';
                    bifPnt1 = startPnt(:,1)+startNorm(:,1)*(bifEndPos(1,2)-round(size(unwrap,1)/2));
                    bifPnt2 = startPnt(:,2)+startNorm(:,2)*(bifEndPos(2,2)-round(size(unwrap,1)/2));
                    plot([bifPnt1(2),bifPnt2(2)],[bifPnt1(1),bifPnt2(1)],'-','LineWidth',2,'color',lineColor);
                    plot(bifPnt1(2),bifPnt1(1),'.r','MarkerSize',12);
                    plot(bifPnt2(2),bifPnt2(1),'.r','MarkerSize',12);
                    ax = gca;
                    ax.Position(3:4) = 1;
                    ax.Position(1:2) = 0;
                    bifImgPnts = [bifImgPnts,bifPnt1,bifPnt2];
                end
                bifurcations(bifIdx).bifImgPnts = bifImgPnts;
                %% show status
                StatusBar(bifIdx, length(bifurcations), 'bifurcation');
            end
            hold off
            outPath = [outDir,'\','bifurcationDetection.png'];
            if flag_outputBifDetection
                print(gcf,outPath,'-dpng');
            end
            %% bifurcation trafo
            fprintf('bifurcation coordinate transformation...');
            close all
            for bifIdx = 1:length(bifurcations)
                allTrafoPnts = [];
                % find points inside bifurcation
                bifImgPnts = bifurcations(bifIdx).bifImgPnts;
                bifBound = boundary(bifImgPnts',0.0);
                [x,y] = find(membraneBW);


                idc = inpolygon(y,x,bifImgPnts(2,bifBound),bifImgPnts(1,bifBound));
                bifInnerPnts = [x(idc),y(idc)];
                bifurcations(bifIdx).bifInnerPnts = bifInnerPnts;
                bifurcations(bifIdx).bifBound = bifBound;

                conVesBdy = bifurcations(bifIdx).conVesBdy;
                allTrafoPnts = [allTrafoPnts;bifInnerPnts];
                for vesIdx = 1:length(conVesBdy)
                    bifEndPos = conVesBdy(vesIdx).bifEndPos;
                    bifEndIdx = conVesBdy(vesIdx).bifEndIdx;
                    unwrap= conVesBdy(vesIdx).unwrap;
                    pnts = conVesBdy(vesIdx).pnts;
                    bifEndPos = min(bifEndPos,size(pnts,1));
                    norms = conVesBdy(vesIdx).normals;
                    startPnt = pnts(bifEndPos(:,1),:)';
                    startNorm = norms(bifEndPos(:,1),:)';
                    bifPnt1 = startPnt(:,1)+startNorm(:,1)*(bifEndPos(1,2)-round(size(unwrap,1)/2));
                    bifPnt2 = startPnt(:,2)+startNorm(:,2)*(bifEndPos(2,2)-round(size(unwrap,1)/2));
                    % first remove parts of bifurcation
                    startY = max(bifEndPos(1,2)-2,1);
                    endY = min(bifEndPos(2,2)+2,size(unwrap,1));

                    yIdc = startY:endY;
                    xLine = max(min(bifEndPos(:,1))-1,2);
                    unwrap(yIdc,ones(1,length(yIdc))*xLine) = 0;
                    unwrap = imcomplement(imfill(imcomplement(unwrap),[round(size(unwrap,1)/2),xLine-1]));
                    % take points in unwrapped
                    [x,y] = find(unwrap);

                    % take tangential vector
                    relVec = bifPnt2-bifPnt1;
                    dirVec = -[relVec(2),-relVec(1)];
                    dirVec = dirVec/norm(dirVec);
                    normVec = relVec'/norm(relVec);
                    if (pnts(3,:)-pnts(1,:))*dirVec' < 0
                        dirVec = -dirVec;
                    end
                    startLength = mean(bifEndPos(:,1));
                    x = x(y>startLength);
                    y = y(y>startLength);

                    connectPnt = mean(startPnt,2)' - startLength*dirVec;
                    vesInnerPnts = connectPnt + y*dirVec + (x-size(unwrap,1)/2)*normVec;
                    allTrafoPnts = [allTrafoPnts;vesInnerPnts];
                end
                % make image
                trafoOffset =  - min(allTrafoPnts);
                allTrafoPnts = allTrafoPnts + trafoOffset;
                bifImgH = ceil(max(allTrafoPnts(:,1)));
                bifImgW = ceil(max(allTrafoPnts(:,2)));
                bifImg = zeros(bifImgH,bifImgW,'logical');
                for imgPntIdx = 1:size(allTrafoPnts,1)
                    try
                        bifImg(ceil(allTrafoPnts(imgPntIdx,1)),ceil(allTrafoPnts(imgPntIdx,2))) = 1;
                        bifImg(floor(allTrafoPnts(imgPntIdx,1)),floor(allTrafoPnts(imgPntIdx,2))) = 1;
                    end
                end
                % fill holes
                bifImg = bwmorph(bifImg,'close');
                close all
                fig = figure;

                imshow(bifImg)

                hold on
                plot(bifImgPnts(2,bifBound)+trafoOffset(2),bifImgPnts(1,bifBound)+trafoOffset(1),'LineWidth',1,'LineStyle','--','color',[0 0.5 0.8])
                hold off
                xlim([1,size(bifImg,2)])
                ylim([1,size(bifImg,1)])
                fig.Units = 'centimeters';
                fig.Position(2) = 3;
                fig.Position(3) = 10;
                fig.Position(4) = fig.Position(3)*bifImgH/bifImgW;
                ax = gca;
                ax.Position(1:2) = 0;
                ax.Position(3:4) = 1;
                if flag_outputBifTrafo
                    outPath = [outDir,'\',sprintf('Bif%03dTrafo.png',bifIdx)];
                    print(gcf,outPath,'-dpng');
                end
                %% show status
                StatusBar(bifIdx, length(bifurcations), 'bifurcation');
            end
            % save variables
            outPath = [outDir,'\','bifurcations.mat'];
            save(outPath,'bifurcations');
        end
    end
    %% show status
    StatusBar(fileIdx, length(filelist), 'file');  
end
