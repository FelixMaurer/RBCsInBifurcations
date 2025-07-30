%--------------------------------------------------------------------------
% Script Name : traceSkel.m
% Authors     : Felix Maurer
% Institution : Saarland University
% Email       : mail@felixmilanmaurer.com
% Date        : 2024
%
% Description :
%   This function traces skeletons of binary images and returns the lines
%   in a struct.
function lines = traceSkel(skelBW)
numberNeighboringPixels = 3;
lut = makelut(@(x)sum(x(:))>=(numberNeighboringPixels+1),3);
lutBW = bwlookup(skelBW,lut) & skelBW;
stats = regionprops(lutBW,'Centroid');
bifPos = [stats.Centroid];
% take out bifurcations
skelBWTrace = skelBW;
bifX = bifPos(1:2:end);
bifY = bifPos(2:2:end);
for bifIdx = 1:length(bifX)
    for idx1 = -1:1
        for idx2 = -1:1
            skelBWTrace(round(bifY(bifIdx))+idx1,round(bifX(bifIdx))+idx2) = 0;
        end
    end
end
% find points of neighbors 1 meaning open end
numberNeighboringPixels = 1;
lut = makelut(@(x)sum(x(:))==(numberNeighboringPixels+1),3);
lutBW = bwlookup(skelBWTrace,lut) & skelBWTrace;
stats = regionprops(lutBW,'Centroid');
endPos = [stats.Centroid];
endX = endPos(1:2:end);
endY = endPos(2:2:end);
endPnt = [endX',endY'];

[y,x] = find(skelBWTrace);
skelPnts = [x,y];
% if there is no endpoint choose random
if isempty(endPnt)
    endPnt = skelPnts(1,:);
end
stats = regionprops(skelBWTrace,'PixelList');
for statsIdx = 1:length(stats)
    pnts = stats(statsIdx).PixelList;
    % find end point that is in this line
    endPntFound = false;
    for endPntIdx = 1:size(endPnt,1)
        thisEndPnt = endPnt(endPntIdx,:);
        relVec =  round(thisEndPnt)  -   round(pnts) ;
        dists = sqrt(relVec(:,1).^2+relVec(:,2).^2);
        minIdx = find(dists==min(dists),1);
        if dists(minIdx) < 1.6
            endPntFound = true;
            break
        end
    end
    % if there was no endpoint in this stat choose random
    if ~endPntFound
        thisEndPnt = pnts(1,:);
        endPntFound = true;
    end
    linePoints = [];
    cPnt = thisEndPnt;

    skelDistPnts = pnts;
    distVec = cPnt-skelDistPnts;
    dists = sqrt(distVec(:,1).^2+distVec(:,2).^2);
    pntIdx = find(dists==min(dists),1,'first');
    % take ot this point
    cond = ones(1,size(skelDistPnts,1),'logical');
    cond(pntIdx) = 0;
    skelDistPnts = skelDistPnts(cond,:);
    % save
    linePoints = [linePoints;cPnt];
    k_nn = 0;
    midPnt = thisEndPnt;
    midPntIdx = 8;
    minDist = 1;
    while minDist < 2
        k_nn = k_nn+1;
        distVec = cPnt-skelDistPnts;
        dists = sqrt(distVec(:,1).^2+distVec(:,2).^2);
        minDist = min(dists);
        pntIdx = find(dists==min(dists),1,'first');
        % determine new point
        cPnt = skelDistPnts(pntIdx,:);
        % take ot this point
        cond = ones(1,size(skelDistPnts,1),'logical');
        cond(pntIdx) = 0;
        skelDistPnts = skelDistPnts(cond,:);
        % save
        linePoints = [linePoints;cPnt];
    end
    lines(statsIdx).pnts = flip(linePoints,2);
end
for lineIdx = 1:length(lines)
    pnts = lines(lineIdx).pnts;
end
end
