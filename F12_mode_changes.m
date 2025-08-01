%--------------------------------------------------------------------------
% Script Name : F11_lingering_bar_plots.m
% Authors     : Felix Maurer, Yazdan Rashidi
% Institution : Saarland University
% Email       : mail@felixmilanmaurer.com
% Date        : 2024
%
% Description :
%   This script analysis modes, specifically wall contact mode leading to
%   the fraction of interacting cells. 
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
rootDir = char(readlines('directory.txt'));
filelist = dir(fullfile([rootDir,'\**\','bifTraj.mat']));
folders = {filelist.folder};
folders = unique(folders);
%%
HRBCthres = 4;
RRBCthres = 4;
zPositions = 0;
for rootIdx = 1:length(folders)
    cfolder = folders{rootIdx};
    load([cfolder,'\bifTraj.mat']);
    for bifIdx = 1:length(bifurcations)
        cBif = bifurcations(bifIdx);
        conVes = cBif.conVesBdy;
        HRBCtraj = cBif.HRBCtraj;
        RRBCtraj = cBif.RRBCtraj;
        bifBW = bifurcations(bifIdx).bifBW;
        for vesIdx = 1:length(conVes)
            cVes = conVes(vesIdx);
            bifBW = bifBW | cVes.vesBW;
        end
        %% find mother branch
        bifAngle = cBif.bifAvgAngle;
        cVes = cBif.conVesBdy;
        angleDif = zeros(1,length(cVes));
        for vesIdx = 1:length(cVes)
            vesStartPnts = flip(cVes(vesIdx).pnts(1:5,:),1);
            dirVec = mean(diff(vesStartPnts));
            vesAngle = -acos(dirVec(:,2)).*(1-2*(-dirVec(:,1)<0))/pi*180;
            angleDif(vesIdx) = abs(bifAngle-vesAngle);
        end
        inVesIdx = find(angleDif == min(angleDif),1,'first');
        %% daughter bw
        bifDaugh = bifurcations(bifIdx).bifBW;
        for vesIdx = 1:length(conVes)
            cVes = conVes(vesIdx);
            if vesIdx ~= inVesIdx
                bifDaugh = bifDaugh | cVes.vesBW;
            end
        end
        bifDaughBounds = bwboundaries(bifDaugh);
        bifDaughBounds = bifDaughBounds{1};
        %% HRBCs
        intThres = mean([cBif.conVesBdy.HRBCVel],'omitnan')./cBif(1).velScaleHRBC*0.3;
        figure;
        imshow(bifBW);
        hold on
        k_Traj = 0;
        chIdc = [];
        HRBCVel = [];
        HRBCModes = [];
        for HRBCidx = 1:length(HRBCtraj)
            pnts = HRBCtraj(HRBCidx).pnts;
            inboundIdc = inpolygon(pnts(:,1),pnts(:,2),bifDaughBounds(:,1),bifDaughBounds(:,2));
            if sum(inboundIdc) > 8
                pnts = pnts(inboundIdc,:);
                velVec = [diff(pnts(:,1)),diff(pnts(:,2))]./diff(pnts(:,3));
                velMag = sqrt(velVec(:,1).^2+velVec(:,2).^2);
                vel = velMag;
                medVel = median(vel(~isnan(vel)));
                vel(isnan(vel)) = medVel;
                HRBCVel = [HRBCVel,medVel];
                intCond = vel < intThres;
                intCond = logical(hampel(double(intCond),3));
                intCond = logical(round(movmean(double(intCond),4)));
                trajMode = 0;
                if sum(intCond) > floor(length(intCond)/3)
                    trajMode = 1;
                end
                HRBCModes = [HRBCModes,trajMode];
                modeCh = diff(intCond);
                modeChIdc = find(modeCh ~= 0);
                if ~isempty(modeChIdc)
                    chIdc = [chIdc,HRBCidx];
                    k_Traj = k_Traj + 1;
                    if ismember(k_Traj,[1 5 6 10 11])
                        offset = -min(pnts(:,3))+sum(zPositions);
                        zPositions = [zPositions,max(pnts(:,3))-min(pnts(:,3))];
                        firstMode = double(intCond(1));
                        numFlo = sum(modeCh < 0) + double(firstMode==0);
                        numRol = sum(modeCh > 0) + double(firstMode==1);

                        floInts = zeros(numFlo,2);
                        interactionInts = zeros(numRol,2);
                        modeChIdc = [0;modeChIdc;length(modeCh)+1];

                        k_int = 1;
                        k_flo = 1;
                        for chIdx = 1:length(modeChIdc)-1
                            if intCond(modeChIdc(chIdx)+1)
                                interactionInts(k_int,:) = [modeChIdc(chIdx)+1,modeChIdc(chIdx+1)];
                                k_int = k_int + 1;
                            else
                                floInts(k_flo,:) = [modeChIdc(chIdx)+1,modeChIdc(chIdx+1)];
                                k_flo = k_flo + 1;
                            end
                        end
                        for intIdx = 1:size(interactionInts,1)
                            intPnts = pnts(interactionInts(intIdx,1):interactionInts(intIdx,2),:);
                            plot3(intPnts(:,2),intPnts(:,1),intPnts(:,3)+offset,'-','color',[0.2 0.6 0.4],'LineWidth',1.3);
                        end
                        for intIdx = 1:size(floInts,1)
                            intPnts = pnts(floInts(intIdx,1):floInts(intIdx,2),:);
                            plot3(intPnts(:,2),intPnts(:,1),intPnts(:,3)+offset,'-','color',[0.2 0.3 0.7],'LineWidth',1.3);
                        end
                        % transition points
                        intChIdc = find(modeCh > 0);
                        floChIdc = find(modeCh < 0);
                        for chIdx = 1:size(intChIdc,1)
                            chPnts = pnts(intChIdc(chIdx):intChIdc(chIdx)+1,:);
                            plot3(chPnts(:,2),chPnts(:,1),chPnts(:,3)+offset,'o-','color',[0.7 0.3 0.3],'MarkerFaceColor',[0.7 0.3 0.3]);
                        end
                        for chIdx = 1:size(floChIdc,1)
                            chPnts = pnts(floChIdc(chIdx):floChIdc(chIdx)+1,:);
                            plot3(chPnts(:,2),chPnts(:,1),chPnts(:,3)+offset,'s-','color',[0.7 0.5 0.9],'MarkerFaceColor',[0.7 0.5 0.9]);
                        end
                    end
                else
                   
                end
            end
        end
        hold off
        title('HRBC')
        bifurcations(bifIdx).HRBCMode = HRBCModes;
        bifurcations(bifIdx).HRBCModeThres = intThres;

        %% RRBCs
        intThres = mean([cBif.conVesBdy.RRBCVel],'omitnan')./cBif(1).velScaleRRBC*0.3;

        hold on
        k_Traj = 0;
        chIdc = [];
        RRBCVel = [];
        RRBCModes = [];
        for RRBCidx = 1:length(RRBCtraj)
            pnts = RRBCtraj(RRBCidx).pnts;
            inboundIdc = inpolygon(pnts(:,1),pnts(:,2),bifDaughBounds(:,1),bifDaughBounds(:,2));
            if sum(inboundIdc) > 8
                pnts = pnts(inboundIdc,:);
                velVec = [diff(pnts(:,1)),diff(pnts(:,2))]./diff(pnts(:,3));
                velMag = sqrt(velVec(:,1).^2+velVec(:,2).^2);
                vel = velMag;
                medVel = median(vel(~isnan(vel)));
                vel(isnan(vel)) = medVel;
                RRBCVel = [RRBCVel,medVel];
                intCond = vel < intThres;
                intCond = logical(hampel(double(intCond),3));
                intCond = logical(round(movmean(double(intCond),4)));
                %%
                trajMode = 0;
                if sum(intCond) > floor(length(intCond)/3)
                    trajMode = 1;
                end
                RRBCModes = [RRBCModes,trajMode];
                %%
                modeCh = diff(intCond);
                modeChIdc = find(modeCh ~= 0);
                if ~isempty(modeChIdc)
                    chIdc = [chIdc,RRBCidx];
                    k_Traj = k_Traj + 1;
                    if ismember(k_Traj,[1 5 6 10 11])
                        offset = -min(pnts(:,3))+sum(zPositions);
                        zPositions = [zPositions,max(pnts(:,3))-min(pnts(:,3))];
                        firstMode = double(intCond(1));
                        numFlo = sum(modeCh < 0) + double(firstMode==0);
                        numRol = sum(modeCh > 0) + double(firstMode==1);

                        floInts = zeros(numFlo,2);
                        interactionInts = zeros(numRol,2);
                        modeChIdc = [0;modeChIdc;length(modeCh)+1];

                        k_int = 1;
                        k_flo = 1;
                        for chIdx = 1:length(modeChIdc)-1
                            if intCond(modeChIdc(chIdx)+1)
                                interactionInts(k_int,:) = [modeChIdc(chIdx)+1,modeChIdc(chIdx+1)];
                                k_int = k_int + 1;
                            else
                                floInts(k_flo,:) = [modeChIdc(chIdx)+1,modeChIdc(chIdx+1)];
                                k_flo = k_flo + 1;
                            end
                        end
                        for intIdx = 1:size(interactionInts,1)
                            intPnts = pnts(interactionInts(intIdx,1):interactionInts(intIdx,2),:);
                        end
                        for intIdx = 1:size(floInts,1)
                            intPnts = pnts(floInts(intIdx,1):floInts(intIdx,2),:);
                        end
                        % transition points
                        intChIdc = find(modeCh > 0);
                        floChIdc = find(modeCh < 0);
                        for chIdx = 1:size(intChIdc,1)
                            chPnts = pnts(intChIdc(chIdx):intChIdc(chIdx)+1,:);
                        end
                        for chIdx = 1:size(floChIdc,1)
                            chPnts = pnts(floChIdc(chIdx):floChIdc(chIdx)+1,:);
                        end
                    end
                else
                    % %plot3(pnts(:,2),pnts(:,1),pnts(:,3),'color',[0.8 0.8 0.8]);
                end
            end
        end
        hold off
        title('RRBC')
        bifurcations(bifIdx).RRBCMode = RRBCModes;
        bifurcations(bifIdx).RRBCModeThres = intThres;

        %% more parameters
        RRBCModes = logical(RRBCModes);
        HRBCModes = logical(HRBCModes);

        RRBCRolVel = RRBCVel(RRBCModes);
        RRBCFloVel = RRBCVel(~RRBCModes);

        HRBCRolVel = HRBCVel(HRBCModes);
        HRBCFloVel = HRBCVel(~HRBCModes);

        bifurcations(bifIdx).RRBCFloVel = RRBCFloVel;
        bifurcations(bifIdx).RRBCRolVel = RRBCRolVel;
        bifurcations(bifIdx).HRBCFloVel = HRBCFloVel;
        bifurcations(bifIdx).HRBCRolVel = HRBCRolVel;

        fprintf('----------------------------------------------------\n');
        fprintf('%.2f percent healthy int.\n',sum(HRBCModes)/length(HRBCModes)*100);
        fprintf('%.2f percent rigid int.\n',sum(RRBCModes)/length(RRBCModes)*100);
        fprintf('%.3f healthy flow vel.\n',mean(HRBCFloVel(~isnan(HRBCFloVel)&~isinf(HRBCFloVel))));
        fprintf('%.3f rigid flow vel.\n',mean(RRBCFloVel(~isnan(RRBCFloVel)&~isinf(RRBCFloVel))));
        fprintf('%.3f healthy flow vel / healthy vel\n',mean(HRBCFloVel(~isnan(HRBCFloVel)&~isinf(HRBCFloVel)))/mean(HRBCVel(~isnan(HRBCVel)&~isinf(HRBCVel))));
        fprintf('%.3f healthy flow vel / rigid flow vel\n',mean(HRBCFloVel(~isnan(HRBCFloVel)&~isinf(HRBCFloVel)))/mean(RRBCFloVel(~isnan(RRBCFloVel)&~isinf(RRBCFloVel))));
        fprintf('----------------------------------------------------\n');
        %% vessel-wise
        if 1
            for vesIdx = 1:length(conVes)
                vesBW = conVes(vesIdx).vesBW;
                vesBounds = bwboundaries(vesBW);
                vesBounds = vesBounds{1};
                %% HRBCs
                intThres = conVes(vesIdx).HRBCVel/cBif.velScaleHRBC*0.3;
                %figure;
                %imshow(bifBW);
                hold on
                k_Traj = 0;
                chIdc = [];
                HRBCVel = [];
                HRBCModes = [];
                for HRBCidx = 1:length(HRBCtraj)
                    pnts = HRBCtraj(HRBCidx).pnts;
                    inboundIdc = inpolygon(pnts(:,1),pnts(:,2),vesBounds(:,1),vesBounds(:,2));
                    if sum(inboundIdc) > 3
                        pnts = pnts(inboundIdc,:);
                        velVec = [diff(pnts(:,1)),diff(pnts(:,2))]./diff(pnts(:,3));
                        velMag = sqrt(velVec(:,1).^2+velVec(:,2).^2);
                        vel = velMag;
                        medVel = median(vel(~isnan(vel)));
                        vel(isnan(vel)) = medVel;
                        HRBCVel = [HRBCVel,medVel];
                        intCond = vel < intThres;
                        intCond = logical(hampel(double(intCond),3));
                        intCond = logical(round(movmean(double(intCond),4)));
                        trajMode = 0;
                        if sum(intCond) > floor(length(intCond)/2)
                            trajMode = 1;
                        end
                        HRBCModes = [HRBCModes,trajMode];
                        modeCh = diff(intCond);
                        modeChIdc = find(modeCh ~= 0);
                        if ~isempty(modeChIdc)
                            chIdc = [chIdc,HRBCidx];
                            k_Traj = k_Traj + 1;
                            if ismember(k_Traj,[1 5 6 10 11])
                                offset = -min(pnts(:,3))+sum(zPositions);
                                zPositions = [zPositions,max(pnts(:,3))-min(pnts(:,3))];
                                firstMode = double(intCond(1));
                                numFlo = sum(modeCh < 0) + double(firstMode==0);
                                numRol = sum(modeCh > 0) + double(firstMode==1);

                                floInts = zeros(numFlo,2);
                                interactionInts = zeros(numRol,2);
                                modeChIdc = [0;modeChIdc;length(modeCh)+1];

                                k_int = 1;
                                k_flo = 1;
                                for chIdx = 1:length(modeChIdc)-1
                                    if intCond(modeChIdc(chIdx)+1)
                                        interactionInts(k_int,:) = [modeChIdc(chIdx)+1,modeChIdc(chIdx+1)];
                                        k_int = k_int + 1;
                                    else
                                        floInts(k_flo,:) = [modeChIdc(chIdx)+1,modeChIdc(chIdx+1)];
                                        k_flo = k_flo + 1;
                                    end
                                end
                                for intIdx = 1:size(interactionInts,1)
                                    intPnts = pnts(interactionInts(intIdx,1):interactionInts(intIdx,2),:);
                                    %plot3(intPnts(:,2),intPnts(:,1),intPnts(:,3)+offset,'-','color',[0.2 0.6 0.4],'LineWidth',1.3);
                                end
                                for intIdx = 1:size(floInts,1)
                                    intPnts = pnts(floInts(intIdx,1):floInts(intIdx,2),:);
                                    %plot3(intPnts(:,2),intPnts(:,1),intPnts(:,3)+offset,'-','color',[0.2 0.3 0.7],'LineWidth',1.3);
                                end
                                % transition points
                                intChIdc = find(modeCh > 0);
                                floChIdc = find(modeCh < 0);
                                for chIdx = 1:size(intChIdc,1)
                                    chPnts = pnts(intChIdc(chIdx):intChIdc(chIdx)+1,:);
                                    %plot3(chPnts(:,2),chPnts(:,1),chPnts(:,3)+offset,'o-','color',[0.7 0.3 0.3],'MarkerFaceColor',[0.7 0.3 0.3]);
                                end
                                for chIdx = 1:size(floChIdc,1)
                                    chPnts = pnts(floChIdc(chIdx):floChIdc(chIdx)+1,:);
                                    %plot3(chPnts(:,2),chPnts(:,1),chPnts(:,3)+offset,'s-','color',[0.7 0.5 0.9],'MarkerFaceColor',[0.7 0.5 0.9]);
                                end
                            end
                        else
                            % %plot3(pnts(:,2),pnts(:,1),pnts(:,3),'color',[0.8 0.8 0.8]);
                        end
                    end
                end
                hold off
                title('HRBC')
                conVes(vesIdx).HRBCMode = HRBCModes;
                conVes(vesIdx).HRBCModeThres = intThres;
                conVes(vesIdx).HRBCVel = HRBCVel;
                %% RRBCs
                intThres = conVes(vesIdx).RRBCVel/cBif.velScaleRRBC*0.3;
                hold on
                k_Traj = 0;
                chIdc = [];
                RRBCVel = [];
                RRBCModes = [];
                for RRBCidx = 1:length(RRBCtraj)
                    pnts = RRBCtraj(RRBCidx).pnts;
                    inboundIdc = inpolygon(pnts(:,1),pnts(:,2),vesBounds(:,1),vesBounds(:,2));
                    if sum(inboundIdc) > 3
                        pnts = pnts(inboundIdc,:);
                        velVec = [diff(pnts(:,1)),diff(pnts(:,2))]./diff(pnts(:,3));
                        velMag = sqrt(velVec(:,1).^2+velVec(:,2).^2);
                        vel = velMag;
                        medVel = median(vel(~isnan(vel)));
                        vel(isnan(vel)) = medVel;
                        RRBCVel = [RRBCVel,medVel];
                        intCond = vel < intThres;
                        intCond = logical(hampel(double(intCond),3));
                        intCond = logical(round(movmean(double(intCond),4)));
                        %%
                        trajMode = 0;
                        if sum(intCond) > floor(length(intCond)/2)
                            trajMode = 1;
                        end
                        RRBCModes = [RRBCModes,trajMode];
                        %%
                        modeCh = diff(intCond);
                        modeChIdc = find(modeCh ~= 0);
                        if ~isempty(modeChIdc)
                            chIdc = [chIdc,RRBCidx];
                            k_Traj = k_Traj + 1;
                            if ismember(k_Traj,[1 5 6 10 11])
                                offset = -min(pnts(:,3))+sum(zPositions);
                                zPositions = [zPositions,max(pnts(:,3))-min(pnts(:,3))];
                                firstMode = double(intCond(1));
                                numFlo = sum(modeCh < 0) + double(firstMode==0);
                                numRol = sum(modeCh > 0) + double(firstMode==1);

                                floInts = zeros(numFlo,2);
                                interactionInts = zeros(numRol,2);
                                modeChIdc = [0;modeChIdc;length(modeCh)+1];

                                k_int = 1;
                                k_flo = 1;
                                for chIdx = 1:length(modeChIdc)-1
                                    if intCond(modeChIdc(chIdx)+1)
                                        interactionInts(k_int,:) = [modeChIdc(chIdx)+1,modeChIdc(chIdx+1)];
                                        k_int = k_int + 1;
                                    else
                                        floInts(k_flo,:) = [modeChIdc(chIdx)+1,modeChIdc(chIdx+1)];
                                        k_flo = k_flo + 1;
                                    end
                                end
                                for intIdx = 1:size(interactionInts,1)
                                    intPnts = pnts(interactionInts(intIdx,1):interactionInts(intIdx,2),:);
                                    %plot3(intPnts(:,2),intPnts(:,1),intPnts(:,3)+offset,'-','color',[0.2 0.6 0.4],'LineWidth',1.3);
                                end
                                for intIdx = 1:size(floInts,1)
                                    intPnts = pnts(floInts(intIdx,1):floInts(intIdx,2),:);
                                    %plot3(intPnts(:,2),intPnts(:,1),intPnts(:,3)+offset,'-','color',[0.2 0.3 0.7],'LineWidth',1.3);
                                end
                                % transition points
                                intChIdc = find(modeCh > 0);
                                floChIdc = find(modeCh < 0);
                                for chIdx = 1:size(intChIdc,1)
                                    chPnts = pnts(intChIdc(chIdx):intChIdc(chIdx)+1,:);
                                    %plot3(chPnts(:,2),chPnts(:,1),chPnts(:,3)+offset,'o-','color',[0.7 0.3 0.3],'MarkerFaceColor',[0.7 0.3 0.3]);
                                end
                                for chIdx = 1:size(floChIdc,1)
                                    chPnts = pnts(floChIdc(chIdx):floChIdc(chIdx)+1,:);
                                    %plot3(chPnts(:,2),chPnts(:,1),chPnts(:,3)+offset,'s-','color',[0.7 0.5 0.9],'MarkerFaceColor',[0.7 0.5 0.9]);
                                end
                            end
                        else
                            % %plot3(pnts(:,2),pnts(:,1),pnts(:,3),'color',[0.8 0.8 0.8]);
                        end
                    end
                end
                hold off
                title('RRBC')
                conVes(vesIdx).RRBCMode = RRBCModes;
                conVes(vesIdx).RRBCModeThres = intThres;
                conVes(vesIdx).RRBCVel = RRBCVel;
            end
        end
        bifurcations(bifIdx).conVesBdy = conVes;
        bifurcations(bifIdx).orVes = inVesIdx;
    end
    %% save
    save([cfolder,'\bifTrajModes.mat']);
    close all;
end