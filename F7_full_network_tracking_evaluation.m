%--------------------------------------------------------------------------
% Script Name : F7_full_network_tracking_evaluation.m
% Author      : Felix Maurer
% Institution : Saarland University
% Email       : mail@felixmilanmaurer.com
% Date        : 2024
%
% Description :
%   Thi script evaluates the initial coarse cell tracking to obtain
%   velocity vectors and average magnitudes and directions.
%
% Usage :
%   - The script requires previous full network tracking
%
% Dependencies :
%
% Reference :
%   This script is associated with the publication
%   Impact of Red Blood Cell Rigidity on in vivo Flow Dynamics and Lingering in Bifurcations
%   by Rashidi et al. 2025
% License :
%   MIT
%%
clc;
cellTypes = {'Healthy_RBCs','Rigid_RBCs'};
cellTypeIdx = 1; % 1: Healthy_RBCs, 2:Rigid_RBCs
maskName = 'Mask.png';
rootDir = char(readlines('directory.txt'));
filelist = dir(fullfile(rootDir, '**\',maskName));  % get list of files and folders in any subfolder
rootFolders = unique({filelist.folder});
timer = tic;
suffix = '_traj';
for rootIdx = 1:length(rootFolders)
    rootdir = rootFolders{rootIdx};
    rootdir = rootdir(1:end-length('\converted_expo'));
    % display
    fprintf('working on -> %s\n',rootdir);
    % load geometry
    load([rootdir,'\converted_expo\geometry\bifurcations.mat']);

    % extract ROI boundaries
    wallImg = imread([rootdir '\converted_expo\',maskName]);
    if length(size(wallImg)) > 2
        wallImg = wallImg(:,:,1);
    end
    wallBW = imbinarize(wallImg);
    ROIbounds = bwboundaries(imfill(wallBW,'holes')); % one boundary for each disconnected branch

    % extract pixellist
    clear pixellist
    for roiIdx = 1:length(ROIbounds)
        roiBdy = ROIbounds{roiIdx};
        [x,y] = find(wallBW - imcomplement(imfill(imcomplement(wallBW),[roiBdy(1,1),roiBdy(1,2)])));
        pixellist(roiIdx).pnts = [x,y];
    end

    % plot boundaries and bif points
    figure;
    hold on
    for roiIdx = 1:length(ROIbounds)
        roiBdy = ROIbounds{roiIdx};
        plot(roiBdy(:,1),roiBdy(:,2),'-b')
    end
    for bifIdx = 1:length(bifurcations)
        bifCtr = bifurcations(bifIdx).bifCtr;
        plot(bifCtr(1),bifCtr(2),'.r')
    end
    hold off
    % load all trajectories
    filelist = dir(fullfile(rootdir, '**\ROI*network_merge.mat'));  % get list of files and folders in any subfolder
    filelist = filelist(~[filelist.isdir]);  % remove folders from list
    clear trajectories;
    if ~isempty(filelist)
        for fileIdx = 1:length(filelist)
            % file location and declaration
            fileFolder = filelist(fileIdx).folder;
            fileName = filelist(fileIdx).name;
            filePath = [fileFolder '\' fileName];
            % determine cell type
            slashIdc = strfind(fileFolder,'\');
            cellTypeStr = fileFolder(max(slashIdc)+1:end);
            if strcmp(cellTypeStr,'Healthy_RBCs')
                cellTypeIdx = 1;
            elseif strcmp(cellTypeStr,'Rigid_RBCs')
                cellTypeIdx = 2;
            end
            % determine roiIdx
            strIdx1 = strfind(fileName,'ROI_')+4;
            strIdx2 = strfind(fileName(strIdx1:end),'_');
            strIdx2 = strIdx2(1)+strIdx1-2;
            roiIdx = str2double(fileName(strIdx1:strIdx2));
            % load data
            load(filePath)
            % save
            trajectories(fileIdx).roiIdx = roiIdx;
            trajectories(fileIdx).cellType = cellTypeIdx;
            trajectories(fileIdx).clu = clu;
        end
        %% analysis
        close all
        clear velocity
        % take out ROIbounds without data
        selIdc = 1:length(ROIbounds);
        for roiIdx = 1:length(ROIbounds)
            selIdc(roiIdx) = sum([trajectories.roiIdx] == roiIdx)>0;
        end
        ROIbounds = {ROIbounds{find(selIdc)}};
        for roiIdx = 1:length(ROIbounds)
            traj = trajectories([trajectories.roiIdx] == roiIdx);
            traj = traj([traj.cellType] == 1);
            clu = traj.clu;
            roiBdy = ROIbounds{roiIdx};
            plot_trajectories = false;
            if plot_trajectories
                figure
                hold on
                plot(roiBdy(:,1),roiBdy(:,2),'-b')
                for cluIdx = 1:length(clu)
                    pnts = clu(cluIdx).points;
                    plot(pnts(:,1),pnts(:,2),'.-')
                end
                hold off
                axis equal
                axis square
            end
            % compute speed
            vs = [];
            ps = [];
            for cluIdx = 1:length(clu)
                pnts = clu(cluIdx).points;
                if ~isempty(pnts)
                    if size(pnts,1)>1
                        v = diff(pnts(:,1:2))./diff(pnts(:,3));
                        p = pnts(1:end-1,1:2);
                        vs = [vs;v];
                        ps = [ps;p];
                    end
                end
            end
            % add boundaries
            ps = [ps;roiBdy];
            vs = [vs;roiBdy*0];
            % query points
            xq = pixellist(roiIdx).pnts(:,1);
            yq = pixellist(roiIdx).pnts(:,2);

            x = ps(:,1);
            y = ps(:,2);
            vx = vs(:,1);
            vy = vs(:,2);

            [Xq, Yq] = meshgrid(min(x):max(x),min(y):max(y));
            % make a mask
            mask = double(wallBW(min(x):max(x),min(y):max(y)));
            mask(mask==0) = nan;
            warning ('off','all');
            % absolute velocity
            absV = sqrt(vx.^2+vy.^2);
            absVq = griddata(x, y, absV, Xq, Yq, 'cubic');
            absVq(isnan(absVq))=0;
            absVq = absVq.*double(wallBW(min(x):max(x),min(y):max(y))');
            absVq = imgaussfilt(absVq,3);
            absVq = absVq.*double(mask');
            % angle
            dirVec = vs;
            dirVec = dirVec./sqrt(dirVec(:,1).^2+dirVec(:,2).^2);
            circX = griddata(x, y, dirVec(:,1), Xq, Yq, 'cubic');
            circY = griddata(x, y, dirVec(:,2), Xq, Yq, 'cubic');
            % smoothing
            circX = smoothdata(circX,1,"gaussian",5);
            circX = smoothdata(circX,2,"gaussian",5);
            circY = smoothdata(circY,1,"gaussian",5);
            circY = smoothdata(circY,2,"gaussian",5);
            % compute angle
            angleq = -acos(circY./sqrt(circX.^2+circY.^2)).*(1-2*(-circX<0))/pi*180;
            angleq = angleq.*double(mask');
            warning ('on','all');
            % save everything
            velocity(roiIdx).X = Xq;
            velocity(roiIdx).Y = Yq;
            velocity(roiIdx).V = absVq;
            velocity(roiIdx).angle = angleq;
        end
        %% combine ROI
        tic
        Vs = zeros(size(wallBW),'double');
        angles = zeros(size(wallBW),'double');
        mask = ones(size(wallBW),'double');
        for roiIdx = 1:length(ROIbounds)
            X = velocity(roiIdx).X;
            Y = velocity(roiIdx).Y;
            V = velocity(roiIdx).V;
            angle = velocity(roiIdx).angle;
            for idx = 1:size(X,1)
                for idy = 1:size(X,2)
                    Vs(X(idx,idy),Y(idx,idy)) = Vs(X(idx,idy),Y(idx,idy))+V(idx,idy);
                    if angles(X(idx,idy),Y(idx,idy)) == 0 || isnan(angles(X(idx,idy),Y(idx,idy)))
                        angles(X(idx,idy),Y(idx,idy)) = angle(idx,idy);
                    end
                    mask(X(idx,idy),Y(idx,idy)) = 0;
                end
            end
        end
        % mask = mask & wallBW;
        angles(mask==1) = nan;
        Vs(mask==1) = nan;
        toc

        % save data
        save([rootdir,'\velocity_magnitude' suffix '.mat'],'Vs');
        save([rootdir,'\velocity_direction' suffix '.mat'],'angles');
        %% plotting
        close all;

        % velocity scale
        frameRate =  99.9255/2;
        dist = 1; % must be set to 1!!
        deltaT = 1/frameRate;
        microScale = 0.65; % micron/px: 20x is 0.65 µm/px, and 50x is 0.26 µm/px
        velScale = microScale * dist / deltaT;

        % velocity asbsolute
        Velocity = Vs*velScale;
        fig = figure;
        surf(Velocity,'EdgeColor','none')
        colormap(jet);
        view(2);
        cbar = colorbar;
        clim([prctile(Velocity(:),5),prctile(Velocity(:),95)]);
        set(cbar,'units','centimeters','position',[10*size(wallBW,2)/size(wallBW,1)+0.5 .5 .5 8])
        cbar.Label.String = 'speed in [microns/s]';
        cbar.Label.Rotation = 90; % Vertical text
        cbar.Label.Position(1) = cbar.Label.Position(1) - 4; % Adjust horizontal offset if needed
        fig.Units = 'centimeters';
        fig.Position = [2 5 10*size(wallBW,2)/size(wallBW,1)+2 10]; % fixed position
        ax = gca;
        ax.Units = 'centimeters';
        ax.Visible = 'off';
        ax.Position(1:2) = 0; % flush
        ax.Position(4) = fig.Position(4); % full coverage in height
        ax.Position(3) = ax.Position(4) * size(wallBW,2)/size(wallBW,1);

        % saving
        print([rootdir,'\converted_expo\','velocity_magnitude' suffix '.png'],'-dpng','-r600');

        % make angle colormap
        Nc = 6*32;
        eps = 1:Nc;
        Nint = round(length(eps)/6);
        % direction 1
        c = zeros(4,3);
        c(5,:) = [1 0 0];
        c(6,:) = [1 1 0];
        c(1,:) = [0 1 0];
        c(2,:) = [0 1 1];
        c(3,:) = [0 0 1];
        c(4,:) = [1 0 1];
        %c = 1-c;
        c = [c(end,:);c];
        idc = 0:length(eps);
        int = idc(mod(idc,(Nint))==0);
        angleCmap = zeros(length(eps),3,'double');
        for idx = 1:size(angleCmap,1)
            color = [0 0 0];
            for intIdx = 1:length(int)-1
                color = color + c(intIdx+1,:) * (idx-int(intIdx))/Nint * double(idx>int(intIdx)) * double(idx<=int(intIdx+1));
            end
            angleCmap(idx,:) = color;
        end
        a1 = angleCmap;
        % direction 2
        c = flip(c,1);
        idc = 0:length(eps);
        int = idc(mod(idc,(Nint))==0);
        angleCmap = zeros(length(eps),3,'double');
        for idx = 1:size(angleCmap,1)
            color = [0 0 0];
            for intIdx = 1:length(int)-1
                color = color + c(intIdx+1,:) * (idx-int(intIdx))/Nint * double(idx>int(intIdx)) * double(idx<=int(intIdx+1));
            end
            angleCmap(idx,:) = color;
        end
        a2 = angleCmap;
        angleCmap = flip(a2,1)+a1;
        % cap
        angleCmap(angleCmap>1) = 1;
        angleCmap(angleCmap<0) = 0;
        % angles
        fig = figure;
        surf(angles,'EdgeColor','none');
        colormap(angleCmap);
        view(2);
        cbar = colorbar('XTick', -180:45:180);
        set(cbar,'units','centimeters','position',[10*size(wallBW,2)/size(wallBW,1)+0.5 .5 .5 8])
        cbar.Label.String = 'angle to x-axis in [°]';
        cbar.Label.Rotation = 90; % Vertical text
        cbar.Label.Position(1) = cbar.Label.Position(1) - 4; % Adjust horizontal offset if needed
        clim([-180,180])
        fig.Units = 'centimeters';
        fig.Position = [25 5 10*size(wallBW,2)/size(wallBW,1)+2 10]; % fixed position
        ax = gca;
        ax.Units = 'centimeters';
        ax.Visible = 'off';
        ax.Position(1:2) = 0; % flush
        ax.Position(4) = fig.Position(4); % full coverage in height
        ax.Position(3) = ax.Position(4) * size(wallBW,2)/size(wallBW,1);

        % saving
        print([rootdir,'\converted_expo\','velocity_direction' suffix '.png'],'-dpng','-r600');

        % status
        time = toc(timer);
        fprintf('status: %.2f percent.\n',rootIdx/length(rootFolders)*100);
        fprintf('remaining time: %.2f sec.\n',time/rootIdx * (1-rootIdx/length(rootFolders)));
    end
end
return
%% debugging: plot boundaries
figure
hold on
for idx = 1:length(ROIbounds)
    pnts = ROIbounds{idx};
    plot(pnts(:,1),pnts(:,2))
end
hold off
