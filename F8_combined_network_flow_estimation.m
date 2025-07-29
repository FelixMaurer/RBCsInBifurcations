%--------------------------------------------------------------------------
% Script Name : F8_combined_network_flow_estimation.m
% Authors     : Felix Maurer, Yazdan Rashidi
% Institution : Saarland University
% Email       : mail@felixmilanmaurer.com
% Date        : 2024
%
% Description :
%   This script combines flow estimate results from the Lucan Kanade method
%   with results from coarse cell tracking to obtain a more robust average
%   flow speed and angle for a predictive search in the refined tracking
%   algorithm
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
%% file location and loop
maskName = 'Mask.png';
rootDir = char(readlines('directory.txt'));
filelist = dir(fullfile(rootDir, '**\',maskName));  % get list of files and folders in any subfolder
rootFolders = unique({filelist.folder});
timer = tic;
suffix = '_traj';
for rootIdx = 1:length(rootFolders)
    rootdir = rootFolders{rootIdx};
    % load wall image
    %try
    wallImg = imread([rootdir,'\',maskName]);
    if length(size(wallImg)) > 2
        wallImg = wallImg(:,:,1);
    end
    wallBW = imbinarize(wallImg);
    % search files
    bif_file = dir(fullfile(rootDir, '**', 'geometry', 'bifurcations.mat'));
    if ~isempty(bif_file)
        load(fullfile(bif_file(1).folder, bif_file(1).name));
    else
        error('bifurcations.mat not found');
    end
    
    % Search and load velocity_magnitude_traj.mat
    vel_file1 = dir(fullfile(rootDir, '**', '*velocity_magnitude_traj.mat'));
    if ~isempty(vel_file1)
        load(fullfile(vel_file1(1).folder, vel_file1(1).name));
        Vs1 = Vs;
    else
        error('velocity_magnitude_traj.mat not found');
    end
    
    % Search and load velocity_magnitude_LK_scale.mat
    vel_file2 = dir(fullfile(rootDir, '**', '*velocity_magnitude_LK_scale.mat'));
    if ~isempty(vel_file2)
        load(fullfile(vel_file2(1).folder, vel_file2(1).name));
        Vs2 = Vs;
    else
        error('velocity_magnitude_LK_scale.mat not found');
    end

    % Search and load velocity_direction_traj.mat
    agl_file1 = dir(fullfile(rootDir, '**', '*velocity_direction_traj.mat'));
    if ~isempty(agl_file1)
        load(fullfile(agl_file1(1).folder, agl_file1(1).name));
        angle1 = angles;
    else
        error('velocity_direction_traj.mat not found');
    end
    
    % Search and load velocity_direction_LK_scale.mat
    agl_file2 = dir(fullfile(rootDir, '**', '*velocity_direction_LK_scale.mat'));
    if ~isempty(agl_file2)
        load(fullfile(agl_file2(1).folder, agl_file2(1).name));
        angle2 = angles;
    else
        error('velocity_direction_LK_scale.mat not found');
    end

    Vs = zeros(size(Vs1),'double');
    Vs1(Vs1==0) = nan;
    Vs2(Vs2==0) = nan;
    for idx = 1:size(Vs1,1)
        for idy = 1:size(Vs1,2)
            k = 0;
            if ~isnan(Vs1(idx,idy))
                Vs(idx,idy) = Vs(idx,idy)+Vs1(idx,idy);
                k = k+1;
            end
            if ~isnan(Vs2(idx,idy))
                Vs(idx,idy) = Vs(idx,idy)+Vs2(idx,idy);
                k = k+1;
            end
            Vs(idx,idy) = Vs(idx,idy)/k;
        end
    end
    xCoor = zeros(size(angle1),'double');
    yCoor = zeros(size(angle1),'double');
    angles = zeros(size(angle1),'double');
    for idx = 1:size(angle1,1)
        for idy = 1:size(angle1,2)
            k = 0;
            if ~isnan(angle1(idx,idy))
                %angles(idx,idy) = angles(idx,idy)+angle1(idx,idy);
                xCoor(idx,idy) = xCoor(idx,idy)+cos(angle1(idx,idy)/180*pi);
                yCoor(idx,idy) = yCoor(idx,idy)+sin(angle1(idx,idy)/180*pi);
                k = k+1;
            end
            if ~isnan(angle2(idx,idy))
                %angles(idx,idy) = angles(idx,idy)+angle2(idx,idy);#
                xCoor(idx,idy) = xCoor(idx,idy)+cos(angle2(idx,idy)/180*pi);
                yCoor(idx,idy) = yCoor(idx,idy)+sin(angle2(idx,idy)/180*pi);
                k = k+1;
            end
            %angles(idx,idy) = angles(idx,idy)/k;
            xCoor(idx,idy) = xCoor(idx,idy)/k;
            yCoor(idx,idy) = yCoor(idx,idy)/k;
        end
    end
    angles = -acos(yCoor./sqrt(xCoor.^2+yCoor.^2)).*(1-2*(-xCoor<0))/pi*180;
    %% go through all bifurcations
    close all;
    % create bifurcation binary image
    allBifBW = zeros(size(wallBW),'logical');
    for bifIdx = 1:length(bifurcations)
        inPnts = bifurcations(bifIdx).bifInnerPnts;
        x = inPnts(:,1);
        y = inPnts(:,2);
        bifBW = zeros(size(wallBW),'logical');
        for pntIdx = 1:size(inPnts,1)
            if inbounds(bifBW,x(pntIdx),y(pntIdx))
                bifBW(x(pntIdx),y(pntIdx)) = 1;
            end
        end
        bifBW = bwmorph(bifBW,'thicken',3);
        bifurcations(bifIdx).bifBW = bifBW;
        allBifBW = allBifBW | bifBW;
    end

    colors = jet(256);
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

    VsMap = zeros(size(wallBW,1),size(wallBW,2),'double');
    AngleMap = zeros(size(wallBW,1),size(wallBW,2),'double');

    for bifIdx = 1:length(bifurcations)
        bifBW = bifurcations(bifIdx).bifBW;
        conVes = bifurcations(bifIdx).conVesBdy;
        for vesIdx = 1:length(conVes)
            pnts = conVes(vesIdx).pnts;
            normals = conVes(vesIdx).normals;
            unwrap = conVes(vesIdx).unwrap;
            bifPos = conVes(vesIdx).bifEndPos;
            bifPos = max(bifPos(:,1));
            vesBW = zeros(size(wallBW),'logical');
            [y,x] = find(unwrap);
            for pntIdx = 1:length(x)
                if x(pntIdx) > bifPos
                    imgPnt = round(pnts(x(pntIdx),:) + normals(x(pntIdx),:) * (-size(unwrap,1)/2+y(pntIdx)));
                    if inbounds(vesBW,imgPnt(1),imgPnt(2))
                        vesBW(imgPnt(1),imgPnt(2)) = 1;
                    end
                end
            end
            % fill
            vesBW = bwmorph(vesBW,'close');
            vesBW = (vesBW - allBifBW)>0;
            stats = regionprops(vesBW,'Area','PixelList');
            [~,sortIdc] = sort([stats.Area]);
            sortIdc = flip(sortIdc);
            if length(sortIdc) > 1 % take out smaller parts
                stats = stats(sortIdc(2:end));
                for idx = 1:length(stats)
                    pnts = flip(stats(idx).PixelList,2);
                    for fillIdx = 1:5:size(pnts,1)
                        vesBW = imcomplement(imfill(imcomplement(vesBW),[pnts(fillIdx,:)]));
                    end
                end
            end
            vesBW = bwmorph(vesBW,'clean');
            % sample all flow values
            [x,y] = find(vesBW);
            flow = [];
            angle = [];
            for pntIdx = 1:length(x)
                if inbounds(vesBW,x(pntIdx),y(pntIdx))
                    flow = [flow,Vs(x(pntIdx),y(pntIdx))];
                    angle = [angle,angles(x(pntIdx),y(pntIdx))];
                end
            end
            flowSelIdc = flow>prctile(flow,10);
            flow = flow(flowSelIdc);
            angle = angle(flowSelIdc);
            avgFlow = median(flow(~isnan(flow)));
            angleVals = angle(~isnan(angle) & ~isinf(angle));
            circlePnts = median([cos(angleVals/180*pi);sin(angleVals/180*pi)]');
            avgAngle = -acos(circlePnts(2)./sqrt(circlePnts(1).^2+circlePnts(2).^2)).*(1-2*(-circlePnts(1)<0))/pi*180;
            % draw on image
            for pntIdx = 1:length(x)
                VsMap(x(pntIdx),y(pntIdx),:) = avgFlow;
            end
            % draw on image
            for pntIdx = 1:length(x)
                AngleMap(x(pntIdx),y(pntIdx),:) = avgAngle;
            end
            %plot(inPnts(:,2),inPnts(:,1))
            % save variable
            conVes(vesIdx).avgFlow = avgFlow;
            conVes(vesIdx).avgAngle = avgAngle;
            conVes(vesIdx).vesBW = vesBW;
        end
        % save
        bifurcations(bifIdx).conVesBdy = conVes;
    end
    for bifIdx = 1:length(bifurcations)
        bifBW = bifurcations(bifIdx).bifBW;
        bifPnts = bifurcations(bifIdx).bifInnerPnts;
        % find values in bifurcation
        [x,y] = find(bifBW);
        flow = [];
        angle = [];
        for pntIdx = 1:size(x,1)
            if inbounds(bifBW,x(pntIdx),y(pntIdx))
                flow = [flow,Vs(x(pntIdx),y(pntIdx))];
                angle = [angle,angles(x(pntIdx),y(pntIdx))];
            end
        end
        flowSelIdc = flow>prctile(flow,10);
        flow = flow(flowSelIdc);
        angle = angle(flowSelIdc);
        avgFlow = median(flow(~isnan(flow)));
        angleVals = angle(~isnan(angle) & ~isinf(angle));
        circlePnts = median([cos(angleVals/180*pi);sin(angleVals/180*pi)]');
        avgAngle = -acos(circlePnts(2)./sqrt(circlePnts(1).^2+circlePnts(2).^2)).*(1-2*(-circlePnts(1)<0))/pi*180;
        % draw on image
        for pntIdx = 1:length(x)
            VsMap(x(pntIdx),y(pntIdx),:) = avgFlow;
        end
        for pntIdx = 1:length(x)
            AngleMap(x(pntIdx),y(pntIdx),:) = avgAngle;
        end
        bifurcations(bifIdx).bifAvgFlow = avgFlow;
        bifurcations(bifIdx).bifAvgAngle = avgAngle;
    end
    %%

    % saving
    AngleMap(VsMap==0) = nan;
    VsMap(VsMap==0) = nan;

    save([rootdir '\bifFlow.mat'],'bifurcations');
    save([rootdir '\FlowMap.mat'],'VsMap');
    save([rootdir '\AngleMap.mat'],'AngleMap');

    % plotting

    % velocity scale
    frameRate =  99.9255/2;
    dist = 1; % always set to 1
    deltaT = 1/frameRate;
    microScale = 0.65; % micron/px: 20x is 0.65 µm/px, and 50x is 0.26 µm/px
    velScale = microScale * dist/deltaT;
    VsMap = VsMap * velScale;
    Vs = Vs * velScale;

    close all
    % speed
    fig = figure;
    surf(VsMap,'EdgeColor','none');view(2);
    colormap(colors);
    cbar = colorbar;
    set(cbar,'units','centimeters','position',[10*size(wallBW,2)/size(wallBW,1)+0.5 .5 .5 8])
    cbar.Label.String = 'speed in [microns/s]';
    cbar.Label.Rotation = 90; % Vertical text
    cbar.Label.Position(1) = cbar.Label.Position(1) - 4; % Adjust horizontal offset if needed
    clim([prctile(Vs(:),10),prctile(Vs(:),90)])

    fig.Units = 'centimeters';
    fig.Position = [2 5 10*size(wallBW,2)/size(wallBW,1)+2 10]; % fixed position
    ax = gca;
    ax.Units = 'centimeters';
    ax.Visible = 'off';
    ax.Position(1:2) = 0; % flush
    ax.Position(4) = fig.Position(4); % full coverage in height
    ax.Position(3) = ax.Position(4) * size(wallBW,2)/size(wallBW,1);

    print([rootdir,'\RGBVs.png'],'-dpng','-r600');

    %% angle
    fig = figure;
    
    surf(AngleMap,'EdgeColor','none');view(2);
    
    colormap(angleCmap)
    cbar = colorbar('XTick', -180:45:180);
    set(cbar,'units','centimeters','position',[10*size(wallBW,2)/size(wallBW,1)+0.5 .5 .5 8])
    cbar.Label.String = 'angle to x-axis [°]';
    cbar.Label.Rotation = 90; % Vertical text
    cbar.Label.Position(1) = cbar.Label.Position(1) - 4; % Adjust horizontal offset if needed
    clim([-180,180])

    % scaling
    fig.Units = 'centimeters';
    fig.Position = [2 5 10*size(wallBW,2)/size(wallBW,1)+2 10]; % fixed position
    ax = gca;
    ax.Units = 'centimeters';
    ax.Visible = 'off';
    ax.Position(1:2) = 0; % flush
    ax.Position(4) = fig.Position(4); % full coverage in height
    ax.Position(3) = ax.Position(4) * size(wallBW,2)/size(wallBW,1);
    ax.Units = 'normalized';
    % add arrows (after scaling)
    drawArrow = @(x,y) quiver3( x(1),y(1),10000,x(2)-x(1),y(2)-y(1),0,0 ,'Color',[0 0 0],'LineWidth',1,'AutoScale','off');    
    hold on
    for bifIdx = 1:length(bifurcations)
        cPnt = bifurcations(bifIdx).bifCtr;
        bifAngle = bifurcations(bifIdx).bifAvgAngle;
        % add bifurcation
        cAngle = bifAngle;
        dirVec = [sin(cAngle/180*pi),cos(cAngle/180*pi)];
        dPnt = cPnt + dirVec*40;
        drawArrow([cPnt(2),dPnt(2)],[cPnt(1),dPnt(1)]);
        %annotation('textarrow',[cPnt(2),dPnt(2)],[cPnt(1),dPnt(1)],'String','y = x ')
        pos = get(gca, 'Position');
        xArrow = [cPnt(2),dPnt(2)];
        yArrow = [cPnt(1),dPnt(1)];
        ta1 = annotation('textarrow',...
        (xArrow-min(xlim))/diff(xlim)*pos(3)+pos(1),(yArrow-min(ylim))/diff(ylim)*pos(4)+pos(2));
        % add vessels
        conVes = bifurcations(bifIdx).conVesBdy;
        for vesIdx = 1:length(conVes)
            vesPnts = conVes(vesIdx).pnts;
            cPnt = mean(vesPnts);
            cAngle = conVes(vesIdx).avgAngle;
            dirVec = [sin(cAngle/180*pi),cos(cAngle/180*pi)];
            dPnt = cPnt + dirVec*40;
            drawArrow([cPnt(2),dPnt(2)],[cPnt(1),dPnt(1)]);
            %annotation('textarrow',[cPnt(2),dPnt(2)],[cPnt(1),dPnt(1)],'String','y = x ')
            pos = get(gca, 'Position');
            xArrow = [cPnt(2),dPnt(2)];
            yArrow = [cPnt(1),dPnt(1)];
            ta1 = annotation('textarrow',...
            (xArrow-min(xlim))/diff(xlim)*pos(3)+pos(1),(yArrow-min(ylim))/diff(ylim)*pos(4)+pos(2));
        end
    end
    hold off

    print([rootdir,'\RGBangles.png'],'-dpng','-r600');
end