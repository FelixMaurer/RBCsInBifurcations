%--------------------------------------------------------------------------
% Script Name : F5_lucas_kanade_flow_evaluation.m
% Author      : Felix Maurer
% Institution : Saarland University
% Email       : mail@felixmilanmaurer.com
% Date        : 2024
%
% Description :
%   This script takes the estimated flow vectors and computes average flow
%   velocity vectors for each vessel.
% 
% Usage :
%   - the Lucas Kanade estimation is required as a previous step as well as
%   the geometry extraction.
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
%%
cellTypes = {'Healthy RBCs','Rigid RBCs'};
cellTypeIdx = 1; % 1: Healthy RBCs, 2:Rigid RBCs
maskFileName = 'Mask.png';
% find data directories
rootDir = char(readlines('directory.txt'));
rootfiles = dir(fullfile(rootDir, '**\',maskFileName));  % get list of files and folders in any subfolder
rootFolders = unique({rootfiles.folder});
timer = tic;
for rootIdx = 1:length(rootFolders)
    rootdir = rootFolders{rootIdx};
    % display
    fprintf('working on -> %s\n',rootdir);

    % extract ROI boundaries
    wallImg = imread([rootdir '\',maskFileName]);
    if length(size(wallImg)) > 2
        wallImg = wallImg(:,:,1);
    end
    wallBW = imbinarize(wallImg);   
    % search flow files
    flowfiles = dir(fullfile(rootdir, '**\*flow_lk_scale.mat'));  
    kFlow = 0;
    flowX = zeros(size(wallBW),'double');
    flowY = zeros(size(wallBW),'double');
    for fileIdx = 1:length(flowfiles)
        % file location and declaration
        fileFolder = flowfiles(fileIdx).folder;
        fileName = flowfiles(fileIdx).name;
        filePath = [fileFolder '\' fileName];
        % load flow
        load(filePath,'flow');
        % also flip
        flowY = flowY+flow.flowX;
        flowX = flowX+flow.flowY;

        kFlow = kFlow + 1;
    end
    % also flip
    flowY = flowY/kFlow;
    flowX = flowX/kFlow;
    %% analysis
    close all
    % velocity magnitude
    absVq = sqrt(flowX.^2+flowY.^2);
    absVq(isnan(absVq))=0;
    absVq = smoothdata(absVq,1,"gaussian",5);
    absVq = smoothdata(absVq,2,"gaussian",5);
    absVq = absVq.*double(wallBW);
    % angle
    flowX = smoothdata(flowX,1,"gaussian",5);
    flowX = smoothdata(flowX,2,"gaussian",5);
    flowY = smoothdata(flowY,1,"gaussian",5);
    flowY = smoothdata(flowY,2,"gaussian",5);
    angleq = -acos(flowY./sqrt(flowX.^2+flowY.^2)).*(1-2*(-flowX<0))/pi*180;
    angleq = angleq.*double(wallBW);
    warning ('on','all');
    [Xq,Yq] = meshgrid(1:size(wallBW,1),1:size(wallBW,2));
    % save everything
    velocity(1).X = Xq;
    velocity(1).Y = Yq;
    velocity(1).V = absVq;
    velocity(1).angle = angleq;

    %% combine ROI
    mask = ~wallBW;
    Vs = absVq; % LK factor
    angles = angleq;
    angles(mask==1) = nan;
    Vs(mask==1) = nan;
    
    suffix = 'LK_scale';
    % save data
    save([rootdir(1:end-length('converted_expo')),'velocity_magnitude_' suffix '.mat'],'Vs');
    save([rootdir(1:end-length('converted_expo')),'velocity_direction_' suffix '.mat'],'angles');

    %% plotting
    close all;
    frameRate = 170; % plasma framerate
    dist = 1; % must be 1 in all cases
    deltaT = 1/frameRate;
    microScale = 0.65; % micron/px: 20x is 0.65 µm/px, and 50x is 0.26 µm/px
    velScale = microScale * dist/deltaT;
    Vs = Vs * velScale;
    Velocity = Vs;
    % velocity asbsolute
    fig = figure;
    surf(Vs,'EdgeColor','none')
    colormap(jet);
    view(2);
    cbar = colorbar;
    cbar.Label.String = 'speed in [microns/s]';
    cbar.Label.Rotation = 90; % Vertical text
    cbar.Label.Position(1) = cbar.Label.Position(1) - 3; % Adjust horizontal offset if needed
    clim([prctile(Velocity(:),5),prctile(Velocity(:),95)]);
    set(cbar,'units','centimeters','position',[10*size(wallBW,2)/size(wallBW,1)+0.5 .5 .5 8])
    clim([prctile(Vs(:),3),prctile(Vs(:),97)])
    fig.Units = 'centimeters';
    fig.Position = [2 5 10*size(wallBW,2)/size(wallBW,1)+2 10]; % fixed position
    ax = gca;
    ax.Units = 'centimeters';
    ax.Visible = 'off';
    ax.Position(1:2) = 0; % flush
    ax.Position(4) = fig.Position(4); % full coverage in height
    ax.Position(3) = ax.Position(4) * size(wallBW,2)/size(wallBW,1);

    % saving
    print([rootdir,'\','velocity_magnitude_traj' suffix '.png'],'-dpng','-r600');

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
    clim([-180,180])
    cbar.Label.String = 'angle to x-axis in [°]';
    cbar.Label.Rotation = 90; % Vertical text
    cbar.Label.Position(1) = cbar.Label.Position(1) - 4; % Adjust horizontal offset if needed
    fig.Units = 'centimeters';
    fig.Position = [25 5 10*size(wallBW,2)/size(wallBW,1)+2 10]; % fixed position
    ax = gca;
    ax.Units = 'centimeters';
    ax.Visible = 'off';
    ax.Position(1:2) = 0; % flush
    ax.Position(4) = fig.Position(4); % full coverage in height
    ax.Position(3) = ax.Position(4) * size(wallBW,2)/size(wallBW,1);

    % saving
    print([rootdir,'\','velocity_direction_' suffix '.png'],'-dpng','-r600');

    % status
    time = toc(timer);
    fprintf('status: %.2f percent.\n',rootIdx/length(rootFolders)*100);
    fprintf('remaining time: %.2f sec.\n',time/rootIdx * (1-rootIdx/length(rootFolders)));
end