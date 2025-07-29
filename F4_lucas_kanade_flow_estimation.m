%--------------------------------------------------------------------------
% Script Name : F4_lucas_kanade_flow_estimation.m
% Author      : Felix Maurer
% Institution : Saarland University
% Email       : mail@felixmilanmaurer.com
% Date        : 2024
%
% Description :
%   This script estimates the optical flow of a video using the Lucas
%   Kanade method.
% 
%   Lucas, B. D. & Kanade, T. (1981). 
%   An Iterative Image Registration Technique with an Application to Stereo Vision. 
%   In Proceedings of the 7th IJCAI (pp. 674–679).
%
%   The flow of intensities in the fluorescent plasma videos is used as an 
%   estimate for RBC velocities for a predictive search.
%
% Usage :
%   - the parent directory to all '_plasma.mj2' files that will be processed
%   should be provided in 'directory.txt'
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
%% find data directories
rootDir = char(readlines('directory.txt'));
filelist = dir(fullfile(rootDir, '**\*_plasma.mj2'));  % get list of files and folders in any subfolder
% go through files and analyse
for fileIdx = 1:length(filelist)
    fileName = filelist(fileIdx).name;
    fileFolder = filelist(fileIdx).folder;
    filePath = [fileFolder,'\',fileName];
    % display
    fprintf('working on -> %s\n',filePath);
    % VideoReader
    vrd = VideoReader(filePath);
    vrd.CurrentTime = 0;
    frame = vrd.readFrame;
    vrd.CurrentTime = 0;
    frameNum = vrd.NumFrames;
    frameNum = min([600,frameNum]); % only take 100 for testing.
    resizeFac = 0.3;
    frame = imresize(frame,resizeFac,'cubic');
    vidH = size(frame,1);
    vidW = size(frame,2);

    displayFlow = false;
    if displayFlow
        h = figure;
        movegui(h);
        hViewPanel = uipanel(h,'Position',[0 0 1 1],'Title','Plot of Optical Flow Vectors');
        hPlot = axes(hViewPanel);
    end
    % video loop
    flowX = zeros(vidH,vidW,'double');
    flowY = zeros(vidH,vidW,'double');
    timer = tic;
    t = zeros(1,frameNum);
    invertFrame = true; % true for bright plasma, dark cells
    fprintf('computing flow.\n');
    % make optical flow object
    opticFlow = opticalFlowLK;
    opticFlow.reset;
    opticFlow.NoiseThreshold = 0.001;
    % loop through all frames
    for frameIdx = 1:frameNum-2
        t(frameIdx) = toc;
        % read frame
        t_current = vrd.CurrentTime;
        if invertFrame
            frame = imcomplement(readFrame(vrd));
        else
            frame = readFrame(vrd);
        end
        % make grayscale
        if length(size(frame)) > 2
            frame = rgb2gray(frame);
        end
        % resize
        frame = imresize(frame,resizeFac,'cubic');

        flow = estimateFlow(opticFlow,double(frame));
        vx = flow.Vx;
        vy = flow.Vy;
        flowX = flowX + vx;
        flowY = flowY + vy;
        % status
        if ~mod(frameIdx-1,50)
            computeTime = toc(timer);
            fprintf('status: %.2f percent.\n',frameIdx/frameNum*100);
            fprintf('remaining: %.2f sec.\n', computeTime/frameIdx * (frameNum-frameIdx));
        end
    end
    flowX = flowX/frameNum;
    flowY = flowY/frameNum;
    % resize
    flowX = imresize(flowX,[vrd.height,vrd.width]);
    flowY = imresize(flowY,[vrd.height,vrd.width]);
    %%
    flowAbs = sqrt(flowX.^2+flowY.^2);
    % plot
    close all
    surf(flowAbs,'EdgeColor','none')
    colormap(jet)
    colorbar;
    clim([prctile(flowAbs(:),0.5),prctile(flowAbs(:),99.5)])
    view(2)
    %% save
    clear flow;
    flow.flowX = flowX;
    flow.flowY = flowY;
    fprintf('')
    save([filePath(1:end-4),'_flow_lk_scale.mat']);
end
