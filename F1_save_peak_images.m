%--------------------------------------------------------------------------
% Script Name : F1_save_peak_images.m
% Author      : Felix Maurer
% Institution : Saarland University
% Email       : mail@felixmilanmaurer.com
% Date        : 2024
%
% Description :
%   This script saves a projection of all peaks (x,y,k) along the k 
%   dimension after running the peak detection algorithm. 
%   Those images are used afterwards for manual alignment of the RBC
%   footage to the plasma footage.
%
% Usage :
%   - the parent directory to all '_peaks.mat' files that will be processed
%   should be provided in 'directory.txt'
%
% Dependencies :
%   - script StatusBar
% Reference :
%   This script is associated with the publication
%   Impact of Red Blood Cell Rigidity on in vivo Flow Dynamics and Lingering in Bifurcations
%   by Rashidi et al. 2025
% License :
%   MIT
%% source
addpath('src');
%% file processing
rootDir = char(readlines('directory.txt'));
filelist = dir(fullfile(rootDir, '**\*.mj2'));  
filelist = filelist(~[filelist.isdir]);
% load plasma image
plasma_filelist = dir(fullfile([rootDir '\**\*plasma*.tif']));
plasmaImg = imread([rootDir, '\', plasma_filelist(1).name]);
% cell type strings
cellTypes = {'Healthy_RBCs','Rigid_RBCs'};
% file processing loop
for fileIdx = 1:length(filelist)
    fileFolder = filelist(fileIdx).folder;
    fileName = filelist(fileIdx).name;
    filePath = [fileFolder '\' fileName];
    if contains(fileFolder,cellTypes)>0
        if contains(fileFolder,cellTypes{1})
            cellColor = [1 0.4 0];
            colorStr = 'red';
        else
            cellColor = [0 0.4 1];
            colorStr = 'blue';
        end
        %% load peaks
        load([filePath(1:end-4) '_peaks.mat']);
        pointImg = zeros(size(plasmaImg,1),size(plasmaImg,2),3,'uint8');
        for pointIdx = 1:length(allpoints)
            peaks = allpoints(pointIdx).peak;
            try
                pointImg = insertShape(pointImg,"circle",[peaks(:,2) peaks(:,1) 1*ones(size(peaks,1),1)],LineWidth=1,Color=colorStr);
            end
        end
        imwrite(pointImg,[filePath(1:end-4) '_peaks_img.png']);
    end
    %% show status
    StatusBar(fileIdx, length(filelist), 'file');
end
