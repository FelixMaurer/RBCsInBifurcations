%--------------------------------------------------------------------------
% Script Name : MergeNetWorkTrajectories
% Authors     : Felix Maurer
% Institution : Saarland University
% Email       : mail@felixmilanmaurer.com
% Date        : 2024
%
% Description :
%   This is a helper script to collect all trajectories from a series 
%   of videos recorded at the same position and merge them into one struct.
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
%% File Loop
cellTypes = {'Healthy_RBCs','Rigid_RBCs'};
clear trajfiles; k_trajfile = 1;
for IDXtype = 1:2
    cFolder = [rootDir '\' cellTypes{IDXtype}];
    filelist = dir(fullfile(cFolder, '\**\*ROI*traj.mat'));  %get list of files and folders in any subfolder
    filelist = filelist(~[filelist.isdir]);  %remove folders from list
    for fileIdx = 1:length(filelist)
        %% file location and identification declaration
        fileFolder = filelist(fileIdx).folder;
        fileName = filelist(fileIdx).name;
        filePath = [fileFolder '\' fileName];
        % extract ROI
        strIdx1 = strfind(fileName,'ROI_')+4;
        strIdx2 = strfind(fileName(strIdx1:end),'_');
        strIdx2 = strIdx2(1)+strIdx1-2;
        roiIdx = str2double(fileName(strIdx1:strIdx2));
        % extract cell type
        cellType = fileFolder(end-2:end);
        trajfiles(k_trajfile).folder = fileFolder;
        trajfiles(k_trajfile).name = fileName;
        trajfiles(k_trajfile).cellType = cellType;
        trajfiles(k_trajfile).roiIdx = roiIdx;
        k_trajfile = k_trajfile+1;
    end
end
%% go through folders and combine data
roiIdc = sort(unique([trajfiles.roiIdx]));
for roiIdx = roiIdc
    for cellTypeIdx = 1:2
        cfolder = [rootDir '\' cellTypes{cellTypeIdx}];
        clear clu;
        k_clu = 1;
        for fileIdx = 1:length(trajfiles)
            fileName = trajfiles(fileIdx).name;
            if strcmp(cfolder,trajfiles(fileIdx).folder) && trajfiles(fileIdx).roiIdx == roiIdx
                % -> they belong together
                filePath = [trajfiles(fileIdx).folder '\' trajfiles(fileIdx).name];
                load(filePath);
                for IDXtraj = 1:length(traj)
                    clu(k_clu).points = traj(IDXtraj).points;
                    k_clu = k_clu+1;
                end
            end
        end
        % save cluster
        try
            save([cfolder '\' 'ROI_' sprintf('%d',roiIdx) '_network_merge.mat'],'clu');
        end
    end
end
%% debugging
close all
if 0 
    figure
    hold on
    for cluIdx = 1:length(clu)
        pnts = clu(cluIdx).points;
        plot3(pnts(:,1),pnts(:,2),pnts(:,3));
        view(2)
    end
    hold off
end

