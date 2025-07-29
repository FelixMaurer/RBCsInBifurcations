%--------------------------------------------------------------------------
% Script Name : F2_extract_and_check_alignment.m
% Author      : Felix Maurer
% Institution : Saarland University
% Email       : mail@felixmilanmaurer.com
% Date        : 2024
%
% This script visualizes the manual alignment of footage from different
% color channels. Usually the RBC fluorescence to the plasma fluorescence
% frame. The manual alignment has to be done previously. 
% The alignment vector for each file is saved for further analysis.
%
% Usage :
%   - after performing manual alignment each video file has the suffix
%   '_offset.png'. The offset is contained in the position of a rectangle 
%   in the frame compared to a rectangle in a reference image
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
%% find files
% parent dir
rootDir = char(readlines('directory.txt'));
% find plasma images
plasmaImgFiles = dir(fullfile([rootDir '\**\*plasma*.tif']));
%% file processing loop
for plasmaFileIdx = 1:length(plasmaImgFiles)
    plasmaFileFolder = plasmaImgFiles(plasmaFileIdx).folder;
    plasmaFileName = plasmaImgFiles(plasmaFileIdx).name;
    plasmaFilePath = [plasmaFileFolder,'\',plasmaFileName];
    % find align images
    clear alignFiles;
    alignFiles(1).files = dir(fullfile([plasmaFileFolder '\Healthy RBCs\*_offset.png']));
    alignFiles(2).files = dir(fullfile([plasmaFileFolder '\Rigid RBCs\*_offset.png']));
    alignFiles(1).peakFiles = dir(fullfile([plasmaFileFolder '\Healthy RBCs\*_peaks_img.png']));
    alignFiles(2).peakFiles = dir(fullfile([plasmaFileFolder '\Rigid RBCs\*_peaks_img.png']));
    for typeIdx = 1:2
        % read plasma image
        plasmaImg = im2gray(imread(plasmaFilePath));
        % find offset align files
        files = alignFiles(typeIdx).files;
        peakFiles = alignFiles(typeIdx).peakFiles;
        for fileIdx = 1:length(files)
            fileFolder = files(fileIdx).folder;
            fileName = files(fileIdx).name;
            filePath = [fileFolder,'\',fileName];
            % read reference image
            try
                refImg = im2gray(imread([fileFolder,'\refImg.png']));
            catch
                fprintf('reference refImg.png not found.\n');
                return
            end
            % read align offset image
            offImg = im2gray(imread(filePath));
            % binarize
            bwRefImg = imbinarize(refImg);
            bwOffImg = imbinarize(offImg);
            % find center of mass
            stats = regionprops(bwRefImg,'centroid');
            refPos = stats.Centroid(1:2);
            stats = regionprops(bwOffImg,'centroid');
            offPos = stats.Centroid(1:2);
            if flag_plot_positions
                plot(refPos(1),refPos(2),'.')
                hold on
                plot(offPos(1),offPos(2),'.')
                hold off
            end
            % load peak images
            peakImg = imread([fileFolder,'\',peakFiles(fileIdx).name]);
            shiftImg = zeros(size(peakImg));
            offset = (offPos-refPos);
            offset = flip(offset);
            % save alignment vector
            alignVec = offset;
            save(sprintf('%s_alignVec.mat',filePath(1:end-4)),'alignVec');
            if offset(1) > 0
                shiftImg(offset(1):size(shiftImg,1),:,:) = peakImg(1:size(shiftImg,1)-offset(1)+1,:,:);
            end
            if offset(1) < 0
                shiftImg(1:size(shiftImg,1)+offset(1)+1,:,:) = peakImg(-offset(1):size(shiftImg,1),:,:);
            end
            if offset(2) > 0
                shiftImg(:,offset(2):size(shiftImg,2),:) = shiftImg(:,1:size(shiftImg,2)-offset(2)+1,:);
            end
            if offset(2) < 0
                shiftImg(:,1:size(shiftImg,2)+offset(2)+1,:) = shiftImg(:,-offset(2):size(shiftImg,2),:);
            end
            plasmaImg = plasmaImg + uint16((2^16-1)*imbinarize(shiftImg));
            %% show status
            StatusBar(fileIdx, length(files), 'file');  
        end
        imwrite(plasmaImg,sprintf('%s_type_%d_alignment.png',plasmaFilePath(1:end-4),typeIdx));
    end
end