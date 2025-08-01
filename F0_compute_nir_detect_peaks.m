%--------------------------------------------------------------------------
% Script Name : F0_compute_nir_detect_peaks.m
% Author      : Felix Maurer
% Institution : Saarland University
% Email       : mail@felixmilanmaurer.com
% Date        : 2024
%
% Description :
%   This script computes the intensity outlier map, i.e., the Normalized
%   Intensity Ratio (NIR), for a sequence of image frames. The NIR at each
%   pixel (i,j) in frame k is defined as:
%
%       NIR(i,j,k) = |V(i,j,k) - mean_k(V)| / std_k(V)
%
%   where V(i,j,k) is the intensity at pixel (i,j) in frame k, and the
%   mean and standard deviation are computed over all pixels in frame k.
%
%   In the NIR image, peaks are detected using a function by Adi Natan,
%   https://stanford.edu/~natan/. The FastPeakFind function applys filters
%   to the frame, and finally a convolution using a filter matrix, which in
%   our case is a square of the approximate size of a cell in pixels. The
%   peak detection method then uses a local maxima detection in the
%   obtained frame D:
%
%       D(i,j) > D(i+a, j+b)   for all (a,b) ∈ {-1, 0, 1} \ (0,0)
%
% Usage :
%   - the parent directory to all .mj2 video files that will be processed
%   should be provided in 'directory.txt'
%   - the output contains all detected 3D points (x,y,k)
%
% Dependencies :
%   - FastPeakFind function
%   mathworks.com/matlabcentral/fileexchange/37388-fast-2d-peak-finder
%   - StatusBar function
%
% Reference :
%   This script is associated with the publication
%   Impact of Red Blood Cell Rigidity on in vivo Flow Dynamics and Lingering in Bifurcations
%   by Rashidi et al. 2025
% License :
%   MIT
%% source
addpath('src');
%% settings
% display peaks while processing
% useful for choosing parameters
flag_show_detection = true;
% output peak images for later visual feedback
% or documentation
flag_output_peaks_img = false;
% additionally save the workspace for later computations
flag_output_data = false;
%% file processing loop
clc; % clear command window
% parent directory with .mj2 files
rootDir = readlines('directory.txt');
filelist = dir(fullfile(rootDir, '**\*.mj2'));  %get list of files and folders in any subfolder
filelist = filelist(~[filelist.isdir]);  %remove folders from list
tic;
for fileIdx = 1:length(filelist)
    fileFolder = filelist(fileIdx).folder;
    fileName = filelist(fileIdx).name;
    filePath = [fileFolder '\' fileName];
    fprintf('working on -> %s.\n',filePath);
    if ~contains(fileName,'_plasma.mj2')
        fprintf('Gathering frames, this takes a while...\n');
        %% (0) initialize video file reader
        video_file = filePath;
        vrd = VideoReader(video_file);
        vrd.CurrentTime = 0;
        nFrames = vrd.NumFrames;
        %% (1) write video into one array
        STARTTIME = 0;
        TIMEINT = 3000;
        vrd.CurrentTime = STARTTIME;
        accframe = zeros(vrd.Height,vrd.Width,min(nFrames,TIMEINT*vrd.FrameRate),'uint16');
        for idxFrame = 1:nFrames
            frame = readFrame(vrd);
            gray_frame = frame(:,:,1);
            accframe(:,:,idxFrame) = gray_frame;
            if vrd.CurrentTime > TIMEINT+STARTTIME-1/vrd.FrameRate
                break;
            end
            %% show status
            StatusBar(idxFrame, nFrames, 'frame', 20);
        end
        accframe = accframe(:,:,1:idxFrame);
        fprintf('Done gathering frames.\n');
        fprintf('Analyzing intensity distribution, this may take a minute...\n');
        %% (1a) detect any larger change
        fprintf('Detect any larger intensity changes during the video...\n');
        % find larger changes in intensity, e.g. due to larger movements
        spatial_int = zeros(1,nFrames);
        for idxFrame = 1:nFrames
            spatial_int(idxFrame) = std(std(double(accframe(:,:,idxFrame))));
            %% show status
            StatusBar(idxFrame, nFrames, 'frame', 20);
        end
        start_frame = 1;
        end_frame = nFrames;
        filtsig = movmedian(spatial_int,50);
        absdiff = abs(diff(filtsig));
        height_thres = 20*mean(absdiff);
        changeIdc = find(absdiff>height_thres);
        if ~isempty(changeIdc)
            [peaks,pos] = findpeaks(absdiff,'MinPeakDistance',100);
            changePos = pos(peaks > height_thres);
            changePos = [1 changePos nFrames];
            % take the maximum interval length
            diffChangePos = diff(changePos);
            maxChangePosIdx = find(diffChangePos == max(diffChangePos),1,'first');
            maxChangePos = changePos(maxChangePosIdx);
            start_frame = maxChangePos+100;
            end_frame = changePos(maxChangePosIdx+1);
        end
        %% (2) Compute average and standard deviation for each pixel
        std_map = zeros(vrd.Height,vrd.Width,'double');
        avg_map = zeros(vrd.Height,vrd.Width,'double');
        for idx = 1:vrd.Height
            for idy = 1:vrd.Width
                int_values = double(accframe(idx,idy,start_frame:end_frame));
                int_values = int_values(int_values > 0);
                if isempty(int_values)
                    int_values = 127;
                end
                std_map(idx,idy) = std(int_values);
                avg_map(idx,idy) = mean(int_values);
            end
        end
        %% (3) propbability map
        STARTTIME = (start_frame-1)/vrd.FrameRate;
        ENDTIME = end_frame/vrd.FrameRate;
        vrd.CurrentTime = STARTTIME;
        clear allpoints; k = 1;
        for idxFrame = start_frame:end_frame
            frame = readFrame(vrd);
            gray_frame = frame(:,:,1);
            % compute the nir outlier probability
            prop = abs(double(gray_frame)-double(avg_map))./double(std_map);
            prop(isinf(prop)) = 0;
            prop(isnan(prop)) = 0;
            prop_map = prop;

            prop_map = 0.2*prop_map;
            prop_map(prop_map>1) = 1;
            prop_map(prop_map<0) = 0;
            %% gaussian filter
            % filter frame
            prop_map_img = uint8(255*prop_map);
            filt_frame = double(imgaussfilt(prop_map_img,2));
            %% find peaks
            d = filt_frame;
            thres = 56; % noise threshold, too small--> to much noise, too high--> less cells
            filt = ones(6,6); % similar to the cell size in pixel (diameter)
            edg = 10; % edge around image to ignore noise, value fo 10 is okay
            res = 1; % 1 is minimum, you can increase to test for enhancement
            peak_points = FastPeakFind(d, thres, filt ,edg, res);
            peak_x = peak_points(2:2:end);
            peak_y = peak_points(1:2:end);
            allpoints(idxFrame).peak = [peak_x, peak_y];
            %% display and output results
            if flag_show_detection
                imshow(frame);
                hold on
                plot(peak_y,peak_x,'.','Marker','*','Color',[1 0 0],'MarkerSize',10)
                hold off
                pause(0.01)
            end
            if flag_output_peaks_img
                outImage = insertShape(uint8(double(frame)/(2^16-1)*255),'circle',[peak_y,peak_x,8*ones(length(peak_y),1)],'LineWidth',1);
                imwrite(outImage,sprintf('%s\\frame%08d.bmp',outFOLDER,idxFrame));
            end

            if vrd.CurrentTime > TIMEINT+STARTTIME-1/vrd.FrameRate
                break;
            end
            %% show status
            StatusBar(idxFrame-start_frame+1, end_frame-start_frame, 'frame', 20);
        end
        %% save results
        save([filePath(1:end-4) '_peaks.mat'],'allpoints');
        %% save all settings
        if flag_output_data
            save([filePath(1:end-4) '_preak_detection.mat']);
        end
        %% show status
        StatusBar(fileIdx, length(filelist), 'file');
    end
end
