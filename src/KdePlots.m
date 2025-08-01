%--------------------------------------------------------------------------
% Script Name : KdePlots.m
% Authors     : Felix Maurer
% Institution : Saarland University
% Email       : mail@felixmilanmaurer.com
% Date        : 2024
%
% Description :
%   This is a tool for making 2D KDEs
%
% Usage :
%   - The script requires previous full network tracking and Lucan Kanade
%   flow estimation
%
% Dependencies :
%
% License :
%   MIT
%% processing
timeScale = 1000/50;

[yHRBCedgeIdc,xHRBCedgeIdc] = find(HRBCedge);
[yRRBCedgeIdc,xRRBCedgeIdc] = find(RRBCedge);

% HRBCs
image = zeros(size(HRBCkde),'logical');
for edgeIdx = 1:length(xHRBCedgeIdc)
    image(round(yHRBCedgeIdc(edgeIdx)),round(xHRBCedgeIdc(edgeIdx))) = 1;
end
bwDist = image;
%bwDist(1,:) = 1;bwDist(end,:) = 1;bwDist(:,1) = 1;bwDist(:,end) = 1;
bwDist = imfill(bwDist,'holes');
[yAbsMax,xAbsMax] = find(HRBCkde == max(HRBCkde(:)),1,'first');
bwDist = imfill(imcomplement(bwDist),[yAbsMax,xAbsMax])-imcomplement(bwDist);
bwDistT = 0*bwDist;
for translateIdx = 0:50
    bwDistT = bwDistT|imtranslate(flip(bwDist,2),[round(size(image,2)/2)-xAbsMax+translateIdx,0]);
end
bwDistT = bwDistT&bwDist;
bwDistT = imclose(bwDistT,strel('disk',40,6));
[xBW,yBW] = find(~bwDist);
HRBCkdeT = HRBCkde;
for pntIdx = 1:length(xBW)
    HRBCkdeT(xBW(pntIdx),yBW(pntIdx)) = 0;
end

% RRBCs
image = zeros(size(RRBCkde),'logical');
for edgeIdx = 1:length(xRRBCedgeIdc)
    image(round(yRRBCedgeIdc(edgeIdx)),round(xRRBCedgeIdc(edgeIdx))) = 1;
end
bwDist = image;
%bwDist(1,:) = 1;bwDist(end,:) = 1;bwDist(:,1) = 1;bwDist(:,end) = 1;
bwDist = imfill(bwDist,'holes');
[yAbsMax,xAbsMax] = find(RRBCkde == max(RRBCkde(:)),1,'first');
bwDist = imfill(imcomplement(bwDist),[yAbsMax,xAbsMax])-imcomplement(bwDist);
bwDistT = 0*bwDist;
for translateIdx = 0:50
    bwDistT = bwDistT|imtranslate(flip(bwDist,2),[xAbsMax-round(size(image,2)/2)+translateIdx,0]);
end
bwDistT = bwDistT&bwDist;
bwDistT = imclose(bwDistT,strel('disk',40,6));

[xBW,yBW] = find(~bwDist);
RRBCkdeT = RRBCkde;
for pntIdx = 1:length(xBW)
    RRBCkdeT(xBW(pntIdx),yBW(pntIdx)) = 0;
end
% 
% surf(RRBCkdeT,'EdgeColor','none');
% return
%% define edge again
edgeThres = mean(0.5*mean(HRBCkde(HRBCkde>0)));
HRBCedge = edge(HRBCkdeT>edgeThres);
[yHRBCedge,xHRBCedge] = find(HRBCedge);
xHRBCedge = x(xHRBCedge);
yHRBCedge = y(yHRBCedge);

edgeThres = mean(0.5*mean(RRBCkde(RRBCkde>0)));
RRBCedge = edge(RRBCkdeT>edgeThres);
[yRRBCedge,xRRBCedge] = find(RRBCedge);
xRRBCedge = x(xRRBCedge);
yRRBCedge = y(yRRBCedge);
%%
minTime = 0;
minDist = 0;

maxTime = max(x(:))*timeScale;
maxDist = max(y(:))*microScale;
%maxDist = 200;
        plot_kde_analysis = true;
        if plot_kde_analysis
            %close all;
            kdefig = figure;
            tiledlayout(3,1);
            nexttile;
            hold on;
            surf(x*timeScale,y*microScale,HRBCkdeT,'EdgeColor','none')
            title(sprintf('healthy N = %d',length(HRBCDist)))
            plot3(HRBCTime*timeScale,HRBCDist*microScale,10000*ones(1,length(HRBCTime)),'.r','MarkerSize',10)
            plot3(xHRBCedge*timeScale,yHRBCedge*microScale,10000*ones(1,length(xHRBCedge)),'.g','MarkerSize',7)
            view(2)
            xlim([minTime,maxTime]);
            xlabel('time [ms]');
            ylim([minDist,maxDist]);
            ylabel(['maximum distance [' char(181) 'm]']);
            %ylim([100 160]);
            nexttile;
            hold on;
            surf(x*timeScale,y*microScale,RRBCkdeT,'EdgeColor','none')
            title(sprintf('rigid N = %d',length(RRBCDist)))
            plot3(RRBCTime*timeScale,RRBCDist*microScale,10000*ones(1,length(RRBCTime)),'.r','MarkerSize',10)
            plot3(xRRBCedge*timeScale,yRRBCedge*microScale,10000*ones(1,length(xRRBCedge)),'.g','MarkerSize',7)
            view(2)
            %xlim([minTime,500]);
            xlabel('time [ms]');
            %ylim([minDist,100]);
            ylabel(['maximum distance [' char(181) 'm]']);
            %ylim([100 160]);
            % figure;hold on;
            % plot(xHRBCedge,yHRBCedge,'.')
            % plot(x(HRBCKdeMaxPos),y(1:length(HRBCKdeMaxPos)),'o-b')
            % %plot(HRBCmaxFit(y(1:length(HRBCKdeMaxPos))),y(1:length(HRBCKdeMaxPos)),'-k')
            % %plot(HRBCMax(1),HRBCMax(2),'xg','LineWidth',3,'MarkerSize',16)
            % 
            % plot(xRRBCedge,yRRBCedge,'.')
            % plot(x(RRBCKdeMaxPos),y(1:length(RRBCKdeMaxPos)),'o-r')
            % %plot(RRBCmaxFit(y(1:length(RRBCKdeMaxPos))),y(1:length(RRBCKdeMaxPos)),'-k')
            % %plot(RRBCMax(1),RRBCMax(2),'xg','LineWidth',3,'MarkerSize',16)
            % 
            xlim([minTime,maxTime]);
            ylim([minDist,maxDist]);
            figkde = gcf;
            

            yline(vesLength*microScale);

            %% some data
            % compute number center of mass (expectation value)
            [X,Y] = meshgrid(x,y);
            HRBCxCtr = (round(sum(sum(X.*HRBCkde))/sum(HRBCkde(:))))*timeScale;
            HRBCyCtr = (round(sum(sum(Y.*HRBCkde))/sum(HRBCkde(:))))*microScale;
            RRBCxCtr = (round(sum(sum(X.*RRBCkde))/sum(RRBCkde(:))))*timeScale;
            RRBCyCtr = (round(sum(sum(Y.*RRBCkde))/sum(RRBCkde(:))))*microScale;

            % compute velocity center of mass (expectation value)
            [X,Y] = meshgrid(x,y);
            Vel = (Y*microScale)./(X*timeScale);
            Vel(isnan(Vel)|isinf(Vel)) = 0;
            HRBCkdeTZeros = HRBCkdeT;
            HRBCkdeTZeros(isnan(HRBCkdeTZeros)) = 0;
            RRBCkdeTZeros = RRBCkdeT;
            RRBCkdeTZeros(isnan(RRBCkdeTZeros)) = 0;
            HRBCVelCtr = ((sum(sum(Vel.*HRBCkdeTZeros))/sum(HRBCkdeTZeros(:))))*1000; % um/s
            RRBCVelCtr = ((sum(sum(Vel.*RRBCkdeTZeros))/sum(RRBCkdeTZeros(:))))*1000; % um/s

            %HRBCVelCtr = HRBCyCtr/HRBCxCtr*1000; % um/s
            %RRBCVelCtr = RRBCyCtr/RRBCxCtr*1000; % um/s

            
            nexttile;
            hold on
            plot3(xHRBCedge*timeScale,yHRBCedge*microScale,9000*ones(1,length(xHRBCedge)),'.r','MarkerSize',7)
            plot3(xRRBCedge*timeScale,yRRBCedge*microScale,10000*ones(1,length(xRRBCedge)),'.b','MarkerSize',7)
            plot3(RRBCxCtr,RRBCyCtr,11000*ones(1,length(HRBCxCtr)),'x','color',[0.1 0.2 0.8],'MarkerSize',12)
            plot3(HRBCxCtr,HRBCyCtr,12000*ones(1,length(RRBCxCtr)),'x','color',[0.8 0.2 0.1],'MarkerSize',12)
            
            view(2)
            hold off
            xlim([minTime,maxTime]);
            xlabel('time [ms]');
            ylim([minDist,maxDist]);
            ylabel(['maximum distance [' char(181) 'm]']);

            title(sprintf('H: v = %.2f, R: v = %.2f',HRBCVelCtr,RRBCVelCtr));
            figkde.Position = [560 428.0000*0.2 253 420.0000*1.7];

            
            print(sprintf('%s\\bif%dves%d-kdes.png',rootdir,bifIdx,vesIdx),'-dpng','-r600');
            close(figkde);

            avoidVariable = 'fig';
            save(sprintf('%s\\bif%dves%d.mat',rootdir,bifIdx,vesIdx),'-regexp', ['^(?!', avoidVariable,'$).']);
        end