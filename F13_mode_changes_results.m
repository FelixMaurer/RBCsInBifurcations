%--------------------------------------------------------------------------
% Script Name : F13_mode_changes_results.m
% Authors     : Felix Maurer, Yazdan Rashidi
% Institution : Saarland University
% Email       : mail@felixmilanmaurer.com
% Date        : 2024
%
% Description :
%   This script processes the mode change analysis results and computes
%   statistics.
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
%% file declarations
rootDir = char(readlines("directory.txt"));
filelist = dir(fullfile([rootDir,'\**\bifTrajModes.mat']));
folders = {filelist.folder};
folders = unique(folders);
%% source
addpath('src');
%% processing
HRBCthres = 2;
RRBCthres = 2;
zPositions = 0;
HRBCModePerc = [];
RRBCModePerc = [];
HRBCModePercErr = [];
RRBCModePercErr = [];
HRBCthres = [];
RRBCthres = [];
HRBCModesVes = [];
RRBCModesVes = [];
HRBCFloDaugh = [];
RRBCFloDaugh = [];
HRBCFloMoth = [];
RRBCFloMoth = [];
for rootIdx = 1:length(folders)
    cfolder = folders{rootIdx};
    load([cfolder,'\bifTrajModes.mat'],'bifurcations');
    cBif = bifurcations;
    for bifIdx = 1:length(cBif)
        HRBCModes = bifurcations(bifIdx).HRBCMode;
        RRBCModes = bifurcations(bifIdx).RRBCMode;
        HRBCthres = [HRBCthres,bifurcations(bifIdx).HRBCModeThres];
        RRBCthres = [RRBCthres,bifurcations(bifIdx).RRBCModeThres];
        HRBCModePerc = [HRBCModePerc,sum(HRBCModes)/length(HRBCModes)*100];
        RRBCModePerc = [RRBCModePerc,sum(RRBCModes)/length(RRBCModes)*100];

        HRBCModePercErr = [HRBCModePercErr,1/sqrt(length(HRBCModes))];
        RRBCModePercErr = [RRBCModePercErr,1/sqrt(length(RRBCModes))];

        conVes = cBif(bifIdx).conVesBdy;
        orVes = cBif(bifIdx).orVes;
        for vesIdx = 1:length(conVes)
            HRBCVel = conVes(vesIdx).HRBCVel;
            HRBCMode = logical(conVes(vesIdx).HRBCMode);

            RRBCVel = conVes(vesIdx).RRBCVel;
            RRBCMode = logical(conVes(vesIdx).RRBCMode);

            HRBCVel = HRBCVel.*cBif(bifIdx).velScaleHRBC;
            RRBCVel = RRBCVel.*cBif(bifIdx).velScaleRRBC;
            if vesIdx==orVes

                HRBCFloMoth = [HRBCFloMoth,mean(HRBCVel(~HRBCMode))];
                RRBCFloMoth = [RRBCFloMoth,mean(RRBCVel(~RRBCMode))];
            else

                HRBCFloDaugh = [HRBCFloDaugh,mean(HRBCVel(~HRBCMode))];
                RRBCFloDaugh = [RRBCFloDaugh,mean(RRBCVel(~RRBCMode))];
            end
        end
    end
end
HRBCModesVes(isnan(HRBCModesVes)) = 0;
RRBCModesVes(isnan(RRBCModesVes)) = 0;
HRBCModesVes = HRBCModesVes * 100;
RRBCModesVes = RRBCModesVes * 100;
%% color
% color definition
linewidth = 1;
lambda = linspace(0,1,100);
R = double(lambda > 0.5).*(lambda-0.5)/0.5;
G = 0*R;
B = flip(R);
smR = fit(linspace(-1,1,100)',R','smoothingspline','smoothingparam',1);
smG = fit(linspace(-1,1,100)',G','smoothingspline','smoothingparam',1);
smB = fit(linspace(-1,1,100)',B','smoothingspline','smoothingparam',1);


%% healthy stds

xData = RRBCModePerc;%/mean(HRBCVels(:,1),'omitnan');
yData = HRBCModePerc;%/mean(HRBCVels(:,2),'omitnan');

xErr = RRBCModePercErr*100;
yErr = HRBCModePercErr*100;


xErr(isnan(xErr)) = median(xErr,'omitnan');
yErr(isnan(yErr)) = median(yErr,'omitnan');


cond = ~isinf(xData)&~isnan(yData)&~isinf(xData)&~isnan(yData)&xData>0&yData>0;

xData = xData(cond);
yData = yData(cond);
xErr = xErr(cond);
yErr = yErr(cond);


axLims = [min([xData,yData]) max([xData,yData])];
axLims = [0 100];

fig = figure();
hold on
x = axLims;


% plot linear model
weights = 1./(xErr+yErr); weights = weights/mean(weights);  
mdl = fitlm(xData,yData,'intercept',true,'RobustOpts','off','weights',weights);
lm_slope = mdl.Coefficients.Estimate(2);
lm_intcp = mdl.Coefficients.Estimate(1);
lm_var = mdl.Coefficients.SE(2);
lm_intvar = mdl.Coefficients.SE(1);
xLine = linspace(-1000,1000,2);
lineMid = (lm_slope)*xLine + lm_intcp;
lineTop = (lm_slope+lm_var)*xLine + lm_intcp+lm_intvar;
lineBot = (lm_slope-lm_var)*xLine + lm_intcp-lm_intvar;
fill([xLine,flip(xLine)],[lineTop, fliplr(lineBot)],[0.85 0.85 0.85],'EdgeColor','none')
lineTop = (lm_slope+lm_var)*xLine + lm_intcp-lm_intvar;
lineBot = (lm_slope-lm_var)*xLine + lm_intcp+lm_intvar;
fill([xLine,flip(xLine)],[lineTop, fliplr(lineBot)],[0.85 0.85 0.85],'EdgeColor','none')


plot(x,x,'--','LineWidth',1.4,'color',[0.2 0.2 0.2])

linewidth = 1;
for pntIdx = 1:length(xData)
    xVal = xData(pntIdx);
    yVal = yData(pntIdx);
    scale = abs(xData-yData);dist = -(xVal-yVal)/sqrt(2)/mean(abs(scale(~isnan(scale))));
    scale = abs(xData-yData);dist = -(xVal-yVal)/sqrt(2)/mean(abs(scale(~isnan(scale))));
    dist
    dist(isnan(dist)) = 0; pointColor = [smR(dist),smG(dist),smB(dist)];
    pointColor(pointColor>1)=1;
    pointColor(pointColor<0)=0;
    errorbar(xData(pntIdx),yData(pntIdx),yErr(pntIdx),yErr(pntIdx),xErr(pntIdx),xErr(pntIdx),'s','color',pointColor,'MarkerSize',6,'MarkerFaceColor','none','MarkerEdgeColor','none','LineWidth',linewidth)
    s = scatter(xData(pntIdx),yData(pntIdx),20,'color',pointColor,'Marker','o','MarkerFaceColor',pointColor,'MarkerFaceAlpha',0.2,'MarkerEdgeColor',pointColor);
end

xlh = xlabel(['\fontname{cmss10}F_{\fontname{cmss10}int}^{\fontname{mwa_cmmi10}R}\fontname{cmss10} [%]']);%,'Interpreter','tex')
ylh = ylabel(['\fontname{cmss10}F_{\fontname{cmss10}int}^{\fontname{mwa_cmmi10}H}\fontname{cmss10} [%]']);%,'Interpreter','tex')

xticks([0 50 100])
yticks([0 50 100])
xticklabels({0,'',100})
yticklabels({0,'',100})



signRankTestLine;

ylh.Position(1) = ylh.Position(1) + 12;
xlh.Position(2) = xlh.Position(2) + 5;


ax.TickLength = [0.03 0.03];

title('migration','interpreter','latex')

exportgraphics(fig,'OUT-mode-changes.png','Resolution',600);
saveas(fig,'OUT-mode-changes.svg');
save('modeChangeResults.mat')
