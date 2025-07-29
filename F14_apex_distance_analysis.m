%--------------------------------------------------------------------------
% Script Name : F14_apex_distance_analysis
% Authors     : Felix Maurer, Yazdan Rashidi
% Institution : Saarland University
% Email       : mail@felixmilanmaurer.com
% Date        : 2024
%
% Description :
%   This script computes statistics on the distance of cell types to the
%   bifurcation apex.
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
%% file processing
parentDir = pwd;
allDirList = dir(fullfile([parentDir(1:end-5),'\**\bifTrajModes.mat']));
allDirList = {allDirList.folder};
%microScale = 0.65; % in case units are microns
microScale = 1; % units are pixels
HRBCApexDist = [];
RRBCApexDist = [];
HRBCApexDistErr = [];
RRBCApexDistErr = [];
for dirIdx = 1:length(allDirList)
    load([allDirList{dirIdx},'\bifLin.mat'],'bifurcations');
    for bifIdx = 1:length(bifurcations)
        try
            HRBCApex = bifurcations(bifIdx).ApexDistanceHealthy;
            RRBCApex = bifurcations(bifIdx).ApexDistanceRigid;

            HRBCApex = HRBCApex(HRBCApex>0)*microScale;
            RRBCApex = RRBCApex(RRBCApex>0)*microScale;

            edges = linspace(0,20,21);

            HRBCcnts = histcounts(HRBCApex,edges);
            RRBCcnts = histcounts(RRBCApex,edges);

            HRBCcnts = HRBCcnts/sum(HRBCcnts)/mean(diff(edges));
            RRBCcnts = RRBCcnts/sum(RRBCcnts)/mean(diff(edges));

            [pdfHRBC,~] = ksdensity(HRBCApex,linspace(0,20,201),'support','positive','boundarycorrection','reflection');
            [pdfRRBC,xPdf] = ksdensity(RRBCApex,linspace(0,20,201),'support','positive','boundarycorrection','reflection');


            HRBCcdf = cumsum(pdfHRBC);
            HRBCcdf = HRBCcdf/max(HRBCcdf);
            RRBCcdf = cumsum(pdfRRBC);
            RRBCcdf = RRBCcdf/max(RRBCcdf);

            maxDist = 4;
            xIdx = find(xPdf>maxDist,1,'first');
            HRBCVal = HRBCcdf(xIdx);
            RRBCVal = RRBCcdf(xIdx);

            HRBCVal =sum(HRBCApex<maxDist)/length(HRBCApex);
            RRBCVal =sum(RRBCApex<maxDist)/length(RRBCApex);
            if length(HRBCApex)>5 && length(RRBCApex)>5
                HRBCApexDist = [HRBCApexDist,HRBCVal];
                RRBCApexDist = [RRBCApexDist,RRBCVal];
                HRBCApexDistErr = [HRBCApexDistErr,HRBCVal/sqrt(length(HRBCApex))];
                RRBCApexDistErr = [RRBCApexDistErr,RRBCVal/sqrt(length(RRBCApex))];
            end
        
       end
    end
end

%%
% color definition
linewidth = 1;
lambda = linspace(0,1,100);
R = double(lambda > 0.5).*(lambda-0.5)/0.5;
G = 0*R;
B = flip(R);
smR = fit(linspace(-1,1,100)',R','smoothingspline','smoothingparam',1);
smG = fit(linspace(-1,1,100)',G','smoothingspline','smoothingparam',1);
smB = fit(linspace(-1,1,100)',B','smoothingspline','smoothingparam',1);

xData = RRBCApexDist; % daughter beginning is closer to bifurcation
yData = HRBCApexDist; % daughter end is further away from bifurcation

xErr = HRBCApexDistErr;
yErr = RRBCApexDistErr;

xErr(isnan(xErr)) = median(xErr,'omitnan');
yErr(isnan(yErr)) = median(yErr,'omitnan');

cond = ~isinf(xData)&~isnan(yData)&~isinf(xData)&~isnan(yData);

xData = xData(cond);
yData = yData(cond);
xErr = xErr(cond);
yErr = yErr(cond);

%axLims = [0 max([xData,yData])];
axLims = [0 1];

fig = figure();
hold on
x = axLims;


% plot linear model
mdl = fitlm(xData,yData,'intercept',true,'RobustOpts','on');
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
    scale = abs(xData-yData);dist = -(xVal-yVal)/sqrt(2)/0.1;
    dist(isnan(dist)) = 0; pointColor = [smR(dist),smG(dist),smB(dist)];
    pointColor(pointColor>1)=1;
    pointColor(pointColor<0)=0;
    errorbar(xData(pntIdx),yData(pntIdx),yErr(pntIdx),yErr(pntIdx),xErr(pntIdx),xErr(pntIdx),'s','color',pointColor,'MarkerSize',6,'MarkerFaceColor','none','MarkerEdgeColor','none','LineWidth',linewidth)
    s = scatter(xData(pntIdx),yData(pntIdx),20,'color',pointColor,'Marker','o','MarkerFaceColor',pointColor,'MarkerFaceAlpha',0.2,'MarkerEdgeColor',pointColor);
end
xlh = xlabel(['\fontname{cmss10}F_{\fontname{mwa_cmmi10}R}^{\fontname{cmss10}int}\fontname{cmss10}']);%,'Interpreter','tex')
ylh = ylabel(['\fontname{cmss10}F_{\fontname{mwa_cmmi10}H}^{\fontname{cmss10}int}\fontname{cmss10}']);%,'Interpreter','tex')
xline([2,5])
yline([2,5])
signRankTestLine;
exportgraphics(fig,'OUT-apexDistace.png','Resolution',600);
saveas(fig,'OUT-apexDistace.svg');