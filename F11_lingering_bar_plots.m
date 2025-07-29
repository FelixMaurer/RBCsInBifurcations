%--------------------------------------------------------------------------
% Script Name : F11_lingering_bar_plots.m
% Authors     : Felix Maurer, Yazdan Rashidi
% Institution : Saarland University
% Email       : mail@felixmilanmaurer.com
% Date        : 2024
%
% Description :
%   This script creates bar plots for comparison of lingering
%   quantification between datasets.
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
%% source
addpath('src');
%% load data and process
load("lingering_analysis_results.mat")
close all
data1 = (allRBCVelMoth-allRBCVelBif)./allRBCVelMoth;
data2 = (allWBCVelMoth-allWBCVelBif)./allWBCVelMoth;
cond2 = data2 > 0 & data1 > 0;
cond1 = cond2;

data1 = allRBCVelBif./allRBCVelMoth;
data2 = allWBCVelBif./allWBCVelMoth;

data(1).dat = data1(cond1);
data(2).dat = data2(cond2);
allData = [];
c = [];
for dataIdx = 1:length(data)
    allData = [allData,data(dataIdx).dat];
    c = [c,dataIdx*ones(1,length(data(dataIdx).dat))];
end
for dataIdx1 = 1:length(data)
    for dataIdx2 = 1:length(data)
        [h,p] = ttest(allData(c==dataIdx1),allData(c==dataIdx2))
    end
end

boxplot(allData,c,'Notch','on','OutlierSize',0.1,'Colors','k','BoxStyle','outline','Widths',.3)
hold off
xticklabels({'healthy', 'rigid'})
xlabel('population')
ylabel('$\frac{t_\mathrm{residence}^R}{t_\mathrm{residence}^H}-1$','interpreter','latex')
hold on
plot(1+(rand(1,sum(c==1))-0.5)*0.25 ...
    ,allData(c==1),'.','MarkerSize',15)
plot(2+(rand(1,sum(c==2))-0.5)*0.25 ...
    ,allData(c==2),'.','MarkerSize',15)
plot(3+(rand(1,sum(c==3))-0.5)*0.25 ...
    ,allData(c==3),'.','MarkerSize',15)
hold off
