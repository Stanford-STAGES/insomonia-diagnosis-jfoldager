%% Initials
startup
% This script houses the sub-categorization of insomnia subjects
% using various clustering approaches.
% Input: genetic data, psg variables
% Output: multiple plots from kmeans, gmms, and pcha
%% Load data
% % % % % % % % % % % % % % % % % % % % % % % % % % % Without mental health
load V1newdata.mat
riskP1      = sum(X.rs6589988,2); 
XsP1        = X;
YsP1        = table();
YsP1.ARI    = Y.ARI;
YsP1.AHI    = Y.AHI;
YsP1.AW25   = Y.AW25;
YsP1.AWN125 = Y.AWN125;
YsP1.PLMI   = Y.PLMI;
YsP1.SE     = Y.SE;
YsP1.SOL    = Y.SOL;
YsP1.TST    = Y.TST; 
CovP1       = [Y.AGE,Y.BMI,Y.SEX]; 
% % % % % % % % % % % % % % % % % % % % % % % % % % % % With mental health
load V2newdata.mat
XsP2        = X;
riskP2      = sum(X.rs6589988,2);
YsP2        = table();
YsP2.ARI    = Y.ARI;
YsP2.AHI    = Y.AHI;
YsP2.AW25   = Y.AW25;
YsP2.AWN125 = Y.AWN125;
YsP2.PLMI   = Y.PLMI;
YsP2.SE     = Y.SE;
YsP2.SOL    = Y.SOL;
YsP2.TST    = Y.TST; 
INS         = Y.INS;
CovP2       = [Y.AGE,Y.BMI,Y.SEX]; 
RISKSP2     = X;
%% PCA
[coeffYP1,scoreYP1,latentYP1]   = pca(zscore(YsP1{:,:}));
[coeffXP1,scoreXP1,latentXP1]   = pca(zscore(XsP1{:,:}));
[coeffYP2,scoreYP2,latentYP2]   = pca(zscore(YsP2{:,:}));
[coeffXP2,scoreXP2,latentXP2]   = pca(zscore(XsP2{:,:}));
Tf1                             = sum(isoutlier(scoreXP1(:,1:2)),2) > 0;
Tf2                             = sum(isoutlier(scoreXP2(:,1:2)),2) > 0;
XsP1(Tf1,:)                     = [];
YsP1(Tf1,:)                     = [];
XsP2(Tf2,:)                     = [];
YsP2(Tf2,:)                     = [];
INS(Tf2,:)                      = [];
RISKSP2(Tf2,:)                  = [];
CovP1(Tf1,:)                    = [];
CovP2(Tf2,:)                    = [];
YzP1                            = zscore(YsP1{:,:});
YzP2                            = zscore(YsP2{:,:});
% For AA (v1)
YsP1arc                         = YzP1;
riskP1arc                       = riskP1;
YsP1arc(sum(isoutlier(YzP1,'mean'),2) > 0,:)   = [];
riskP1arc(sum(isoutlier(YzP1,'mean'),2) > 0,:) = [];
% For AA (v2)
YsP2arc     = YzP2;
riskP2arc   = riskP2;
YsP2arc(sum(isoutlier(YzP2,'mean'),2) > 0,:)    = [];
riskP2arc(sum(isoutlier(YzP2,'mean'),2) > 0,:)  = [];
INS(sum(isoutlier(YzP2,'mean'),2) > 0)          = [];
RISKSP2(sum(isoutlier(YzP2,'mean'),2) > 0,:)    = [];
%% Kmeans
clusterData = YsP2arc(INS == 0,:);
M           = size(clusterData,2);
N           = size(clusterData,1);
R           = 20;
K           = 20;
optC        = nan(R,1);
vals        = nan(R,K);
% Estimate number of clusters
for r = 1:R
    [Train, Test]   = crossvalind('HoldOut', size(clusterData,1), 0.5);
    eva             = evalclusters(clusterData(Train,:),'gmdistribution','DaviesBouldin','KList',1:K);
    optC(r)         = eva.OptimalK;
    vals(r,:)       = eva.CriterionValues;
end
% plot 1
figure, plot(vals')
box off, grid on
xlabel('Number of clusters')
ylabel('Davies Bouldin index')
% plot 2
figure, errorbar(nanmean(vals),nanstd(vals))
box off, grid on
xlabel('Number of clusters')
ylabel('Davies Bouldin index')

K = 7;
% Final clustering:
[labels,centroids] = kmeans(clusterData,K);
figure,  bar(centroids), xlabel('Cluster'), ylabel('Coordinate'), 
grid on, box off
legend(YsP1.Properties.VariableNames,'Location','northeastoutside');
figure, histogram(labels), xlabel('Cluster'), ylabel('Count')
grid on, box off
% GMM
[coeffYP2,scoreYP2,latentYP2]   = pca(clusterData);
GMModel = fitgmdist(scoreYP2(:,1:8),5,'CovarianceType','diagonal','Replicates',100);

figure,  bar(GMModel.mu), xlabel('Cluster'), ylabel('Coordinate'), 
grid on, box off
legend(YsP1.Properties.VariableNames,'Location','northeastoutside');
figure, histogram(labels), xlabel('Cluster'), ylabel('Count')
grid on, box off

%% GMM
R           = 100;
K           = 20;
trainNLog   = nan(R,K);
testNLog    = nan(R,K);
gmmData     = YsP2arc(INS == 0,:);
optC        = nan(R,1);
vals        = nan(R,K);
for r = 1:R
    [Train, Test] = crossvalind('HoldOut', size(gmmData,1), 0.25);
    eva = evalclusters(gmmData(Train,:),'gmdistribution','DaviesBouldin','KList',1:K);
    optC(r)         = eva.OptimalK;
    vals(r,:)       = eva.CriterionValues;
end
figure, plot(vals')
box off, grid on
xlabel('Number of clusters')
ylabel('Davies Bouldin index')

figure, errorbar(nanmean(vals),nanstd(vals))
box off, grid on
xlabel('Number of clusters')
ylabel('Davies Bouldin index')
%% PCHA
aaData      = YsP2arc;
M           = size(aaData,2);
N           = size(aaData,1);
R           = 10;
K           = 10;
sigmaKmeans = nan(R,K);
SSE         = nan(R,K);
RSS         = nan(R,K);
varexpl     = nan(R,K);
for r = 1:R
    for k = 1:K
        [~,~,~,SSE(r,k),varexpl(r,k)]=PCHA(aaData',k);
        RSS(r,k)            = sqrt(SSE(r,k));
        sigmaKmeans(r,k)    = SSE(r,k)/N;
    end
end

aaData          = YsP1arc;
[XC,S,~,~,~]    = PCHA(aaData',7);
AAplot(XC,S,3,XsP1.rs6589988((~sum(isoutlier(YzP1,'mean'),2) > 0))...
    ,{YsP1.Properties.VariableNames(1:end),{'GS','Low TST, SE','High ARI, AHI', ...
    'Long arousals','LS, low SE, TST','Long SOL','High PLMI'}})
title('Dosage of rs6589988')

YsP1.AA = nan(size(YsP1,1),size(S,1));
YsP1.AA(~sum(isoutlier(YzP1,'mean'),2) > 0,:) = S';
writetable(YsP1,'AAdata.csv')

aaData      = YsP2arc(INS == 1,:);
[XC,S,C,~,varexpl] = PCHA(aaData',7);
AAplot(XC,S,7,RISKSP2.rs6589988(INS == 1)...
    ,{YsP2.Properties.VariableNames(1:end),{'GS','Light sleeper','Long arousals', ...
    'Short sleeper','High ARI, AHI','Light sleep, PLMI ','Low TST, SE. High SOL'}})
title('Dosage of rs6589988')

YsP2.AA = nan(size(YsP2,1),size(S,1));
YsP2.AA(~sum(isoutlier(YzP2,'mean'),2) > 0,:) = S';

% Statistical test for archetypes
for i = 1:size(XsP1,2)
    for ii = 1:size(YsP1.AA,2)
        mdl = fitglm([XsP1{:,i},CovP1],YsP1.AA(:,ii),'CategoricalVars',...
            logical([0 0 0 1]));
        p(i,ii) = mdl.Coefficients.pValue(2);
        e(i,ii) = mdl.Coefficients.Estimate(2);
    end    
end

[p, idx] = sort(p);
SNPs = XsP1.Properties.VariableNames(idx);

%% Case study:
% ins  = find(INS == 1);
% con  = find(INS == 0);
% inss = ins(randi(length(ins),2,1));
% cons = con(randi(length(con),2,1));
% Genetic profile along with sleep parameters
figure
hb = bar(sum(RISKSP2{[inss,cons],:},2)/(186*2));
xticklabels({'Insomniac 1','Insomniac 2','Control 1','Control 2'})
xt = get(gca, 'XTick');
set(gca, 'FontSize', 16)
xtickangle(90)
ylabel('Polygenic risk score')

% Archetypes
hb = bar(S(:,[inss,cons]));
legend({'Insomniac 1','Insomniac 2','Control 1','Control 2'},'Location','northeastoutside');
xlabel('Archetype')
xticklabels({'GS','Low TST, SE','High ARI, AHI', ...
    'Long arousals','LS, low SE, TST','Long SOL','High PLMI'})
xt = get(gca, 'XTick');
set(gca, 'FontSize', 16)
xtickangle(90)

% Model predictions taken from PostProcessing.m: 
%  [ 0.66  0.47  0.63  0.45 ]  \\
%  [ 0.47  0.42  0.83  0.92 ]  \\
%  [ 0.21  0.53  0.37  0.52 ]  \\
%  [ 0.51  0.45  0.61  0.66 ]  \\


