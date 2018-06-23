startup
%% Load
Xorg   = readtable('GenotypesV1.csv'); 
Yorg   = readtable('PhenotypesV1.csv'); 
Yorg.Cohorts(strcmp(Yorg.Cohort,'WSC'),1) = 1;
Yorg.Cohorts(strcmp(Yorg.Cohort,'MrOS'),1) = 2;
Yorg.Cohorts(strcmp(Yorg.Cohort,'SSC'),1) = 3;
Yorg = Yorg(:,3:end);
Xorg = Xorg(:,3:end);

Z   = [Xorg,Yorg];
Z   = rmmissing(Z);

Z.TST = Z.TST/60;
Z.SOL = Z.SOL/60;
Z.REML = Z.REML/60;

cohortGrpZ = Yorg.Cohorts;
X   = Z(:,1:239);
Y   = Z(:,240:end);
Xa  = table2array(X);
Xz  = zscore(Xa);
Ya  = table2array(Y(:,1:end-3));
Yz  = zscore(Ya);
showplots = false;
keepplots = false;
saveplots = false;
Zz  = Z; % Zz.INS = [];
Zz  = table2array(Zz);
Zz  = zscore(Zz);
%% ANOVA
% [p,~,~] = anova1(Yorg.AGE,Yorg.Cohort);
% [p,~,~] = anova1(Yorg.BMI,Yorg.Cohort);
% [p,~,~] = anova1(Yorg.AHI,Yorg.Cohort);
% [p,~,~] = anova1(Yorg.PLMI,Yorg.Cohort);
% [p,~,~] = anova1(Yorg.ESS,Yorg.Cohort);
% [p,~,~] = anova1(Yorg.AGE,Yorg.Cohort);
% [p,~,~] = anova1(Yorg.TST,Yorg.Cohort);
% [p,~,~] = anova1(Yorg.SOL,Yorg.Cohort);
% [p,~,~] = anova1(Yorg.SE,Yorg.Cohort);
% [p,~,~] = anova1(Yorg.REM,Yorg.Cohort);
% [p,~,~] = anova1(Yorg.REML,Yorg.Cohort);
%%
% MEIS1 = round(Xorg.rs4547518);
% MEIS1 = round(Xorg.rs113851554);

if showplots
figure,
histogram(Zins.TST,'Normalization','probability')
hold on
histogram(Znonins.TST,'Normalization','probability')
figure,
histogram(Zins.ARI,'Normalization','probability')
hold on
histogram(Znonins.ARI,'Normalization','probability')
figure,
histogram(Zins.AW25,'Normalization','probability')
hold on
histogram(Znonins.AW25,'Normalization','probability')
figure,
histogram(Zins.AW5,'Normalization','probability')
hold on
histogram(Znonins.AW5,'Normalization','probability')
figure,
histogram(Zins.ESS,'Normalization','probability')
hold on
histogram(Znonins.ESS,'Normalization','probability')
figure,
histogram(Zins.PLMI,'Normalization','probability')
hold on
histogram(Znonins.PLMI,'Normalization','probability')
figure,
histogram(Zins.SE,'Normalization','probability')
hold on
histogram(Znonins.SE,'Normalization','probability')
end
%% Initial probabilities
isWSC = Z.WSC == 1; nWSC = sum(isWSC);
isMrOS = Z.MrOS == 1; nMrOS = sum(isMrOS);
isSSC = Z.SSC == 1; nSSC = sum(isSSC);
% isINSwsc = Z.INS(isWSC); fINSwsc = sum(isINSwsc)/nWSC;
% isINSmros = Z.INS(isMrOS); fINSmros = sum(isINSmros)/nMrOS;
% isINSssc = Z.INS(isSSC); fINSssc = sum(isINSssc)/nSSC;

nMalWSC = sum(Z.SEX(isWSC)); nFemWSC = nWSC - nMalWSC;
nMalMrOS = sum(Z.SEX(isMrOS)); nFemMrOS = nMrOS - nMalMrOS;
nMalSSC = sum(Z.SEX(isSSC)); nFemSSC = nSSC - nMalSSC;

muAgeWSC = mean(Z.AGE(isWSC)); stdAgeWSC = std(Z.AGE(isWSC));
muAgeMrOS = mean(Z.AGE(isMrOS)); stdAgeMrOS = std(Z.AGE(isMrOS));
muAgeSSC = mean(Z.AGE(isSSC)); stdAgeSSC = std(Z.AGE(isSSC));

muBMIWSC = mean(Z.BMI(isWSC)); stdBMIWSC = std(Z.BMI(isWSC));
muBMIMrOS = mean(Z.BMI(isMrOS)); stdBMIMrOS = std(Z.BMI(isMrOS));
muBMISSC = mean(Z.BMI(isSSC)); stdBMISSC = std(Z.BMI(isSSC));

muAHIWSC = mean(Z.AHI(isWSC)); stdAHIWSC = std(Z.AHI(isWSC));
muAHIMrOS = mean(Z.AHI(isMrOS)); stdAHIMrOS = std(Z.AHI(isMrOS));
muAHISSC = mean(Z.AHI(isSSC)); stdAHISSC = std(Z.AHI(isSSC));

muARIWSC = mean(Z.ARI(isWSC)); stdARIWSC = std(Z.ARI(isWSC));
muARIMrOS = mean(Z.ARI(isMrOS)); stdARIMrOS = std(Z.ARI(isMrOS));
muARISSC = mean(Z.ARI(isSSC)); stdARISSC = std(Z.ARI(isSSC));

muPLMIWSC = mean(Z.PLMI(isWSC)); stdPLMIWSC = std(Z.PLMI(isWSC));
muPLMIMrOS = mean(Z.PLMI(isMrOS)); stdPLMIMrOS = std(Z.PLMI(isMrOS));
muPLMISSC = mean(Z.PLMI(isSSC)); stdPLMISSC = std(Z.PLMI(isSSC));

i = 1;
muPSG(i,1) = mean(Z.TST(isWSC)); stdPSG(i,1) = std(Z.TST(isWSC));
muPSG(i,2) = mean(Z.TST(isMrOS)); stdPSG(i,2) = std(Z.TST(isMrOS));
muPSG(i,3) = mean(Z.TST(isSSC)); stdPSG(i,3) = std(Z.TST(isSSC));
i = i + 1;
muPSG(i,1) = mean(Z.SOL(isWSC)); stdPSG(i,1) = std(Z.SOL(isWSC));
muPSG(i,2) = mean(Z.SOL(isMrOS)); stdPSG(i,2) = std(Z.SOL(isMrOS));
muPSG(i,3) = mean(Z.SOL(isSSC)); stdPSG(i,3) = std(Z.SOL(isSSC));
i = i + 1;
muPSG(i,1) = mean(100*Z.SE(isWSC)); stdPSG(i,1) = std(100*Z.SE(isWSC));
muPSG(i,2) = mean(100*Z.SE(isMrOS)); stdPSG(i,2) = std(100*Z.SE(isMrOS));
muPSG(i,3) = mean(100*Z.SE(isSSC)); stdPSG(i,3) = std(100*Z.SE(isSSC));
i = i + 1;
muPSG(i,1) = mean(100*Z.REM(isWSC)); stdPSG(i,1) = std(100*Z.REM(isWSC));
muPSG(i,2) = mean(100*Z.REM(isMrOS)); stdPSG(i,2) = std(100*Z.REM(isMrOS));
muPSG(i,3) = mean(100*Z.REM(isSSC)); stdPSG(i,3) = std(100*Z.REM(isSSC));
i = i + 1;
muPSG(i,1) = mean(Z.REML(isWSC)); stdPSG(i,1) = std(Z.REML(isWSC));
muPSG(i,2) = mean(Z.REML(isMrOS)); stdPSG(i,2) = std(Z.REML(isMrOS));
muPSG(i,3) = mean(Z.REML(isSSC)); stdPSG(i,3) = std(Z.REML(isSSC));
i = i + 1;
muPSG(i,1) = mean(Z.AW25(isWSC)); stdPSG(i,1) = std(Z.AW25(isWSC));
muPSG(i,2) = mean(Z.AW25(isMrOS)); stdPSG(i,2) = std(Z.AW25(isMrOS));
muPSG(i,3) = mean(Z.AW25(isSSC)); stdPSG(i,3) = std(Z.AW25(isSSC));
% i = i + 1;
% muPSG(i,1) = mean(Z.DEP(isWSC)); stdPSG(i,1) = std(Z.DEP(isWSC));
% muPSG(i,2) = mean(Z.DEP(isMrOS)); stdPSG(i,2) = std(Z.DEP(isMrOS));
% muPSG(i,3) = mean(Z.DEP(isSSC)); stdPSG(i,3) = std(Z.DEP(isSSC));
% i = i + 1;
% muPSG(i,1) = mean(Z.ANX(isWSC)); stdPSG(i,1) = std(Z.ANX(isWSC));
% muPSG(i,2) = mean(Z.ANX(isMrOS)); stdPSG(i,2) = std(Z.ANX(isMrOS));
% muPSG(i,3) = mean(Z.ANX(isSSC)); stdPSG(i,3) = std(Z.ANX(isSSC));
% i = i + 1;
% muPSG(i,1) = mean(Z.INS(isWSC)); stdPSG(i,1) = std(Z.INS(isWSC));
% muPSG(i,2) = mean(Z.INS(isMrOS)); stdPSG(i,2) = std(Z.INS(isMrOS));
% muPSG(i,3) = mean(Z.INS(isSSC)); stdPSG(i,3) = std(Z.INS(isSSC));

muPSG = round(muPSG,2);
stdPSG = round(stdPSG,2);
%% PCA and tSNE
[coeffX,scoreX,latentX] = pca(Xz);
[coeffY,scoreY,latentY] = pca(Yz);
if showplots
figure,% subplot(3,2,1)
    g = gscatter(scoreX(:,1),scoreX(:,2),cohortGrpZ,'brg','xos');  
    xlabel('PC 1'), grid on
    ylabel('PC 2'), box off
    hLeg = legend({'WSC','MrOS','SSC'}); % set(hLeg,'visible','off');
    set(gca,'color','none')
export_fig 'PCA_SNPs.png' -transparent 
if saveplots,print(strcat(plotPath2,'PCA_SNPs'),'-depsc'); end
figure,%subplot(3,2,3)
    biplot(coeffX(:,1:2),'VarLabels',X.Properties.VariableNames), xlabel('PC 1'), 
    ylabel('PC 2'), 
if saveplots,print(strcat(plotPath2,'Biplot_SNPs'),'-depsc'); end
figure,%subplot(3,2,5)
    stem(100*latentX/(sum(latentX))), xlabel('Components'), 
    ylabel('Explained Variance (\%)'), box off, grid on, xlim([0 length(latentX)]),
if saveplots,print(strcat(plotPath2,'Explained_Variance_SNPs'),'-depsc'); end
figure,%subplot(3,2,2)
    g = gscatter(scoreY(:,1),scoreY(:,2),cohortGrpZ,'brg','xos'); 
    xlabel('PC 1'), grid on
    ylabel('PC 2'), box off
    hLeg = legend({'WSC','MrOS','SSC'}); % set(hLeg,'visible','off');
if saveplots,print(strcat(plotPath2,'PCA_PSG'),'-depsc'); end
figure,%subplot(3,2,4)
    biplot(coeffY(:,1:2),'VarLabels',Y.Properties.VariableNames), 
    xlabel('PC 1'), 
    ylabel('PC 2'), 
if saveplots,print(strcat(plotPath2,'Biplot_PSG'),'-depsc'); end
figure,%subplot(3,2,6)
    stem(100*latentY/(sum(latentY))), xlabel('Components'), 
    ylabel('Explained Variance (\%)'), box off, xlim([0 length(latentY)]), grid on
if saveplots,print(strcat(plotPath2,'Explained_Variance_PSG'),'-depsc'); end
figure,
hSub = subplot(5,1,1); gscatter([1, 1],[nan, nan],[0,1], 'br','xo'); set(hSub, 'Visible', 'off'); 
subplot(5,1,[2 3]), g = gscatter(scoreX(:,1),scoreX(:,2),isInsomniac,'br','xo'); 
ylabel('PC 2','FontSize',16), grid on, box off, %title('PCA of allel dosages')
hLeg = legend({'Non-Insomnia','Insomnia'}); set(hLeg,'visible','off');
xt = get(gca, 'XTick'); set(gca, 'FontSize', 16);
set(gca,'color','none')
subplot(5,1,[4 5]), g = gscatter(scoreY(:,1),scoreY(:,2),isInsomniac,'br','xo'); 
xlabel('PC 1','FontSize',16), grid on, %title('PCA of PSG variables')
ylabel('PC 2','FontSize',16), box off,
hLeg = legend({'Non-Insomnia','Insomnia'}); set(hLeg,'visible','off');
legend(hSub, 'Non-Insomnia','Insomnia', 'Location', 'east');
set(gca,'color','none')
export_fig 'PCA_INS.png' -transparent 
if saveplots,print(strcat(plotPath2,'PCA_INS'),'-depsc'); end
end
if ~keepplots, close all; end
%% Outliers
TFsnps = isoutlier(scoreX,'mean'); TFsnps = ~TFsnps(:,1);
TFpsg = isoutlier(scoreY,'mean'); TFpsg = ~TFpsg(:,1);
if showplots
figure,gscatter(scoreX(:,1),scoreX(:,2),~TFsnps,'br','xx'); 
hLeg = legend('Included','Outlier'); 
xlabel('PC 1'), grid on, 
ylabel('PC 2'), box off,
set(gca,'color','none')
figure,gscatter(scoreY(:,1),scoreY(:,2),~TFpsg,'br','xx'); 
hLeg = legend('Included','Outlier'); 
xlabel('PC 1'), grid on, 
ylabel('PC 2'), box off,
set(gca,'color','none')
end
Xnew = X(TFsnps,:); Ynew = Y(TFsnps,:);
Yz = Yz(TFsnps,:); Xz = Xz(TFsnps,:);
scoreX = scoreX(TFsnps,:);
scoreY = scoreY(TFsnps,:);
if showplots
    figure,
    hSub = subplot(5,1,1); gscatter([1, 1],[nan, nan],[0,1], 'br','xo'); set(hSub, 'Visible', 'off'); 
    subplot(5,1,[2 3]), g = gscatter(scoreX(:,1),scoreX(:,2),isInsomniac,'br','xo'); 
    ylabel('PC 2'), grid on, box off, %title('PCA of allel dosages')
    hLeg = legend({'Non-Insomnia','Insomnia'}); set(hLeg,'visible','off');
    xt = get(gca, 'XTick'); set(gca, 'FontSize', 16);
    set(gca,'color','none')
    subplot(5,1,[4 5]), g = gscatter(scoreY(:,1),scoreY(:,2),isInsomniac,'br','xo'); 
    xlabel('PC 1'), grid on, %title('PCA of PSG variables')
    ylabel('PC 2'), box off,
    hLeg = legend({'Non-Insomnia','Insomnia'}); set(hLeg,'visible','off');
    legend(hSub, 'Non-Insomnia','Insomnia', 'Location', 'east');
    set(gca,'color','none')
    export_fig 'PCA_INS.png' -transparent 
end


% for k = 1:20
%     [Train, Test] = crossvalind('HoldOut', size(Yz,1), 0.25);
%     for c = 1:20
%         GMModel = fitgmdist(Yz(Train,:),2,'CovarianceType','diagonal');
%         trainNLog(k,c) = GMModel.NegativeLogLikelihood;
%         [~, testNLog(k,c)] = posterior(GMModel,Yz(Test,:));
%     end
% end
% a = corrcoef(Xa);
% imagesc(abs(a)>0.6), colorbar
K = 6;
figure,
yyaxis left
plot(median(trainNLog));
yyaxis right
plot(median(testNLog));
GMModel = fitgmdist(scoreX,6,'CovarianceType','diagonal');
[P, ~] = posterior(GMModel,scoreX);
[~,label] = max(P,[],2);
histogram(label)

[label,b] = kmeans(Xa(TFsnps,:),6);
histogram(label)
gscatter(scoreX(:,1),scoreX(:,2),label,'brgym','osxox')
set(gca,'color','none')


import bioma.data.*

alpha = 0.01;
for ii = 1:6
    A = DataMatrix(Xa(label == ii,:)',Xorg.Properties.VariableNames,[]);
    idx = 1;
    for i = 1:6, if i == ii, continue, end
        B = DataMatrix(Xa(label == i,:)',Xorg.Properties.VariableNames,[]);
        PValues = mattest(A, B);
        ps{ii}(:,idx) = mafdr(PValues.double)<alpha;
        idx = idx + 1;
    end
end


c1 = [find(prod(ps{1},2) == 1)';...
    median(Xa(label == 1,prod(ps{1},2) == 1))];
X.Properties.VariableNames(c1(1,:))
% c2 = [find(prod(ps{2},2) == 1)';...
%     median(Xa(label == 2,prod(ps{3},2) == 1))];
c3 = [find(prod(ps{3},2) == 1)';...
    median(Xa(label == 3,prod(ps{3},2) == 1))];
X.Properties.VariableNames(c3(1,:))
c4 = [find(prod(ps{4},2) == 1)';...
    median(Xa(label == 4,prod(ps{4},2) == 1))];
X.Properties.VariableNames(c4(1,:))
c5 = [find(prod(ps{5},2) == 1)';...
    median(Xa(label == 5,prod(ps{5},2) == 1))];
c6 = [find(prod(ps{6},2) == 1)';...
    median(Xa(label == 6,prod(ps{6},2) == 1))];

histogram(label)

Ynew.Label = label;
Ynew.Label = [];

% if showplots
% try load PSG_GMM, catch
% GMModel = fitgmdist(scoreY,2,'CovarianceType','diagonal',...
%             'Replicates',200);
% save('PSG_GMM','GMModel')
% end
% 
% [~,label] = max(posterior(GMModel,scoreY),[],2);
% X.Label = label;
% X.Label = [];
% figure,gscatter(scoreY(:,1),scoreY(:,2),label,'br','xo'); 
% ylabel('PC 2','FontSize',16), grid on, box off
% xlabel('PC 1','FontSize',16), 
% hLeg = legend({'Cluster 1','Cluster 2','Cluster 3'},'Location','northwest'); 
% set(gca,'color','none')
% export_fig 'PCA_PSG_GMM.png' -transparent 
% figure,
% gscatter(scoreY(:,1),scoreY(:,2),isInsomniac,'br','xo'); 
% ylabel('PC 2','FontSize',16), grid on, box off
% xlabel('PC 1','FontSize',16), 
% legend({ 'Non-Insomnia','Insomnia'}, 'Location', 'northwest');
% set(gca,'color','none')
% export_fig 'PCA_PSG_INS.png' -transparent 
% end
%% tSNE
s = tsne(Xz);
figure
    g = gscatter(s(:,1),s(:,2),cohortGrpZ,'brg','xos'); grid on
    legend({'WSC','MrOS','SSC'})
s = tsne(Yz);
figure
    g = gscatter(s(:,1),s(:,2),cohortGrpZ,'brg','xos'); grid on
    legend({'WSC','MrOS','SSC'})
%% Not including covariates
X       = Xorg;
Y       = Yorg;
mdls    = cell(1,size(Y,2));
pVals   = ones(size(X,2),size(Y,2));
effects = ones(size(X,2),size(Y,2));
SNPs    = cell(size(X,2),size(Y,2));
allPredictors = X.Properties.VariableNames;
for i = 1:size(Y,2)
    mdls{i}         = fitglm([X,Y(:,i)]);
    [pVal,idx]      = sort(mdls{i}.Coefficients.pValue(2:end));
    effect          = mdls{i}.Coefficients.Estimate(2:end);
    pVals(:,i)      = pVal;
    effects(:,i)    = effect(idx);
    SNPs(:,i)       = allPredictors(idx);
end
% PLMI: rs113851554
%% Alpha correction
% Alpha
alpha = 0.05;
% Bonferroni
alpha_bonf = alpha/(length(mdls)*size(X,2));
% Benjamini and Hochberg 
k = size(Xorg,2)*size(Yorg,2);
ordpVals = sort(pVals(:),'descend');
alpha_bhIdx = find(ordpVals<=((1:k)*(alpha/k))');
alpha_bh = ordpVals(alpha_bhIdx(1));
%  Benjamini and Yekutieli
alpha_by = alpha/sum(1./(1:k));

SNPs(pVals<alpha_by);
% pVals(pVals<alpha_by))

%% Including covariates
X       = Xnew;
Y       = Ynew;
% X.PLMI = Y.PLMI; Y.PLMI = [];
X.AGE = Y.AGE; Y.AGE = [];
% X.AHI = Y.AHI; Y.AHI = [];
% X.ARI = Y.ARI; Y.ARI = [];
X.BMI = Y.BMI; Y.BMI = [];
% X.DEP = Y.DEP; Y.DEP = [];
% X.ANX = Y.ANX; Y.ANX = [];
X.SEX = Y.SEX; Y.SEX = [];
X.WSC = Y.WSC; Y.WSC = [];
X.MrOs = Y.MrOS; Y.MrOS = [];
X.SSC = Y.SSC; Y.SSC = [];
mdls    = cell(1,size(Y,2));
pVals   = ones(size(X,2),size(Y,2));
effects = ones(size(X,2),size(Y,2));
SNPs    = cell(size(X,2),size(Y,2));
allPredictors = X.Properties.VariableNames;
for i = 1:size(Y,2)
    mdls{i}         = fitglm([X,Y(:,i)]);
    [pVal,idx]      = sort(mdls{i}.Coefficients.pValue(2:end));
    effect          = mdls{i}.Coefficients.Estimate(2:end);
    pVals(:,i)      = pVal;
    effects(:,i)    = effect(idx);
    SNPs(:,i)       = allPredictors(idx);
end
%% Per cohort
mdls    = cell(3,1);
pVals   = cell(3,1);
effects = cell(3,1);
SNPs    = cell(3,1);
allPredictors = X.Properties.VariableNames;
for c = 0:2 % cohorts
    idxs = X{:,end-c} > 0;
for i = 1:size(Y,2)
    mdls{c+1,i}         = fitglm([X(idxs,1:end-3),Y(idxs,i)]);
    [pVal,idx]          = sort(mdls{c+1,i}.Coefficients.pValue(2:end));
    effect              = mdls{c+1,i}.Coefficients.Estimate(2:end);
    pVals{c+1}(:,i)      = pVal;
    effects{c+1}(:,i)    = effect(idx);
    SNPs{c+1}(:,i)       = allPredictors(idx);
end
end
%% Logistic regression covariates
% X       = [Xorg,Yorg]; X.INS = [];
% Y       = Yorg(:,12);
% B = fitglm([X, Y], 'Distribution', 'binomial', 'Link', 'logit');


[beta,Sigma,E,CovB,logL] = mvregress(table2cell(Xorg),table2array(Yorg));

X       = Xorg;
Y       = Yorg;
X.PLMI = Y.PLMI; Y.PLMI = [];
X.AGE = Yorg.AGE; Y.AGE = [];
X.AHI = Yorg.AHI; Y.AHI = [];
X.ARI = Y.ARI; Y.ARI = [];
X.BMI = Yorg.BMI; Y.BMI = [];
X.DEP = Yorg.DEP; Y.DEP = [];
X.SEX = Yorg.SEX; Y.SEX = [];
X.WSC = Y.WSC; Y.WSC = [];
X.MrOs = Y.MrOS; Y.MrOS = [];
X.SSC = Y.SSC; Y.SSC = [];
[beta,Sigma,E,CovB,logL] = mvregress(table2array(X),table2array(Y));
%% Cross validation
% kFolds  = 20;
mdls    = cell(3,size(Y,2));
pVals   = cell(3,1);
effects = cell(3,1);
SNPs    = cell(3,1);
allPredictors = X.Properties.VariableNames;
for c = 0:2 % cohorts
    idxs = Y{:,end-c}>0;
for i = 1:size(Y,2)
    mdls{c+1,i}         = fitglm([X(idxs,:),Y(idxs,i)]);
    [pVal,idx]      = sort(mdls{c+1,i}.Coefficients.pValue(2:end));
    effect          = mdls{c+1,i}.Coefficients.Estimate(2:end);
    pVals{c+1}(:,i)      = pVal;
    effects{c+1}(:,i)    = effect(idx);
    SNPs{c+1}(:,i)       = allPredictors(idx);
end
end
% FDR = mafdr(pVals{1,1}(:,1))
% N       = 3000;
% for k = 1:kFolds
%     Indices = crossvalind('Kfold', N, size(Yall,1));
%     pVals{k}   = nan(size(Xall,2),size(Yall,2));
%     effects{k} = nan(size(Xall,2),size(Yall,2));
% for i = 1:size(Yall,2)
%     mdls{k,i}         = fitglm([Xall(Indices,:),Yall(Indices,i)]);
%     pVals{k}(:,i)   = mdls{k,i}.Coefficients.pValue(2:end);
%     effects{k}(:,i) = mdls{k,i}.Coefficients.Estimate(2:end);
% end
% disp(k)
% end


%% Hypnogram 
nHypFeatures    = 9;
snps            = cell(size(pVals,1),nHypFeatures);
phenotypes      = strcat('',strrep(strrep(curFields,'h_',''),'aw','Aw'),'');
A               = cell(size(pVals,1),nHypFeatures); A(:)={''};
maxHeight       = 1;
for i = 1:nHypFeatures
    snps(1:length(SNPs(pVals(:,i)<alpha_bh,i)),i) = SNPs(pVals(:,i)<alpha_bh,i);
    A(1:length(SNPs(pVals(:,i)<alpha_bh,i)),i) = snps(1:length(SNPs(pVals(:,i)<alpha_bh,i)),i); 
    if sum(pVals(:,i)<alpha_by)>0,A(pVals(:,i)<alpha_by,i)= strcat(A(pVals(:,i)<alpha_by,i),{'*'}); end
    if sum(pVals(:,i)<alpha_bonf)>0,A(pVals(:,i)<alpha_bonf,i)= strcat(A(pVals(:,i)<alpha_bonf,i),{'*'}); end
    if length(SNPs(pVals(:,i)<alpha_bh,i))>maxHeight,maxHeight=length(SNPs(pVals(:,i)<alpha_bh,i)); end
end
A = A(1:maxHeight,:);
AA = strrep(A,'*','');
ismember(AA,minAllelFreqSNP)
matrix2latex(A,'Data/Hyp_Significant_SNPs.tex',...
    'columnLabels',strrep(phenotypes(1:nHypFeatures),'_',' '), ...
    'longtable', true, 'lines', false);
%% EEG 
nEEGFeatures    = length(curFields) - nHypFeatures;
snps            = cell(size(pVals,1),nEEGFeatures);
phenotypes      = strcat('',strrep(strrep(curFields,'C3Ref_',''),'mean',''),'');
phenotypes      = strrep(strrep(phenotypes,'C3Ref',''),'mu','all');
curNumSnps      = 0;
A               = cell(size(pVals,1),nEEGFeatures); A(:)={''};
maxHeight       = 1;
for i = nHypFeatures+1:nHypFeatures+nEEGFeatures
    idx = i - nHypFeatures;
    snps(1:length(SNPs(pVals(:,i)<alpha_bh,i)),idx) = SNPs(pVals(:,i)<alpha_bh,i);
    A(1:length(SNPs(pVals(:,i)<alpha_bh,i)),idx) = snps(1:length(SNPs(pVals(:,i)<alpha_bh,i)),idx); 
    if sum(pVals(:,i)<alpha_by)>0,A(pVals(:,i)<alpha_by,idx)= strcat(A(pVals(:,i)<alpha_by,idx),{'*'}); end
    if sum(pVals(:,i)<alpha_bonf)>0,A(pVals(:,i)<alpha_bonf,idx)= strcat(A(pVals(:,i)<alpha_bonf,idx),{'*'}); end
    if length(SNPs(pVals(:,i)<alpha_bh,i))>maxHeight,maxHeight=length(SNPs(pVals(:,i)<alpha_bh,i)); end
end
A = A(1:maxHeight,:);
AA = strrep(A,'*','');
ismember(AA,minAllelFreqSNP)
matrix2latex(A,'Data/EEG_Significant_SNPs.tex',...
    'columnLabels',strrep(phenotypes(nHypFeatures+1:nHypFeatures+nEEGFeatures),'_',' '), ...
    'longtable', true, 'lines', false);
%%
significance = struct();
T = readtable('SNP2Gene.csv');
snpss = [];
for i = 1:length(mdls)
    pVals = mdls{i}.Coefficients.pValue(2:end);
    snps = mdls{i}.Coefficients.Properties.RowNames(2:end);
    snpss = [snpss; snps(pVals<alpha_bh)];
    f = strrep(mdls{i}.ResponseName,'h_','');
    significance.(f).pBH = pVals(pVals<alpha_bh);
    significance.(f).snpBH = snps(pVals<alpha_bh);
    significance.(f).eBH = mdls{i}.Coefficients.Estimate(pVals<alpha_bh);
    significance.(f).pBY = pVals(pVals<alpha_by);
    significance.(f).snpBY = snps(pVals<alpha_by);
%     significance.(f).eBY = mdls_hyp{i}.Coefficients.Estimate(pVals<alpha_by);
end
featureNames = fields(significance);
columnNames = [phenotypes; 'Nearest Gene'];
rowNames = unique(snps);
genes = T.NearestGene(contains(T.rsID,rowNames));
A = cell(length(rowNames),length(columnNames)-1);
A(:) = {''};
for i = 1:length(featureNames)
    curSig = significance.(featureNames{i});
    pBH = curSig.pBH;
    eBH = curSig.eBH;
    for ii = 1:length(pBH)
        curSNP = curSig.snpBH{ii};
        if eBH(ii) < 0
            A{contains(rowNames,curSNP),i} = sprintf('$%.2f$*',eBH(ii));
        else
            A{contains(rowNames,curSNP),i} = sprintf('$+%.2f$*',eBH(ii));
        end
        if sum(contains(curSig.snpBY,curSNP)) == 1 
            if eBH(ii) < 0        
                A{contains(rowNames,curSNP),i} = sprintf('$%.2f$**',eBH(ii));
            else
                A{contains(rowNames,curSNP),i} = sprintf('$+%.2f$**',eBH(ii));
            end
        end
    end
end
A = [A, genes];
matrix2latex(A,'Data/WSC/All_Significant_Effects.tex',...
    'rowLabels',rowNames,'columnLabels',columnNames, ...
    'longtable', true, 'lines', false);

matrix2latex(SNPs(1:5,15:end),'Data/Test.tex',...
    'longtable', true, 'lines', false);

