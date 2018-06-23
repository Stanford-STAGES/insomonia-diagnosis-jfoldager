startup
%% Load
load V2newdata.mat
shouldplot = false;
Xa      = table2array(X);
% Eliminate large cross correlations
Y.AW5   = [];
Y.AWN15 = [];
Y.WA    = [];
Ya      = table2array(Y(:,1:end-1));
% Load SNP to gene map
snp2gene  = readtable('SNP2GeneNew.csv');
%% PCA plot
[coeffX,scoreX,latentX] = pca(zscore(Xa));
[coeffY,scoreY,latentY] = pca(zscore(Ya));
figure,% subplot(3,2,1)
    g = gscatter(scoreX(:,1),scoreX(:,2),Y.Cohorts,'brg','xos');  
    xlabel('PC 1'), grid on
    ylabel('PC 2'), box off
    hLeg = legend({'WSC','MrOS','SSC'}); % set(hLeg,'visible','off');
    set(gca,'color','none')
print(strcat(plotPath2,'PCA_SNPs'),'-depsc');
figure,%subplot(3,2,3)
    biplot(coeffX(:,1:2),'VarLabels',X.Properties.VariableNames), xlabel('PC 1'), 
    ylabel('PC 2'), 
print(strcat(plotPath2,'Biplot_SNPs'),'-depsc');
figure,%subplot(3,2,5)
    stem(100*latentX/(sum(latentX))), xlabel('Components'), 
    ylabel('Explained Variance (\%)'), box off, grid on, xlim([0 length(latentX)]),
print(strcat(plotPath2,'Explained_Variance_SNPs'),'-depsc');
figure,%subplot(3,2,2)
    g = gscatter(scoreY(:,1),scoreY(:,2),Y.Cohorts,'brg','xos'); 
    xlabel('PC 1'), grid on
    ylabel('PC 2'), box off
    hLeg = legend({'WSC','MrOS','SSC'}); % set(hLeg,'visible','off');
print(strcat(plotPath2,'PCA_PSG'),'-depsc');
figure,%subplot(3,2,4)
    biplot(coeffY(:,1:2),'VarLabels',Y.Properties.VariableNames), 
    xlabel('PC 1'), 
    ylabel('PC 2'), 
print(strcat(plotPath2,'Biplot_PSG'),'-depsc'); 
figure,%subplot(3,2,6)
    stem(100*latentY/(sum(latentY))), xlabel('Components'), 
    ylabel('Explained Variance (\%)'), box off, xlim([0 length(latentY)]), grid on
print(strcat(plotPath2,'Explained_Variance_PSG'),'-depsc'); 
figure,
%     hSub = subplot(5,1,1); 
%         gscatter([1, 1],[nan, nan],[0,1], 'br','xo'); set(hSub, 'Visible', 'off'); 
%     subplot(5,1,[2 3]), 
figure,
    g = gscatter(scoreX(:,1),scoreX(:,2),Y.INS,'br','xo'); 
    ylabel('PC 2','FontSize',16), grid on, box off, %title('PCA of allel dosages')
    hLeg = legend({'Non-Insomnia','Insomnia'}); % set(hLeg,'visible','off');
    set(gca,'color','none')
    subplot(5,1,[4 5]), g = gscatter(scoreY(:,1),scoreY(:,2),Y.INS,'br','xo'); 
        xlabel('PC 1'), grid on, %title('PCA of PSG variables')
        ylabel('PC 2'), box off,
        hLeg = legend({'Non-Insomnia','Insomnia'}); set(hLeg,'visible','off');
        legend(hSub, 'Non-Insomnia','Insomnia', 'Location', 'east');
        set(gca,'color','none')
print(strcat(plotPath2,'PCA_INS'),'-depsc'); 
%% Cross correlation
% [rhoY,pvalY] = corrcoef(Ya);
rhoY = eye(size(Ya,2));
cats = {'ANX','DEP','INS','SEX'};
for i = 1:size(Y,2)-1
   for ii = 1:(size(Y,2)-1), if i == ii, continue; end
       if sum(contains(Y.Properties.VariableNames{i},cats)) == 1
           if sum(contains(Y.Properties.VariableNames{ii},cats)) == 1
                rho = corr([Ya(:,i),Ya(:,ii)],'Type','Spearman');
                rhoY(i,ii) = rho(2,1);
           else
                [rhoY(i,ii),~,~,~] = pointbiserial(Ya(:,i),Ya(:,ii));
           end 
       elseif sum(contains(Y.Properties.VariableNames{ii},cats)) == 1
           if sum(contains(Y.Properties.VariableNames{i},cats)) == 1
                rho = corr([Ya(:,i),Ya(:,ii)],'Type','Spearman');
                rhoY(i,ii) = rho(2,1);
           else
                [rhoY(i,ii),~,~,~] = pointbiserial(Ya(:,ii),Ya(:,i));
           end 
       else
           a = corrcoef(Ya(:,i),Ya(:,ii));
           rhoY(i,ii) =  a(2,1);
       end
   end
end
%     Plot cross correlations
figure,imagesc(rhoY)
ax = gca; 
names = Y.Properties.VariableNames;
ax.YTick = 1:length(rhoY); ax.YTickLabel = names;
ax.XTick = 1:length(rhoY); ax.XTickLabel = names;
xtickangle(ax,90), c = colorbar; 
print(strcat(plotPath3,'PSG_CrossCor'),'-depsc');
%     Plot cross correlations > 0.85
figure,imagesc(abs(rhoY) > 0.85)
ax = gca; 
names = Y.Properties.VariableNames;
ax.YTick = 1:length(rhoY); ax.YTickLabel = names;
ax.XTick = 1:length(rhoY); ax.XTickLabel = names;
xtickangle(ax,90), c = colorbar; c.Ticks = [0,1]; colormap(gray(2)) 
c.TickLabels = {'<0.85','\geq 0.85'} ; 
print(strcat(plotPath3,'PSG_CrossCor085'),'-depsc');
%%
% Insomnia / non insomnia statistics
idx = 1;
for i = 1:length(Y.Properties.VariableNames)-1
   if strcmp(Y.Properties.VariableNames{i},'INS'), continue, end
    [~,t(i,1)] = ttest2(Y.(Y.Properties.VariableNames{i})(Y.INS == 1),...
        Y.(Y.Properties.VariableNames{i})(Y.INS == 0));
    uniTables(i,1) = mean(Y.(Y.Properties.VariableNames{i})(Y.INS == 1));
    uniTables(i,2) = std(Y.(Y.Properties.VariableNames{i})(Y.INS == 1));
    uniTables(i,3) = mean(Y.(Y.Properties.VariableNames{i})(Y.INS == 0));
    uniTables(i,4) = std(Y.(Y.Properties.VariableNames{i})(Y.INS == 0));
    if sum(strcmp(Y.Properties.VariableNames{i},{'ANX','DEP','SEX'})) == 1
        [~,p(idx)] = chi2gof(Y.(Y.Properties.VariableNames{i})(Y.INS == 1));
    idx = idx + 1;
    end
end
uniTables(strcmp(Y.Properties.VariableNames,'INS'),:) = [];
t(strcmp(Y.Properties.VariableNames,'INS'),:) = [];
%% LOGISTIC REGRESSION predict insomnia yes/no
load V2newdata.mat
Ys      = Y;
% Ys.N1   = [];
% Ys.N2   = [];
% Ys.N3   = [];
% Ys.REM  = [];
% Ys.REML = [];
resps    = table();
resps.INS = Ys.INS;
Ys.INS      = [];
categ = {'ANX','DEP','SEX','Cohorts'};
clear p beta or
for i = 1:size(Ys,2)
    if sum(strcmp(Ys.Properties.VariableNames(i), categ)) == 1 
        logmdl = fitglm([Ys(:,i),resps],'CategoricalVars',...
            Ys.Properties.VariableNames(i),'Distribution','binomial'); 
    else
        logmdl = fitglm([Ys(:,i),resps],'Distribution','binomial');  
    end
    p(i,1)  = logmdl.Coefficients.pValue(2);
    beta(i,1)  = logmdl.Coefficients.Estimate(2);
    or(i,1)  = exp(beta(i,1));
end
% Next iteration
Covariates = table();
Covariates.DEP = Ys.DEP; Ys.DEP = [];
all = table();
for i = 1:size(Ys,2)
    if sum(strcmp(Ys.Properties.VariableNames(i), categ)) == 1 
        logmdl = fitglm([Ys(:,i),Covariates,resps],'CategoricalVars',...
            [Ys.Properties.VariableNames(i),'DEP'],'Distribution','binomial'); 
    else
        logmdl = fitglm([Ys(:,i),Covariates,resps],'CategoricalVars',...
            'DEP','Distribution','binomial');  
    end
    p(i,1) = logmdl.Coefficients.pValue(2);
    t(i,1) = logmdl.Coefficients.tStat(2);
    se(i,1) = logmdl.Coefficients.SE(2);
    e(i,1) = logmdl.Coefficients.Estimate(2);
    all = [all;logmdl.Coefficients(2,:)];
end
% Next iteration
Covariates.Cohorts = Ys.Cohorts; Ys.Cohorts = [];
all = table();
for i = 1:size(Ys,2)
    if sum(strcmp(Ys.Properties.VariableNames(i), categ)) == 1 
        logmdl = fitglm([Ys(:,i),Covariates,resps],'CategoricalVars',...
            [Ys.Properties.VariableNames(i),'DEP','Cohorts'],'Distribution','binomial'); 
    else
        logmdl = fitglm([Ys(:,i),Covariates,resps],'CategoricalVars',...
            {'DEP','Cohorts'},'Distribution','binomial');  
    end
    p(i,1) = logmdl.Coefficients.pValue(2);
    t(i,1) = logmdl.Coefficients.tStat(2);
    se(i,1) = logmdl.Coefficients.SE(2);
    e(i,1) = logmdl.Coefficients.Estimate(2);
    all = [all;logmdl.Coefficients(2,:)];
end
% Next iteration
Covariates.TST = Ys.TST; Ys.TST = [];
all = table();
for i = 1:size(Ys,2)
    if sum(strcmp(Ys.Properties.VariableNames(i), categ)) == 1 
        logmdl = fitglm([Ys(:,i),Covariates,resps],'CategoricalVars',...
            [Ys.Properties.VariableNames(i),'DEP','Cohorts','TST'],'Distribution','binomial'); 
    else
        logmdl = fitglm([Ys(:,i),Covariates,resps],'CategoricalVars',...
            {'DEP','Cohorts','TST'},'Distribution','binomial');  
    end
    p(i,1) = logmdl.Coefficients.pValue(2);
    t(i,1) = logmdl.Coefficients.tStat(2);
    se(i,1) = logmdl.Coefficients.SE(2);
    e(i,1) = logmdl.Coefficients.Estimate(2);
    all = [all;logmdl.Coefficients(2,:)];
end
Covariates.ANX = Ys.ANX; Ys.ANX = [];
all = table();
for i = 1:size(Ys,2)
    if sum(strcmp(Ys.Properties.VariableNames(i), categ)) == 1 
        logmdl = fitglm([Ys(:,i),Covariates,resps],'CategoricalVars',...
            [Ys.Properties.VariableNames(i),'DEP','Cohorts','TST','ANX'],'Distribution','binomial'); 
    else
        logmdl = fitglm([Ys(:,i),Covariates,resps],'CategoricalVars',...
            {'DEP','Cohorts','TST','ANX'},'Distribution','binomial');  
    end
    p(i,1) = logmdl.Coefficients.pValue(2);
    t(i,1) = logmdl.Coefficients.tStat(2);
    se(i,1) = logmdl.Coefficients.SE(2);
    e(i,1) = logmdl.Coefficients.Estimate(2);
    all = [all;logmdl.Coefficients(2,:)];
end
%% MODEL 1 age, sex, cohort, pca
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%%%%% UNIVARIATE LOGISTIC REGRESSION model predict insomnia yes/no WITH SNPS 
% + covariates
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
startup
categ = {'SEX','Cohorts'};
load V2newdata.mat
% Predictors
Xs      = X;
% Covariates
covars          = table();
covars.AGE      = Y.AGE; 
covars.SEX      = Y.SEX; 
covars.Cohorts  = Y.Cohorts; 
for i = 1:10, covars.(strcat('PC',num2str(i))) = nan(size(covars,1),1); end
PCs             = readtable('PCA1to10.csv');
for i = 1:size(IDs,1)
    if sum(strcmp(PCs.ID,IDs.ID{i})) == 1
        for ii = 1:10, covars.(strcat('PC',num2str(ii)))(i) = ...
                PCs{strcmp(PCs.ID,IDs.ID{i}),ii+1}; end
    end
end
% Response
resps       = table();
resps.INS   = Y.INS; 
% Model
pVals = nan(size(Xs,2),1);
effects = nan(size(Xs,2),1);
for i = 1:size(Xs,2)
    linmdl      = fitglm([Xs(:,i),covars,resps], 'CategoricalVars',...
        categ,'Distribution','binomial');
    pVals(i)    = linmdl.Coefficients.pValue(2);
    effects(i)  = linmdl.Coefficients.Estimate(2);
end
[pValsSort,idxs] = sort(pVals);
effectsSort      = effects(idxs);
predictorsSort   = Xs.Properties.VariableNames(idxs)';
nearestgene      = cell(size(predictorsSort));
snp2gene         = readtable('SNP2GeneNew.csv');
for i = 1:length(predictorsSort)
    nearestgene{i} = snp2gene.NearestGene{ismember(snp2gene.rsID,predictorsSort{i})};
end
pValsSortFDR =  mafdr(pValsSort);
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%%%%%%%% UNIVARIATE LINEAR REGRESSION model predict PSG WITH SNPS + covariates
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
load V1newdata.mat
% Predictors
Xs      = X;
% Responses
Ys      = table();
Ys.ARI  = Y.ARI;
Ys.AW25 = Y.AW25;
Ys.AWN125  = Y.AWN125;
Ys.PLMI = Y.PLMI;
Ys.SE   = Y.SE;
Ys.SOL  = Y.SOL;
Ys.TST  = Y.TST;
% Covariates
covars          = table();
covars.AGE      = Y.AGE; 
covars.SEX      = Y.SEX; 
covars.Cohorts  = Y.Cohorts; 
categ = {'SEX','Cohorts'};
for i = 1:10, covars.(strcat('PC',num2str(i))) = nan(size(covars,1),1); end
PCs             = readtable('PCA1to10.csv');
for i = 1:size(IDs,1)
    if sum(strcmp(PCs.ID,IDs.ID{i})) == 1
        for ii = 1:10, covars.(strcat('PC',num2str(ii)))(i) = ...
                PCs{strcmp(PCs.ID,IDs.ID{i}),ii+1}; end
    end
end
pVals = nan(size(Xs,2),size(Ys,2));
effects = nan(size(Xs,2),size(Ys,2));
alpha_bonf = 0.05/(size(Xs,2)*size(Ys,2));
predictors = repmat(Xs.Properties.VariableNames',1,size(Ys,2));
% Predict each Ys with each snp
for i = 1:size(Xs,2)
    for ii = 1:size(Ys,2) 
        linmdl          = fitglm([Xs(:,i),covars,Ys(:,ii)],...
            'CategoricalVars',categ);
        pVals(i,ii)     = linmdl.Coefficients.pValue(2);
        effects(i,ii)   = linmdl.Coefficients.Estimate(2);
    end
end
[pValsSort,idxs]= sort(pVals);
effectsSort     = effects(idxs);
predictorsSort  = predictors(idxs);
nearestgene      = cell(size(predictorsSort));
% Load SNP to gene map
snp2gene         = readtable('SNP2GeneNew.csv');
for i = 1:size(predictorsSort,1)
    for ii = 1:size(predictorsSort,2)
        nearestgene{i,ii} = snp2gene.NearestGene{ismember(snp2gene.rsID,predictorsSort{i,ii})};
    end
end
uniTables = {};
for i = 1:size(predictorsSort,2)
    uniTables{i} = table();
    uniTables{i}.rsID = predictorsSort(:,i);
    uniTables{i}.NearestGene = nearestgene(:,i);
    uniTables{i}.Effect = num2str(round(effectsSort(:,i),4));
    uniTables{i}.pValue = num2str(round(pValsSort(:,i),4));
end
[corrected] = BHcorrection(pValsSort, 0.1);
%% MODEL 2 age, anxiety, depression, sex, cohort, pca
% startup
% Load SNP to gene map
snp2gene         = readtable('SNP2GeneNew.csv');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%%%%% UNIVARIATE LOGISTIC REGRESSION model predict insomnia yes/no WITH SNPS 
% + covariates
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
load V2newdata.mat
% Predictors
Xs              = X;
% Covariates
covars          = table();
covars.AGE      = Y.AGE; 
covars.SEX      = Y.SEX; 
covars.DEP      = Y.DEP; 
covars.ANX      = Y.ANX; 
categ = {'SEX','Cohorts','DEP','ANX'};
covars.Cohorts  = Y.Cohorts; 
for i = 1:10, covars.(strcat('PC',num2str(i))) = nan(size(covars,1),1); end
PCs             = readtable('PCA1to10.csv');
for i = 1:size(IDs,1)
    if sum(strcmp(PCs.ID,IDs.ID{i})) == 1
        for ii = 1:10, covars.(strcat('PC',num2str(ii)))(i) = ...
                PCs{strcmp(PCs.ID,IDs.ID{i}),ii+1}; end
    end
end
% Response
resps       = table();
resps.INS   = Y.INS; 
pVals = nan(size(Xs,2),1);
effects = nan(size(Xs,2),1);
for i = 1:size(Xs,2)
    linmdl      = fitglm([Xs(:,i),covars,resps], 'CategoricalVars',...
        categ,'Distribution','binomial');
    pVals(i)    = linmdl.Coefficients.pValue(2);
    effects(i)  = linmdl.Coefficients.Estimate(2);
end
[pValsSort,idxs] = sort(pVals);
effectsSort      = effects(idxs);
predictorsSort   = Xs.Properties.VariableNames(idxs)';
nearestgene      = cell(size(predictorsSort));
for i = 1:length(predictorsSort)
    nearestgene{i} = snp2gene.NearestGene{ismember(snp2gene.rsID,predictorsSort{i})};
end
insTables = {};
for i = 1:size(predictorsSort,2)
    insTables{i} = table();
    insTables{i}.rsID = predictorsSort(:,i);
    insTables{i}.NearestGene = nearestgene(:,i);
    insTables{i}.Effect = num2str(round(effectsSort(:,i),4));
    insTables{i}.pValue = num2str(round(pValsSort(:,i),4));
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%%%%%%%% UNIVARIATE LINEAR REGRESSION model predict PSG WITH SNPS + covariates
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% Predictors
Xs      = X;
% Responses
Ys      = table();
Ys.ARI  = Y.ARI;
Ys.AW25 = Y.AW25;
Ys.AWN125  = Y.AWN125;
Ys.PLMI = Y.PLMI;
Ys.SE   = Y.SE;
Ys.SOL  = Y.SOL;
Ys.TST  = Y.TST;
% Covariates
covars          = table();
covars.AGE      = Y.AGE; 
covars.SEX      = Y.SEX; 
covars.Cohorts  = Y.Cohorts; 
covars.DEP      = Y.DEP; 
covars.ANX      = Y.ANX; 
for i = 1:10, covars.(strcat('PC',num2str(i))) = nan(size(covars,1),1); end
PCs             = readtable('PCA1to10.csv');
for i = 1:size(IDs,1)
    if sum(strcmp(PCs.ID,IDs.ID{i})) == 1
        for ii = 1:10, covars.(strcat('PC',num2str(ii)))(i) = ...
                PCs{strcmp(PCs.ID,IDs.ID{i}),ii+1}; end
    end
end
pVals = nan(size(Xs,2),size(Ys,2));
effects = nan(size(Xs,2),size(Ys,2));
alpha_bonf = 0.05/(size(Xs,2)*size(Ys,2));
predictors = repmat(Xs.Properties.VariableNames',1,size(Ys,2));
% Predict each Ys with each snp
for i = 1:size(Xs,2)
    for ii = 1:size(Ys,2) 
        linmdl          = fitglm([Xs(:,i),covars,Ys(:,ii)],...
            'CategoricalVars',categ);
        pVals(i,ii)     = linmdl.Coefficients.pValue(2);
        effects(i,ii)   = linmdl.Coefficients.Estimate(2);
    end
end
[pValsSort,idxs]= sort(pVals);
effectsSort     = effects(idxs);
predictorsSort  = predictors(idxs);
nearestgene      = cell(size(predictorsSort));
for i = 1:size(predictorsSort,1)
    for ii = 1:size(predictorsSort,2)
        nearestgene{i,ii} = snp2gene.NearestGene{ismember(snp2gene.rsID,predictorsSort{i,ii})};
    end
end
for i = 1:size(predictorsSort,2)
    for ii = 1:size(predictorsSort,1)
        uniTables{i}.ANXDEPeff{ismember(uniTables{i}.rsID,predictorsSort{ii,i})} = ...
            num2str(round(effectsSort(ii,i),4));
        uniTables{i}.ANXDEPpval{ismember(uniTables{i}.rsID,predictorsSort{ii,i})} = ...
            num2str(round(pValsSort(ii,i),4));
    end
end
%% MODEL 3 age, anxiety, depression, sex, cohort, pca, plmi
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%%%%%%%% UNIVARIATE LINEAR REGRESSION model predict PSG WITH SNPS + covariates
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% startup
% Load SNP to gene map
snp2gene         = readtable('SNP2GeneNew.csv');
load V2data.mat
% Predictors
Xs      = X;
% Responses
Ys      = table();
Ys.ARI  = Y.ARI;
Ys.AW25 = Y.AW25;
Ys.AWN125  = Y.AWN125;
Ys.SE   = Y.SE;
Ys.SOL  = Y.SOL;
Ys.TST  = Y.TST;
% Covariates
covars          = table();
covars.AGE      = Y.AGE; 
covars.SEX      = Y.SEX; 
covars.Cohorts  = Y.Cohorts; 
covars.DEP      = Y.DEP; 
covars.ANX      = Y.ANX; 
categ           = {'SEX','Cohorts','DEP','ANX'};
covars.PLMI     = Y.PLMI; 
for i = 1:10, covars.(strcat('PC',num2str(i))) = nan(size(covars,1),1); end
PCs             = readtable('PCA1to10.csv');
for i = 1:size(IDs,1)
    if sum(strcmp(PCs.ID,IDs.ID{i})) == 1
        for ii = 1:10, covars.(strcat('PC',num2str(ii)))(i) = ...
                PCs{strcmp(PCs.ID,IDs.ID{i}),ii+1}; end
    end
end
pVals = nan(size(Xs,2),size(Ys,2));
effects = nan(size(Xs,2),size(Ys,2));
alpha_bonf = 0.05/(size(Xs,2)*size(Ys,2));
predictors = repmat(Xs.Properties.VariableNames',1,size(Ys,2));
% Predict each Ys with each snp
for i = 1:size(Xs,2)
    for ii = 1:size(Ys,2) 
        linmdl          = fitglm([Xs(:,i),covars,Ys(:,ii)],...
            'CategoricalVars',categ);
        pVals(i,ii)     = linmdl.Coefficients.pValue(2);
        effects(i,ii)   = linmdl.Coefficients.Estimate(2);
    end
end
[pValsSort,idxs]= sort(pVals);
effectsSort     = effects(idxs);
predictorsSort  = predictors(idxs);
nearestgene      = cell(size(predictorsSort));
for i = 1:size(predictorsSort,1)
    for ii = 1:size(predictorsSort,2)
        nearestgene{i,ii} = snp2gene.NearestGene{ismember(snp2gene.rsID,predictorsSort{i,ii})};
    end
end
for i = 1:size(predictorsSort,2)
    for ii = 1:size(predictorsSort,1)
        uniTables{i}.PLMI{ismember(uniTables{i}.rsID,predictorsSort{ii,i})} = ...
            num2str(round(effectsSort(ii,i),4));
        uniTables{i}.ANXDEPpval{ismember(uniTables{i}.rsID,predictorsSort{ii,i})} = ...
            num2str(round(pValsSort(ii,i),4));
    end
end
%% Multivariate linear model predicting PSG
startup
snp2gene    = readtable('SNP2GeneNew.csv');
categ       = {'SEX','Cohorts','DEP','ANX'};
load V1newdata.mat

[~,score,~] = pca(table2array(X));
PCs         = array2table(score);
% Predictors
Xs      = X;
% Covariates
Covariates = table();
Covariates.AGE  = Y.AGE;
Covariates.SEX  = Y.SEX;
Covariates.BMI  = Y.BMI;
Covariates.Cohorts  = Y.Cohorts;
for i = 1:10, Covariates.(strcat('PC',num2str(i))) = nan(size(Covariates,1),1); end
PCs             = readtable('PCA1to10.csv');
for i = 1:size(IDs,1)
    if sum(strcmp(PCs.ID,IDs.ID{i})) == 1
        for ii = 1:10, Covariates.(strcat('PC',num2str(ii)))(i) = ...
                PCs{strcmp(PCs.ID,IDs.ID{i}),ii+1}; end
    end
end
% Responses
Ys      = table();
Ys.ARI  = Y.ARI;
Ys.AW25 = Y.AW25;
Ys.AWN125  = Y.AWN125;
Ys.PLMI = Y.PLMI;
Ys.SE   = Y.SE;
Ys.SOL  = Y.SOL;
Ys.TST  = Y.TST;

categ = {'SEX','Cohorts'};
% For each response
for i = 1:size(Ys,2)
        if strcmp(Ys.Properties.VariableNames(i), {'PLMI','INS'}) == 1 
            mdls  = fitglm([Xs,Covariates,Ys(:,i)],'CategoricalVars',categ,'Distribution','binomial');
        else 
            mdls  = fitglm([Xs,Covariates,Ys(:,i)],'CategoricalVars',categ); 
        end
        [pVal,idx]      = sort(mdls.Coefficients.pValue(2:end));
        effect          = mdls.Coefficients.Estimate(2:end);
        predictor       = mdls.CoefficientNames(2:end);
        pVals(:,i)      = pVal;
        effects(:,i)    = effect(idx);
        predictors(:,i) = predictor(idx);
end
pFDR            = mafdr(pVals(:),'BHFDR', 0.1);
pFDR            = reshape(pFDR,[],size(pVals,2));
pFDR            = array2table(pFDR,'VariableNames',Ys.Properties.VariableNames);

newSNPs = table2array(pFDR) < 0.05;

newPredictors = predictors(newSNPs);
newPredictors = newPredictors(contains(newPredictors,'rs'));
Xs = Xs(:,ismember(Xs.Properties.VariableNames,newPredictors));
clear pVals effects predictors
Ys.PLMI = Ys.PLMI >=15; 
for i = 1:size(Ys,2)
        if strcmp(Ys.Properties.VariableNames(i), {'PLMI','INS'}) == 1 
            mdls  = fitglm([Xs,Covariates,Ys(:,i)],'CategoricalVars',categ,'Distribution','binomial');
        else 
            mdls  = fitglm([Xs,Covariates,Ys(:,i)],'CategoricalVars',categ); 
        end
        [pVal,idx]      = sort(mdls.Coefficients.pValue(2:end));
        effect          = mdls.Coefficients.Estimate(2:end);
        predictor       = mdls.CoefficientNames(2:end);
        pVals(:,i)      = pVal;
        effects(:,i)    = effect(idx);
        predictors(:,i) = predictor(idx);
end
pFDR            = mafdr(pVals(:),'BHFDR', 0.1);
pFDR            = reshape(pFDR,[],size(pVals,2));
pFDR            = array2table(pFDR,'VariableNames',Ys.Properties.VariableNames);
%% Total risk and PSG variables
load RisksV1.mat
load V1data.mat
categ = {'SEX','Cohorts'};
% Predictors
Xs      = array2table(risks,'VariableNames',{'Risk'});
% Responses
Ys      = table();
Ys.ARI  = Y.ARI;
Ys.AW25 = Y.AW25;
Ys.AWN125  = Y.AWN125;
Ys.PLMI = Y.PLMI;
Ys.SE   = Y.SE;
Ys.SOL  = Y.SOL;
Ys.TST  = Y.TST;
% Covariates
covars          = table();
covars.AGE      = Y.AGE; 
covars.SEX      = Y.SEX; 
covars.Cohorts  = Y.Cohorts; 
for i = 1:10, covars.(strcat('PC',num2str(i))) = nan(size(covars,1),1); end
PCs             = readtable('PCA1to10.csv');
for i = 1:size(IDs,1)
    if sum(strcmp(PCs.ID,IDs.ID{i})) == 1
        for ii = 1:10, covars.(strcat('PC',num2str(ii)))(i) = ...
                PCs{strcmp(PCs.ID,IDs.ID{i}),ii+1}; end
    end
end
pVals = nan(size(Xs,2),size(Ys,2));
effects = nan(size(Xs,2),size(Ys,2));
alpha_bonf = 0.05/(size(Xs,2)*size(Ys,2));
predictors = repmat(Xs.Properties.VariableNames',1,size(Ys,2));
% Predict each Ys with each snp
for i = 1:size(Xs,2)
    for ii = 1:size(Ys,2) 
        linmdl          = fitglm([Xs(:,i),covars,Ys(:,ii)],...
            'CategoricalVars',categ);
        pVals(i,ii)     = linmdl.Coefficients.pValue(2);
        effects(i,ii)   = linmdl.Coefficients.Estimate(2);
    end
end
[pValsSort,idxs]= sort(pVals);
effectsSort     = effects(idxs);
predictorsSort  = predictors(idxs);

plot(risks,Ys.TST,'*')
load V2data.mat
[coeffY,scoreY,latentY] = pca(Y{:,2:end-1});
[coeffX,scoreX,latentX] = pca(X{:,:});
figure,
    g = gscatter(scoreY(:,1),scoreY(:,2),Y.Cohorts,'brg','xos');  
    xlabel('PC 1'), grid on
    ylabel('PC 2'), box off
    hLeg = legend({'WSC','MrOS','SSC'}); 
    
figure,
    g = gscatter(scoreY(:,1),scoreY(:,2),Y.INS,'br','xo');  
    xlabel('PC 1'), grid on
    ylabel('PC 2'), box off
    hLeg = legend({'Non-Insomnia','Insomnia'});
    
figure,
    g = gscatter(scoreX(:,1),scoreX(:,2),Y.Cohorts,'brg','xos');  
    xlabel('PC 1'), grid on
    ylabel('PC 2'), box off
    hLeg = legend({'WSC','MrOS','SSC'}); 
Tf = isoutlier(scoreX(:,1));
scoreX(Tf,:) = [];
cohort  = Y.Cohorts;
cohort(Tf,:) = [];
ins     = Y.INS;
ins(Tf,:) = [];
scoreY(Tf,:) = [];
figure,
    g = gscatter(scoreX(:,1),scoreX(:,2),cohort,'brg','xos');  
    xlabel('PC 1'), grid on
    ylabel('PC 2'), box off
    hLeg = legend({'WSC','MrOS','SSC'}); 
figure,
    g = gscatter(scoreY(:,1),scoreY(:,2),cohort,'brg','xos');  
    xlabel('PC 1'), grid on
    ylabel('PC 2'), box off
    hLeg = legend({'WSC','MrOS','SSC'}); 
figure,
    g = gscatter(scoreX(:,1),scoreX(:,2),ins,'br','xo');  
    xlabel('PC 1'), grid on
    ylabel('PC 2'), box off
    hLeg = legend({'Non-Insomnia','Insomnia'});
figure,
    g = gscatter(scoreY(:,1),scoreY(:,2),ins,'br','xo');  
    xlabel('PC 1'), grid on
    ylabel('PC 2'), box off
    hLeg = legend({'Non-Insomnia','Insomnia'});

