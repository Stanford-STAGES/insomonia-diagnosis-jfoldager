startup
shouldplot = false;

data = readtable('AllData.csv');
% load negativeloglikelihood.mat
if shouldplot, plot(nanmean(nLogLs)), end
nGm = 6;
labels = data.Properties.VariableNames(:,3:end);
cohorts = strcmp(data.Cohort,'WSC');
cohorts = cohorts + strcmp(data.Cohort,'MrOS')*2;
cohorts = cohorts + strcmp(data.Cohort,'SSC')*3;

dataArray = table2array(data(:,3:end));
dataZ = zscore(dataArray);
yCols = 21;
xCols = 1:(size(dataArray,2)-yCols); 
yCols = (xCols(end)+1):(xCols(end)+yCols);
xLabels = labels(xCols);
yLabels = strrep(labels(yCols),'h_','');
X = dataArray(:,xCols);
Xz= dataZ(:,xCols);
Y = dataArray(:,yCols);
Yz= dataZ(:,yCols);

[pcaAllC,pcaAllS,pcaAllZ] = pca(dataZ);
[pcaXzC,pcaXzS,pcaXzZ] = pca(Xz);
[pcaYzC,pcaYzS,pcaYzZ] = pca(Yz);
if shouldplot
    figure
    gscatter(pcaXzS(:,1),pcaXzS(:,2),cohorts)
    figure
    gscatter(pcaYzS(:,1),pcaYzS(:,2),cohorts)
end

yGroup = nan(size(Y));
idx = 1;
for i = 1:length(yCols)
    q = quantile(dataArray(:,yCols(i)),[0.05,0.25,0.50,0.75,0.95]);
    yGroup(:,i) = (dataArray(:,yCols(i)) < q(1))*1;
    yGroup(:,i) = yGroup(:,i) + ((dataArray(:,yCols(i)) >= q(1)).*(dataArray(:,yCols(i)) < q(2)))*2;
    yGroup(:,i) = yGroup(:,i) + ((dataArray(:,yCols(i)) >= q(2)).*(dataArray(:,yCols(i)) < q(3)))*3;
    yGroup(:,i) = yGroup(:,i) + ((dataArray(:,yCols(i)) >= q(3)).*(dataArray(:,yCols(i)) < q(4)))*4;
    yGroup(:,i) = yGroup(:,i) + ((dataArray(:,yCols(i)) >= q(4)).*(dataArray(:,yCols(i)) < q(5)))*5;
    yGroup(:,i) = yGroup(:,i) + ((dataArray(:,yCols(i)) >= q(5)))*6;
end

xApriori = nan(size(dataArray,1),length(xCols)*3);
newXlabels = strcat(repmat(xLabels,3,1),repmat({'|0';'|1';'|2'},1,size(xLabels,2)));
newXlabels = newXlabels(:);
idx = 1;
for i = 1:length(xCols)
    xApriori(:,idx) = dataArray(:,xCols(i)) < 0.5;
    xApriori(:,idx+1) = (dataArray(:,xCols(i)) >= 0.5).*(dataArray(:,xCols(i)) < 1.5);
    xApriori(:,idx+2) = dataArray(:,xCols(i)) >= 1.5;
    idx = idx + 3;
end

classifierData = [array2table(X,'VariableNames',xLabels) ...
    ,array2table(yGroup,'VariableNames',yLabels)];


yApriori = nan(size(dataArray,1),length(yCols)*6);
newYlabels = strcat(repmat(yLabels,6,1),repmat({'q1';'q2';'q3';'q4';'q5';'q6'},1,size(yLabels,2)));
newYlabels = newYlabels(:)';
idx = 1;
qs = nan(length(yCols),5);
for i = 1:length(yCols)
    q = quantile(dataArray(:,yCols(i)),[0.05,0.25,0.50,0.75,0.95]);
    qs(i,:) = q;
    yApriori(:,idx) = dataArray(:,yCols(i)) < q(1);
    yApriori(:,idx+1) = (dataArray(:,yCols(i)) >= q(1)).*(dataArray(:,yCols(i)) < q(2));
    yApriori(:,idx+2) = (dataArray(:,yCols(i)) >= q(2)).*(dataArray(:,yCols(i)) < q(3));
    yApriori(:,idx+3) = (dataArray(:,yCols(i)) >= q(3)).*(dataArray(:,yCols(i)) < q(4));
    yApriori(:,idx+4) = (dataArray(:,yCols(i)) >= q(4)).*(dataArray(:,yCols(i)) < q(5));
    yApriori(:,idx+5) = (dataArray(:,yCols(i)) >= q(5));
    idx = idx + 6;
end
% Only two groups
% yApriori = Y>quantile(Y,0.5);

xProbabilities = sum(xApriori)./size(xApriori,1);
xHomoHetro = reshape(xProbabilities,3,[]);
MinAF = xHomoHetro(2,:)/2 + xHomoHetro(1,:);
MajAF = xHomoHetro(2,:)/2 + xHomoHetro(3,:);
nMinThresh = sum(MinAF < 0.05); % should be zero
if nMinThresh > 0
    minFreqSNPs = xLabels(MinAF < 0.05);
    xLabels = xLabels(~(MinAF < 0.05));
    X = X(:,~(MinAF < 0.05));
    Xz = Xz(:,~(MinAF < 0.05));
    xCols = xCols(~(MinAF < 0.05));
    data(:,ismember(data.Properties.VariableNames,minFreqSNPs)) = [];
    xApriori = nan(size(dataArray,1),(length(xCols))*3);
    newXlabels = strcat(repmat(xLabels,3,1),repmat({'|0';'|1';'|2'},1,size(xLabels,2)));
    newXlabels = newXlabels(:);
    idx = 1;
    for i = 1:length(xCols)
        xApriori(:,idx) = dataArray(:,xCols(i)) < 0.5;
        xApriori(:,idx+1) = (dataArray(:,xCols(i)) >= 0.5).*(dataArray(:,xCols(i)) < 1.5);
        xApriori(:,idx+2) = dataArray(:,xCols(i)) >= 1.5;
        idx = idx + 3;
    end
end

P = nan(size(xApriori,2),size(yApriori,2));
allels = [];
for allel = 1:size(xApriori,2)
    for psg = 1:size(yApriori,2)
        if sum(xApriori(:,allel)) > 0
            P(allel,psg) = sum((yApriori(:,psg)).*(xApriori(:,allel)))/sum(xApriori(:,allel));
        else
            allels = [allels, allel];
        end
    end
end
allels = unique(allels);
% P(any(isnan(P')),:) = [];
[P,idx] = sort(P,'descend');

extremesHigh    = P(:,6:6:size(P,2));
idxHigh         = idx(:,6:6:size(idx,2));
extremesLow     = P(:,1:6:size(P,2));
idxLow          = idx(:,1:6:size(idx,2));

%%
nMdls   = 9; %size(Y,2)
mdls    = cell(1,nMdls);
lassomdls    = cell(1,nMdls);
ridgemdls    = cell(1,nMdls);
lassoInfo    = cell(1,nMdls);
ridgeInfo    = cell(1,nMdls);
pVals   = ones(size(X,2),nMdls);
effects = ones(size(X,2),nMdls);
SNPs    = cell(size(X,2),nMdls);
for i = 1:nMdls
%     mdls{i}         = fitglm([X,Y(:,i)]);
%     [lassomdls{i},lassoInfo{i}] =  lassoglm(X,Y(:,i),'normal','CV',10);
%     [ridgemdls{i},ridgeInfo{i}] =  lassoglm(X,Y(:,i),'normal','CV',10,'Alpha',0.01);
%     [pVal,idx]      = sort(mdls{i}.Coefficients.pValue(2:end));
%     effect          = mdls{i}.Coefficients.Estimate(2:end);
%     pVals(:,i)      = pVal;
%     effects(:,i)    = effect(idx);
%     SNPs(:,i)       = allSNPs(idx);
disp(i)
end

% save('lassoresults','lassomdls','FitInfo')
% save('ridgeresults','ridgemdls','ridgeInfo')
load lassoresults 
load ridgeresults

yPhe = 1;
figure,lassoPlot(lassomdls{yPhe},FitInfo{yPhe},'plottype','CV');
legend('show') % Show legend
figure,boxplot(ridgemdls{yPhe}')
xticklabels({''}), ylabel('TST'), xlabel('SNPs')
for yPhe = 1:nMdls
    [m,idx] = min(FitInfo{yPhe}.Deviance);
    betas(:,yPhe) = lassomdls{yPhe}(:,idx);
end
FitInfo{yPhe}.Lambda(idx)
intercept = FitInfo{yPhe}.Intercept(idx);
lassos = (abs(betas) > 0);
lassoSNPs = xLabels(lassos);
lassoX = X(:,lassos);
lassoXz = Xz(:,lassos);

Yy = intercept + sum(X.*repmat(betas',size(X,1),1),2);
plot(Y(:,1),Yy,'.')

%% 
% gmmData = Yz;
gmmData = Yz;
GMM = fitgmdist(gmmData,6,'CovarianceType','diagonal');
[ps,nLogs] = posterior(GMM,gmmData);
[m, i] = max(ps,[],2);
classifierData2 = [array2table(X,'VariableNames',xLabels) ...
    ,array2table(yApriori(:,1),'VariableNames',{yLabels{1}})];

if shouldplot, figure, histogram(m); figure, histogram(i); end
if shouldplot, figure, scatter3(pcaXzS(:,1),pcaXzS(:,2),pcaXzS(:,3),[],i); end
% SNPs, Subjects in group n, Dosages
figure,
subplot(2,3,1),imagesc(X(i == 1,:)), colorbar, colormap(jet(5)), title('1')
subplot(2,3,2),imagesc(X(i == 2,:)), colorbar, colormap(jet(5)), title('2')
subplot(2,3,3),imagesc(X(i == 3,:)), colorbar, colormap(jet(5)), title('3')
subplot(2,3,4),imagesc(X(i == 4,:)), colorbar, colormap(jet(5)), title('4')
subplot(2,3,5),imagesc(X(i == 5,:)), colorbar, colormap(jet(5)), title('5')
subplot(2,3,6),imagesc(X(i == 6,:)), colorbar, colormap(jet(5)), title('6')
% Each clusters mean dosage
figure,
Xmeans = [mean(X(i == 1,:)); ...
            mean(X(i == 2,:)); ...
            mean(X(i == 3,:)); ...
            mean(X(i == 4,:)); ...
            mean(X(i == 5,:)); ...
            mean(X(i == 6,:))];
imagesc(Xmeans), colorbar, colormap(jet(10))
% Range/variation in X data (SNPs)
figure,
subplot(2,1,1),imagesc(range(Xmeans)), colorbar, colormap(jet(10))
subplot(2,1,2),imagesc(var(Xmeans)), colorbar, colormap(jet(10))
% Each clusters mean dosage
figure, yPhe = 2;
Xmeans = [mean(X(yGroup(:,yPhe) == 1,:)); ...
            mean(X(yGroup(:,yPhe) == 2,:)); ...
            mean(X(yGroup(:,yPhe) == 3,:)); ...
            mean(X(yGroup(:,yPhe) == 4,:)); ...
            mean(X(yGroup(:,yPhe) == 5,:)); ...
            mean(X(yGroup(:,yPhe) == 6,:))];
imagesc(Xmeans), colorbar, colormap(jet(10))




figure,imagesc(xApriori(i == 1,:)), colorbar, colormap(gray(2))
figure,imagesc(xApriori(i == 2,:)), colorbar, colormap(gray(2))
figure,imagesc(xApriori(i == 3,:)), colorbar, colormap(gray(2))
figure,imagesc(xApriori(i == 4,:)), colorbar, colormap(gray(2))
figure,imagesc(xApriori(i == 5,:)), colorbar, colormap(gray(2))
figure,imagesc(xApriori(i == 6,:)), colorbar, colormap(gray(2))

figure,plot(mean(xApriori(i == 1,:))), colorbar, colormap(gray(2))
figure,plot(mean(xApriori(i == 2,:))), colorbar, colormap(gray(2))
figure,plot(mean(xApriori(i == 3,:))), colorbar, colormap(gray(2))
figure,plot(mean(xApriori(i == 4,:))), colorbar, colormap(gray(2))
figure,plot(mean(xApriori(i == 5,:))), colorbar, colormap(gray(2))
figure,plot(mean(xApriori(i == 6,:))), colorbar, colormap(gray(2))

figure,plot(mean(xApriori(i == 1,:)) - mean(xApriori(i == 2,:))), colorbar, colormap(gray(2))
figure,plot(mean(xApriori(i == 1,:)) - mean(xApriori(i == 3,:))), colorbar, colormap(gray(2))
figure,plot(mean(xApriori(i == 1,:)) - mean(xApriori(i == 4,:))), colorbar, colormap(gray(2))
figure,plot(mean(xApriori(i == 1,:)) - mean(xApriori(i == 5,:))), colorbar, colormap(gray(2))
figure,plot(mean(xApriori(i == 1,:)) - mean(xApriori(i == 6,:))), colorbar, colormap(gray(2))







s