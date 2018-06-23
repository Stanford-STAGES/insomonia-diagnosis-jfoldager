startup
shouldplot = false;

data = readtable('AllData.csv');
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
%%
if shouldplot
    for i = 1:length(yCols)
        figure,
        gscatter(b(:,1),b(:,2),yGroup(:,i))
        legend({'[0,0.05] ','[0.05,0.25]',...
            '[0.25,0.50]',...
            '[0.50,0.75]',...
            '[0.75,0.95]',...
            '[0.95,1]'})
        title(yLabels{i})
    end
end
%% K-means
[idx,C] = kmeans(dataZ,5);

%% Gaussian Mixture Model
% Cross validation to choose K
K = 10;
C = 20;
nLogLsTrain = nan(K,C);
nLogLsTest = nan(K,C);
gmData = Xz;
for k = 1:K
    fprintf('Fold nr. %i\n',k);
    [Train, Test] = crossvalind('HoldOut', size(dataZ,1), 0.25);    
    for c = 1:C
        try
            GMModel = fitgmdist(gmData(Train,:),c,'Replicates',100, 'RegularizationValue',0.1,'CovarianceType','diagonal');
            nLogLsTrain(k,c) = GMModel.NegativeLogLikelihood;

            [~,nlogL] = posterior(GMModel,gmData(Test,:));
            nLogLsTest(k,c) = nlogL;
        catch
            
        end
    end
end


yyaxis left
plot(nanmean(nLogLsTrain),'LineWidth',3)
xlabel('Mixture components')
ylabel('Neg. Log-Likelihood')
yyaxis right
plot(nanmean(nLogLsTest),'LineWidth',3)
ylabel('Neg. Log-Likelihood')
legend('Train','Test'), grid on

save('gmmclusterselection','nLogLsTest','nLogLsTrain')

% CGobj = clustergram(dataZ,'Cluster','row', 'RowPDist', 'jaccard');
% get(cgo)
% set(gco,'Linkage','complete','Dendrogram',3)
%% Association
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

xProbabilities = sum(xApriori)./size(xApriori,1);
xHomoHetro = reshape(xProbabilities,3,[]);
MinAF = xHomoHetro(2,:)/2 + xHomoHetro(3,:);
MajAF = xHomoHetro(2,:)/2 + xHomoHetro(1,:);
nMinThresh = sum(MinAF < 0.05); % should be zero
nMajThresh = sum(MajAF > 0.95); % should be zero
 
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



allApriori = [xApriori,yApriori] > 0.5;
% Apriori algorithm
csvwrite('apri.csv',allApriori);
[AssocRules, FrqItemsets, Summary1, Summary2] =...
    apriori('apri.csv',0.1, 0.5);



