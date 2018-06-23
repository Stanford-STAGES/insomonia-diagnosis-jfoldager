%% Clustering PSG
startup
load V1data.mat
% Responses
Ys      = table();
Ys.ARI  = Y.ARI;
Ys.AW25 = Y.AW25;
Ys.AWN125  = Y.AWN125;
Ys.PLMI = Y.PLMI;
Ys.SE   = Y.SE;
Ys.SOL  = Y.SOL;
Ys.TST  = Y.TST;
Yz      = zscore(table2array(Ys));
%% PCA
[coeffY,scoreY,latentY] = pca(zscore(table2array(Ys)));
expVar = cumsum(latentY/sum(latentY));
% Explained variance
figure, stem(expVar)
% Biplot
figure,biplot(coeffY(:,1:2),'VarLabels',Ys.Properties.VariableNames)
%% GMM
% nCom = find(expVar >= 0.8); com = 7; %1:nCom(1);
% for k = 1:10
%    [Train, Test] = crossvalind('HoldOut', size(scoreY,1), 0.25);
%    for c = 1:20
%        mdl = fitgmdist(scoreY(Train,com),c,'CovarianceType','diagonal',...
%            'Replicates',10);
%        nLogLiTrain(k,c) = mdl.NegativeLogLikelihood;
%        [~,nLogLiTest(k,c)] = posterior(mdl,scoreY(Test,com));
%        BIC(k,c) = mdl.BIC;
%        AIC(k,c) = mdl.AIC;
%    end    
% end
% 
% % Plot log likelihood
% figure,grid on
% yyaxis left
% plot(mean(nLogLiTrain))
% ylabel('Neg. Log-Likelihood'), 
% yyaxis right
% plot(mean(nLogLiTest)), legend({'Train','Test'})
% xlabel('Components'), %ylim([min(mean(nLogLiTest)) - 1000,max(mean(nLogLiTest))])
% box off, 
% set(gca,'color','none')
% % print(strcat(plotPath3,'GMM_PSG_CV_LogLi'),'-depsc');
% % Plot BIC, AIC
% figure, grid on
% yyaxis left
% plot(mean(AIC)), 
% yyaxis right
% plot(mean(BIC)), legend({'AIC','BIC'})
% xlabel('Components')
% box off, 
% set(gca,'color','none')
% print(strcat(plotPath3,'GMM_PSG_AIC_BIC'),'-depsc');
% Pick K
K = 4;
mdl = fitgmdist(scoreY,K,'CovarianceType','diagonal');
[P,~] = posterior(mdl,scoreY);
[~,cLabel] = max(P,[],2);
cDist = histcounts(cLabel)./sum(histcounts(cLabel));
[~,outlierC] = min(cDist);
% Remove 'outliers'
scoreY(cLabel == outlierC,:) = []; 
Yz(cLabel == outlierC,:) = []; 
Ys(cLabel == outlierC,:) = []; 
Y(cLabel == outlierC,:) = []; 
X(cLabel == outlierC,:) = []; 
cLabel(cLabel == outlierC,:) = []; 
% Plot clusters
figure,gscatter(scoreY(:,1),scoreY(:,2),cLabel,'brgy','oxsx')
box off, grid on
xlabel('PC 1')
ylabel('PC 2')
set(gca,'color','none')
% legend({num2str(round(cDist(1),4)),num2str(round(cDist(2),4)),...
%     num2str(round(cDist(3),4)),num2str(round(cDist(4),4))});
% print(strcat(plotPath3,'PSG_PCA_Clusters'),'-depsc');
% Cluster distributions
figure, histogram(cLabel)
% By cohort
figure,gscatter(scoreY(:,1),scoreY(:,2),Y{:,end},'brg','xos')
hLeg = legend({'WSC','MrOS','SSC'}, 'Location', 'southwest');
box off, grid on
xlabel('PC 1')
ylabel('PC 2')
set(gca,'color','none')
% print(strcat(plotPath3,'PSG_PCA_Cohort_Clusters'),'-depsc');
figure,scatter(scoreY(:,1),scoreY(:,2),[],Y.AGE), colorbar

% Bar chart (each cohort)
data = nan(length(unique(cLabel)),length(unique(Y{:,end})));
for i = 1:length(unique(cLabel))
    for ii = 1:length(unique(Y{:,end}))
        data(i,ii) = sum((cLabel == i).*(Y{:,end} == ii));
    end
end
hb = bar(data./(sum(data)));
set(hb(1), 'FaceColor','b')
set(hb(2), 'FaceColor','r')
set(hb(3), 'FaceColor','g')


