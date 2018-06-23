startup
pth         = 'Notebooks\Results\SHHS\Both\';
lsting      = extractfield(dir(pth),'name')';
spectral    = lsting(contains(lsting,'Ys'));
temporal    = lsting(contains(lsting,'Yt'));
%% 
temporalVal = temporal(contains(temporal,'val'));    
spectralVal = spectral(contains(spectral,'val'));    

preds       = nan(500,2);
trues       = nan(500,1);
idx1        = 0;
idx2        = 0;
for i = 1:size(spectralVal,1)
        curSpec = csvread(strcat(pth,spectralVal{i}));
        curTemp = csvread(strcat(pth,temporalVal{i}));
    if  sum(contains(temporalVal{i},'Pred')) == 1
        preds(1+idx1:idx1+size(curSpec,1),1) = curSpec;
        preds(1+idx1:idx1+size(curSpec,1),2) = curTemp;
        idx1 = idx1 + size(curSpec,1);
    elseif isequal(curSpec,curTemp)
        trues(1+idx2:idx2+size(curSpec,1),1) = curSpec;
        idx2 = idx2 + size(curSpec,1);
    end
end
preds(~any(~isnan(preds), 2),:)=[];
trues(~any(~isnan(trues), 2),:)=[];
all                             = [preds,trues];
% classification learner
%% Test
temporalTest    = temporal(contains(temporal,'test'));    
spectralTest    = spectral(contains(spectral,'test'));    
preds           = nan(500,2);
trues           = nan(500,1);
idx1            = 0;
idx2            = 0;
for i = 1:size(temporalTest,1)
        curSpec = csvread(strcat(pth,temporalTest{i}));
        curTemp = csvread(strcat(pth,spectralTest{i}));
    if  sum(contains(spectralTest{i},'Pred')) == 1
        preds(1+idx1:idx1+size(curSpec,1),1) = curSpec;
        preds(1+idx1:idx1+size(curSpec,1),2) = curTemp;
        idx1 = idx1 + size(curSpec,1);
    elseif isequal(curSpec,curTemp)
        trues(1+idx2:idx2+size(curSpec,1),1) = curSpec;
        idx2 = idx2 + size(curSpec,1);
    end
end
preds(~any(~isnan(preds), 2),:)=[];
trues(~any(~isnan(trues), 2),:)=[];
all = [all;[preds,trues]];
% all = [preds,trues];
% classification learner application launc using 'all'
% view tree
% figure,view(trainedModel.ClassificationTree,'Mode','graph');

x = [0.673185185	0.62962963	0.576271186	0.528301887	0.641860465	0.666666667	0.599508197	0.518269231	0.537037037	0.624];
mean(x(1:2:end))

%% Looking into raw features
cohort      = 'WSC';
pth         = strcat('Features\',cohort,'\');
lsting      = extractfield(dir(pth),'name')';
lst         = lsting(contains(lsting,'h5'));
data        = [];
spec        = [];
temp        = [];
for i = 1:size(lst,1)
    data = [data;h5read(strcat(pth,lst{i}),'/demographics')];
    spec = [spec;mean(h5read(strcat(pth,lst{i}),'/spectral'))];
    temp = [temp;mean(h5read(strcat(pth,lst{i}),'/temporal'))];
end

ins = data(data(:,end) == 1,:);
con = data(data(:,end) == 0,:);

all     = [spec,temp,data(:,end)];
spec    = [spec,data(:,end)];
temp    = [temp,data(:,end)];

% Logistic regression
mdl = fitglm(temp(:,1:end-1),data(:,end),'Distribution','binomial');
% age, agediff, ahi, anx, ari, bmi, dep, ess, sex, ins












