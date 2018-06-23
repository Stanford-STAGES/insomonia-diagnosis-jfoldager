startup
% % % % % % % % % % % % % % % % % % % % % % % % % % % % 
% 1: Never (0)
% 2: Rarely (1x/month or less)
% 3: Sometimes (2-4x/month)
% 4: Often (5-15x/month)
% 5: Almost Always (16-30x/month)
% % % % % % % % % % % % % % % % % % % % % % % % % % % % 
% 1: All of the time
% 2: Most of the time
% 3: A good bit of the time
% 4: Some of the time
% 5: A little of time
% 6: None of the time
% % % % % % % % % % % % % % % % % % % % % % % % % % % % 
insThresholdV1 = 1;
insThresholdV2 = 5;
V1 = readtable('shhs1-dataset-0.13.0.csv');
fallingasleepV1 = (V1.TFA02<=insThresholdV1)*1;
wakingupnightV1 = (V1.WUDNRS02<=insThresholdV1)*1;
morningawakeV1  = (V1.WU2EM02<=insThresholdV1)*1;
INSv1           = fallingasleepV1+wakingupnightV1+morningawakeV1;

V2 = readtable('shhs2-dataset-0.13.0.csv');
fallingasleepV2 = (V2.sh308a>=insThresholdV2)*1;
wakingupnightV2 = (V2.sh308b>=insThresholdV2)*1;
morningawakeV2  = (V2.sh308c>=insThresholdV2)*1;
INSv2           = fallingasleepV2+wakingupnightV2+morningawakeV2;

warning ('off','all');
DATA            = table();
idx             = 1;
for i = 1:size(V1,1)
    if sum(V1.nsrrid(i) == V2.nsrrid) == 1
        DATA.ID(idx)  = V1.nsrrid(i);
        DATA.AGE(idx) = V1.age_s1(i);
        DATA.AGEDIFF(idx) = V2.age_s2(V1.nsrrid(i) == V2.nsrrid) - V1.age_s1(i);
        DATA.AHI(idx) = V1.ahi_a0h4(i);
        DATA.ANX(idx) = ((V1.NRVOUS25(i)<=3)*1) > 0;
        DATA.ARI(idx) = V1.ai_all(i);
        DATA.BMI(idx) = V1.bmi_s1(i);
        DATA.DEP(idx) = ((V1.DOWN25(i)<=3)*1+(V1.BLUE25(i)<=3)*1) > 0;
        DATA.ESS(idx) = V1.ESS_s1(i);
        DATA.INSv1(idx) = INSv1(i)>0;
        DATA.INSv2(idx) = INSv2(V1.nsrrid(i) == V2.nsrrid)>0;
        DATA.INSFAv1(idx) = fallingasleepV1(i)>0;
        DATA.INSFAv2(idx) = fallingasleepV2(V1.nsrrid(i) == V2.nsrrid)>0;
        DATA.INSWUv1(idx) = wakingupnightV1(i)>0;
        DATA.INSWUv2(idx) = wakingupnightV2(V1.nsrrid(i) == V2.nsrrid)>0;
        DATA.INSEMv1(idx) = morningawakeV1(i)>0;
        DATA.INSEMv2(idx) = morningawakeV2(V1.nsrrid(i) == V2.nsrrid)>0;
        DATA.SEX(idx) = V1.gender(i) == 1;
        idx = idx + 1;
    end   
end
insRepData = DATA((DATA.INSv2-DATA.INSv1)==1,:);
insRepData = rmmissing(insRepData);
% Subtypes
criterion   = (((DATA.INSFAv2-DATA.INSFAv1)==1)*1 + ...
    ((DATA.INSFAv2-DATA.INSFAv1)==1)*1 + ...
    ((DATA.INSFAv2-DATA.INSFAv1)==1)*1)>0;
insMulData      = DATA(criterion,:);
insfaRepData    = DATA((DATA.INSFAv2-DATA.INSFAv1)==1,:);
insfaRepData    = rmmissing(insfaRepData);
inswuRepData    = DATA((DATA.INSWUv2-DATA.INSWUv1)==1,:);
inswuRepData    = rmmissing(inswuRepData);
insemRepData    = DATA((DATA.INSEMv2-DATA.INSEMv1)==1,:);
insemRepData    = rmmissing(insemRepData);
a               = strcat('*',string(insRepData.ID),'*');
DATA            = rmmissing(DATA);
%% Find categorical distances
A       = insRepData;
A(:,10:17) = [];
C       = DATA(sum(DATA{:,10:17},2) == 0,:);
C(:,10:17) = [];
CATexp  = [A.DEP,A.ANX,A.SEX]*1; A.DEP = []; A.ANX = []; A.SEX = []; A.AGEDIFF = [];
CATcon  = [C.DEP,C.ANX,C.SEX]*1; C.DEP = []; C.ANX = []; C.SEX = []; C.AGEDIFF = [];
catDist = pdist2(CATexp,CATcon,'jaccard');
catDist(isnan(catDist)) = 0;
%% Find euclidean distances
covariatesExp   = zscore(table2array(A(:,2:end)));
covariatesCon   = zscore(table2array(C(:,2:end)));
dists           = pdist2(covariatesExp,covariatesCon)/size(covariatesCon,2) ...
    + catDist/size(catDist,2);
[mins,idxss]    = sort(dists');
usedIdx         = [];
controlGrpIDs   = nan(size(covariatesExp,1),1);
for i = 1:size(idxss,2)
    if sum(ismember(usedIdx,idxss(1,i))) == 0
        controlGrpIDs(i) = C.ID(idxss(1,i));
        usedIdx = [usedIdx;idxss(1,i)];
    else
        found = false;
        ii = 2;
        while ~found
            if sum(ismember(usedIdx,idxss(ii,i))) == 0
            controlGrpIDs(i) = C.ID(idxss(ii,i));
            usedIdx         = [usedIdx;idxss(ii,i)];
            found           = true;
            end
            ii = ii + 1;
        end
    end
end
experimentGrpIDs = A.ID;
% sum(ismember(experimentGrpIdx,controlGrpIdx)) % should be zero!
a = strcat('*',string(experimentGrpIDs),'*');
b = strcat('*',string(controlGrpIDs),'*');
writetable([insRepData(:,1:9),insRepData(:,end)],'SHHSfutureInsomniacs.csv');
conRepData = DATA(ismember(DATA.ID,controlGrpIDs),:);
writetable([conRepData(:,1:9),conRepData(:,end)],'SHHSfutureControls.csv');
ALL = [[conRepData(:,1:9),conRepData(:,end)];[insRepData(:,1:9),insRepData(:,end)]];
ALL.INS = [zeros(size(conRepData,1),1);ones(size(insRepData,1),1)];
ALL.ID = string(ALL.ID);
writetable(ALL,'SHHSfutureAll.csv');
%% Other plots
T       = insRepData;
% Anxiety, depression and insomnia a.f.o age
data    = [T.AGE, T.ANX, T.DEP, T.INSv1];
sex     = [T.SEX];
% Males
data(sex == 0,:) = [];
[~,idxs] = sort(data(:,1));
data = data(idxs,:);
data(any(isnan(data), 2), :) = [];
data = [data(:,1), movmean(data(:,2:end),100)];
figure,
hold on
plot(data(:,1),data(:,2),'-r')
plot(data(:,1),data(:,3),'-g')
plot(data(:,1),data(:,4),'-b')
grid on, box off
xlabel('Age'),ylabel('Probability')
% Females
data    = [T.AGE, T.ANX, T.DEP, T.INSv1];
data(sex == 1,:) = [];
[~,idxs] = sort(data(:,1));
data = data(idxs,:);
data(any(isnan(data), 2), :) = [];
data = [data(:,1), movmean(data(:,2:end),100)];
hold on
plot(data(:,1),data(:,2),'--r')
plot(data(:,1),data(:,3),'--g')
plot(data(:,1),data(:,4),'--b')
legend('Anxiety','Depression','Insomnia',...
    'Location','best')
title('--- (males), -\,-\,- (females)')








