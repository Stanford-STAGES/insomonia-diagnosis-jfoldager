startup
% Input: csv files to be gathered
% Objective: cleaning up all data and extract selected attributes
% Output: WSC.csv which contains the objective
%%
WSCDem = readtable('plm_t.csv');
WSCObj = readtable('objectiveAndsubjective.csv');
WSCNEW = readtable('WSC_NEWEST.csv',...
    detectImportOptions('WSC_NEWEST.csv')); 
WSCNEW.VISIT = [];
WSCNEW.SEX = contains(WSCNEW.SEX,'M');
WSCARI = readtable('WSC_ARI.csv');
WSC = WSCObj;
% Extract only specific attributes
names = WSC.Properties.VariableNames;
nameIdxs = ismember(names,{'SUBJ_ID','anxiety','bmi','age','AHI4_ADJUSTED_V2',...
    'SEX','plmi_d140_sw','depressed','insom_at_lab','doze_sc'});
WSC = WSC(:,nameIdxs);
WSC.SUBJ_ID = strcat(WSC.SUBJ_ID,'_',num2str(WSCObj.VISIT_NUMBER));
WSC.Properties.VariableNames = upper(WSC.Properties.VariableNames);
WSC = [WSC(:,1),WSC(:,sort(WSC.Properties.VariableNames(2:end)))];
WSC.Properties.VariableNames = {'ID','AGE','AHI','ANX',...
    'BMI','DEP','ESS','INS','PLMI','SEX'};
ARIs = [];
WSC.ARI = nan(size(WSC,1),1);
for i = 1:size(WSCARI,1)
    if sum(ismember(WSC.ID,WSCARI.Record_ID{i})) == 1
        WSC.ARI(ismember(WSC.ID,WSCARI.Record_ID{i})) = ...
            WSCARI.ARI(i);
    end
end
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
pth                 = 'Data\wsc\psgs\';
listing             = extractfield(dir(pth), 'name')';
events              = listing(contains(lower(listing),'.sta'));
nSamples            = length(events);
T                   = table();
T.ID                = cell(nSamples,1);
T.TST               = nan(nSamples,1);
T.SOL               = nan(nSamples,1);
T.SE                = nan(nSamples,1);
T.AW25              = nan(nSamples,1);
T.AW5               = nan(nSamples,1);
T.AWN125            = nan(nSamples,1);
T.AWN15             = nan(nSamples,1);
T.WA                = nan(nSamples,1);
T.N1                = nan(nSamples,1);
T.N2                = nan(nSamples,1);
T.N3                = nan(nSamples,1);
T.REM               = nan(nSamples,1);
for i = 1:nSamples
    ii      = split(events{i}, ' ');
    T.ID{i} = ii{1};
    fileID  = fopen(strcat(pth,events{i}));
        STA     = textscan(fileID, '%f%f%f', 'Delimiter', '\t', 'EmptyValue', -1);
    fclose(fileID);
    hypnogram           = STA{2}; hypnogram(hypnogram==4) = 3;
    lightsOffEpoch      = find(diff(hypnogram == 7) == -1); 
    if isempty(lightsOffEpoch), lightsOffEpoch = 1; else, lightsOffEpoch = lightsOffEpoch(1)+1; end
    if hypnogram(end,1) == 7
        lightsOnEpoch   = find(diff(hypnogram == 7) == 1); 
        lightsOnEpoch   = lightsOnEpoch(end);
    else
        lightsOnEpoch   = length(hypnogram);
    end
    hypnogram        = hypnogram(lightsOffEpoch:lightsOnEpoch);
    SOP            = find((hypnogram>0).*(hypnogram<7) == 1); SOP=SOP(1);
    sleephyp        = hypnogram(SOP:end);
    T.TST(i)       = sum(hypnogram>0)*30;
    T.SOL(i)       = (SOP-1)*30;
    T.SE(i)        = sum((sleephyp > 0).*(sleephyp<7) == 1)/length(sleephyp);
    T.WA(i)        = sum(sleephyp == 0)/length(sleephyp);
    T.N1(i)        = sum(sleephyp == 1)/length(sleephyp);
    T.N2(i)        = sum(sleephyp == 2)/length(sleephyp);
    T.N3(i)        = sum(sleephyp == 3)/length(sleephyp);
    T.REM(i)       = sum(sleephyp == 5)/length(sleephyp);
    try a  = find(sleephyp == 5); T.REML(i) = (a(1)-1)*30; catch, end
    a = diff(find(diff(diff(cumsum((sleephyp == 0) + (sleephyp == 7)))) ~= 0));
    a = a(1:2:end);
    T.AW25(i) = sum(a>=(2.5*60)/30)/(T.TST(i)/60/60);
    T.AW5(i)  = sum(a>=(5*60)/30)/(T.TST(i)/60/60);
    a = diff(find(diff(diff(cumsum((sleephyp == 0) + (sleephyp == 1) + (sleephyp == 7)))) ~= 0));
    a = a(2:2:end);
    T.AWN125(i)  = sum(a>=(2.5*60)/30)/(T.TST(i)/60/60);
    T.AWN15(i)    = sum(a>=(5*60)/30)/(T.TST(i)/60/60);      
    fprintf('%i out of %i is done.\n',i,nSamples);
end

WSCObj =  T;
names = WSCObj.Properties.VariableNames;
for i = 2:length(names), WSC.(names{i}) = nan(size(WSC,1),1);end

for ii = 1:size(WSC,1)
    if sum(contains(WSCObj.ID,WSC.ID{ii})) == 1
for i = 2:length(names), WSC.(names{i})(ii) = WSCObj.(names{i})(contains(WSCObj.ID,WSC.ID{ii})); end
    end
end

WSC.SEX = strcmp(WSC.SEX,'M');
WSC = [WSC(:,1),WSC(:,sort(WSC.Properties.VariableNames(2:end)))];
writetable(WSC,'WSCall.csv')


