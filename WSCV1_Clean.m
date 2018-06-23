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
% Extract only first visit
usedids = [];
INSANXDEP = [];
for i = 1:size(WSCObj,1)
    id = WSCObj.SUBJ_ID{i};
    if sum(ismember(usedids,id)) == 0
        INSANXDEP = [INSANXDEP;WSCObj(i,:)];
        usedids = [usedids;{id}];
    end
end
% Extract only specific attributes
names = INSANXDEP.Properties.VariableNames;
nameIdxs = ismember(names,{'SUBJ_ID','anxiety','bmi','depressed','insom_at_lab','doze_sc'});
INSANXDEP = INSANXDEP(:,nameIdxs);
INSANXDEP = rmmissing(INSANXDEP);
% Extract only first visit
usedids = [];
INS = [];
WSCNEW = rmmissing(WSCNEW);
for i = 1:size(WSCNEW,1)
    id = WSCNEW.ID{i};
    if sum(ismember(usedids,id)) == 0
        INS = [INS;WSCNEW(i,:)];
        usedids = [usedids;{id}];
    end
end

usedids = [];
ARIs = [];
for i = 1:size(WSCARI,1)
    id = split(WSCARI.Record_ID{i},'_'); id = id{1};
    if sum(ismember(usedids,id)) == 0
        ARIs = [ARIs;WSCARI(i,:)];
        usedids = [usedids;{id}];
    end
end

usedids = [];
WSCNEWW = [];
for i = 1:size(WSCNEW,1)
    id = WSCNEW.ID{i};
    if sum(ismember(usedids,id)) == 0
        WSCNEWW = [WSCNEWW;WSCNEW(i,:)];
        usedids = [usedids;{id}];
    end
end
WSCNEWW.ANX = nan(size(WSCNEWW,1),1);
WSCNEWW.DEP = nan(size(WSCNEWW,1),1);
for i = 1:length(INSANXDEP.SUBJ_ID)
    id = INSANXDEP.SUBJ_ID{i};
    if sum(contains(WSCNEWW.ID,id)) == 1
        WSCNEWW.ANX(contains(WSCNEWW.ID,id)) = INSANXDEP.anxiety(i);      
        WSCNEWW.DEP(contains(WSCNEWW.ID,id)) = INSANXDEP.depressed(i);     
        WSCNEWW.ESS(contains(WSCNEWW.ID,id)) = INSANXDEP.doze_sc(i);      
    end
end

WSCNEWW.ARI = nan(size(WSCNEWW,1),1);
for i = 1:length(ARIs.Record_ID)
    id = regexprep(ARIs.Record_ID{i},'_\d','');
    if sum(contains(WSCNEWW.ID,id)) == 1
        WSCNEWW.ARI(contains(WSCNEWW.ID,id)) = ARIs.ARI(i);        
    end
end

WSC = WSCNEWW;
WSC.GWAS = [];
WSC.TST  = []; % gets calculated later
pth = 'wsc\hypnograms\';
listing             = extractfield(dir(pth), 'name')';
evtListing          = listing(contains(lower(listing),'.sta'));
ids                 = split(evtListing,{' ','_'}); ids = ids(:,1);
[ids,ia,ic]         = unique(ids,'stable');
events              = evtListing(ia);
nSubjects           = length(events);
T                   = table();
T.ID                = cell(nSubjects,1);
T.TST               = nan(nSubjects,1);
T.SOL               = nan(nSubjects,1);
T.SE                = nan(nSubjects,1);
T.AW25              = nan(nSubjects,1);
T.AW5               = nan(nSubjects,1);
T.AWN125            = nan(nSubjects,1);
T.AWN15             = nan(nSubjects,1);
T.WA                = nan(nSubjects,1);
T.N1                = nan(nSubjects,1);
T.N2                = nan(nSubjects,1);
T.N3                = nan(nSubjects,1);
T.REM               = nan(nSubjects,1);
for i = 1:length(ids)
    curSubj = i;
    T.ID{i} = ids{i};
    fileID  = fopen(strcat(pth,evtListing{i}));
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
    T.TST(curSubj)       = sum(hypnogram>0)*30;
    T.SOL(curSubj)       = (SOP-1)*30;
    T.SE(curSubj)        = sum((sleephyp > 0).*(sleephyp<7) == 1)/length(sleephyp);
    T.WA(curSubj)        = sum(sleephyp == 0)/length(sleephyp);
    T.N1(curSubj)        = sum(sleephyp == 1)/length(sleephyp);
    T.N2(curSubj)        = sum(sleephyp == 2)/length(sleephyp);
    T.N3(curSubj)        = sum(sleephyp == 3)/length(sleephyp);
    T.REM(curSubj)       = sum(sleephyp == 5)/length(sleephyp);
    try a  = find(sleephyp == 5); T.REML(curSubj) = (a(1)-1)*30; catch, end
    a = diff(find(diff(diff(cumsum((sleephyp == 0) + (sleephyp == 7)))) ~= 0));
    a = a(1:2:end);
    T.AW25(curSubj) = sum(a>=(2.5*60)/30)/(T.TST(curSubj)/60/60);
    T.AW5(curSubj)  = sum(a>=(5*60)/30)/(T.TST(curSubj)/60/60);
    a = diff(find(diff(diff(cumsum((sleephyp == 0) + (sleephyp == 1) + (sleephyp == 7)))) ~= 0));
    a = a(2:2:end);
    T.AWN125(curSubj)  = sum(a>=(2.5*60)/30)/(T.TST(curSubj)/60/60);
    T.AWN15(curSubj)    = sum(a>=(5*60)/30)/(T.TST(curSubj)/60/60);      
    fprintf('%i out of %i is done.\n',curSubj,nSubjects);
end

WSCObj =  T;
names = WSCObj.Properties.VariableNames;
for i = 2:length(names), WSC.(names{i}) = nan(size(WSC,1),1);end

for ii = 1:size(WSC,1)
    if sum(contains(WSCObj.ID,WSC.ID{ii})) == 1
for i = 2:length(names), WSC.(names{i})(ii) = WSCObj.(names{i})(contains(WSCObj.ID,WSC.ID{ii})); end
    end
end



v2              = readtable('Data/wsc/data_all5_11_18.csv');
[ids,ia,ic]     = unique(v2.SUBJ_ID,'stable');
v2              = v2(ia,:);
names  = {'ID','AHI4_ADJUSTED_V2','doze_sc','bmi','SEX','age',...
    'insom_at_lab','anxiety','depressed'};
v2data = v2(:,ismember(v2.Properties.VariableNames,names));
v2data = [v2data(:,2),v2data(:,sort([v2data.Properties.VariableNames(1),v2data.Properties.VariableNames(3:end)]))];
v2data.Properties.VariableNames{2} = 'AHI'; v2data.AHI = round(str2double(v2data.AHI),1);
v2data.Properties.VariableNames{4} = 'AGE'; v2data.AGE = round(str2double(v2data.AGE),1);
v2data.Properties.VariableNames{5} = 'ANX'; v2data.ANX = str2double(v2data.ANX)==1;
v2data.Properties.VariableNames{6} = 'BMI'; v2data.BMI = round(str2double(v2data.BMI),1);
v2data.Properties.VariableNames{7} = 'DEP'; v2data.DEP = str2double(v2data.DEP)==1;
v2data.Properties.VariableNames{8} = 'ESS'; v2data.ESS = round(str2double(v2data.ESS),1);
v2data.Properties.VariableNames{9} = 'INS'; v2data.INS = str2double(v2data.INS)==1;
v2data.SEX = strcmp(v2data.SEX,'M');


WSC = [WSC(:,1),WSC(:,sort(WSC.Properties.VariableNames(2:end)))];
writetable(WSC,'WSCv1.csv')

