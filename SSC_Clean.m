startup
pth1 = 'Data\ssc\BMI_ESS_INS_DEP_ANX.csv';
pth2 = 'Data\ssc\GWAS_SEX_AGE.csv';
%% Load relevant tables
SSC1  =   readtable(pth1); %length(unique(SSC1.DbID)) == length(SSC1.DbID)
SSC2  =   readtable(pth2); %length(unique(SSC2.DbID)) == length(SSC2.DbID)
SSC2.GWAS = strrep(SSC2.GWAS,'.','');
SSC = table();
for i = 1:size(SSC1,1)
    if ismember(SSC1.DbID(i),SSC2.DbID) == 1
        SSC = [SSC;SSC1(i,:),SSC2(SSC2.DbID==SSC1.DbID(i),2:end)];
    end
end
SSC.Height = [];
SSC.Weight = [];
SSC.BMI = str2double(SSC.BMI);
SSC = rmmissing(SSC,1);
SSC.BMI = round(SSC.BMI,2);
SSC.AGE = round(SSC.AGE,2);
SSC.Properties.VariableNames{1} = 'ID';
%% Hypnograms
pth1 = 'Data\ssc\Events.csv';
pth2 = 'Data\ssc\Hypnograms.csv';
SSCObj =   readtable(pth2);
names = SSCObj.Properties.VariableNames;
for i = 2:length(names), SSC.(names{i}) = nan(size(SSC,1),1); end
for i = 1:size(SSCObj,1)
    id = SSCObj.ID(i);
    if sum(SSC.ID==id) == 1
        for ii = 2:length(names)
            SSC.(names{ii})(SSC.ID==id) = SSCObj.(names{ii})(i); 
        end
    end
end
%% Events
SSCObj      = readtable(pth1);
SSC.ARI     = nan(size(SSC,1),1);
SSC.AHI     = nan(size(SSC,1),1);
SSC.PLMI    = nan(size(SSC,1),1);
for i = 1:size(SSCObj,1)
    id = SSCObj.ID(i);
    if sum(SSC.ID==id) == 1
        SSC.ARI(SSC.ID==id) = SSCObj.ARI(i);
        SSC.ARI(SSC.ID==id) = SSC.ARI(SSC.ID==id)/(SSC.TST(SSC.ID==id)/60/60);
        SSC.AHI(SSC.ID==id) = SSCObj.AHI(i);
        SSC.AHI(SSC.ID==id) = SSC.AHI(SSC.ID==id)/(SSC.TST(SSC.ID==id)/60/60);
        SSC.PLMI(SSC.ID==id) = SSCObj.PLMI(i);
        SSC.PLMI(SSC.ID==id) = SSC.PLMI(SSC.ID==id)/(SSC.TST(SSC.ID==id)/60/60);
    end
end
SSC         = rmmissing(SSC,1);
SSC = [SSC(:,1),SSC(:,sort(SSC.Properties.VariableNames(2:end)))];
%% GWAS
pth2    = 'Data\ssc\genetics\IDConversionTable.csv';
SSCObj  =  readtable(pth2);
SSC2    = table();
for i = 1:size(SSC,1)
    if sum(ismember(SSCObj.DbID,SSC.ID(i))) == 1
        SSC2 = [SSC2;SSC(i,:)];
        SSC2.GWAS{end} = SSCObj.GWASLocation{i};
    end
end
SSC2.SEX = strcmp(SSC2.SEX,'M');
writetable(SSC2,'SSC.csv');
