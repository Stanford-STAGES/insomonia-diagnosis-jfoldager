startup
MrOSv2 = readtable('mros-visit2-dataset-0.3.0.csv',...
    detectImportOptions('mros-visit2-dataset-0.3.0.csv'));
names = MrOSv2.Properties.VariableNames;
nameIdxs = ismember(names,{'nsrrid','vs2age1','gender','hwbmi','poai_all',...
    'slisiscr','epepwort' ...
%   ,'poxslpmn' ,'potmremp', 'pqpeffcy','potmst1p','potmst1p','potmst2p','potms34p','poslprdp' ...
    });
MrOSv2 = MrOSv2(:,nameIdxs);
MrOSv2 = [MrOSv2(:,1),MrOSv2(:,sort(MrOSv2.Properties.VariableNames(2:end)))];
MrOSv2.Properties.VariableNames = {'ID','ESS','SEX','BMI','ARI' ...
    'ISI','AGE'};
%     'TST','REM','N3','N1','N2','SRST','SE',...
MrOSv2 = [MrOSv2(:,1),MrOSv2(:,sort(MrOSv2.Properties.VariableNames(2:end)))];
MrOSv2.SEX = MrOSv2.SEX == 2;
MrOSv2.ISI = MrOSv2.ISI >= 10;
MrOSv2.Properties.VariableNames{13} = 'INS';
MrOSv2 = rmmissing(MrOSv2,1);
Hyps    = readtable('MrOSHypnograms.csv');
names = Hyps.Properties.VariableNames;
for i = 2:length(names), MrOSv2.(names{i}) = nan(size(MrOSv2,1),1);end

for ii = 1:size(MrOSv2,1)
    if sum(contains(Hyps.ID,MrOSv2.ID{ii})) == 1
for i = 2:length(names), MrOSv2.(names{i})(ii) = Hyps.(names{i})(contains(Hyps.ID,MrOSv2.ID{ii})); end
    end
end

Hyps = readtable('mros\psgv2\GOLDBERG_VARS_MROS.csv');
Hyps.AXANXSC = [];
Hyps.AXDEPSC = [];
Hyps.Properties.VariableNames{1} = 'ANX';
Hyps.Properties.VariableNames{2} = 'DEP';
Hyps.Properties.VariableNames{3} = 'ID';
names = Hyps.Properties.VariableNames;
for i = 1:length(names)-1, MrOSv2.(names{i}) = nan(size(MrOSv2,1),1);end

for ii = 1:size(MrOSv2,1)
    if sum(contains(Hyps.ID,MrOSv2.ID{ii})) == 1
for i = 1:length(names)-1, MrOSv2.(names{i})(ii) = Hyps.(names{i})(contains(Hyps.ID,MrOSv2.ID{ii})); end
    end
end
MrOSv2 = [MrOSv2(:,1),MrOSv2(:,sort(MrOSv2.Properties.VariableNames(2:end)))];
writetable(MrOSv2,'MrOS.csv')
