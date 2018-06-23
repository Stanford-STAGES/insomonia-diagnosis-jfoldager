startup
MrOS = readtable('mros-visit1-dataset-0.3.0.csv',...
    detectImportOptions('mros-visit1-dataset-0.3.0.csv'));
names = MrOS.Properties.VariableNames;
nameIdxs = ismember(names,{'nsrrid','vsage1','gender','hwbmi','poai_all',...
    'epepwort' }); % PQPTMWAK
MrOS = MrOS(:,nameIdxs);
MrOS.Properties.VariableNames = {'ID', 'AGE','ESS','BMI','ARI' ...
    ,'SEX'};
MrOS = [MrOS(:,1),MrOS(:,sort(MrOS.Properties.VariableNames(2:end)))];
MrOS.SEX = MrOS.SEX == 2;
MrOS = rmmissing(MrOS,1);
Hyps    = readtable('MrOSHypnogramsV1.csv');
names = Hyps.Properties.VariableNames;
for i = 2:length(names), MrOS.(names{i}) = nan(size(MrOS,1),1);end

for ii = 1:size(MrOS,1)
    if sum(contains(Hyps.ID,MrOS.ID{ii})) == 1
for i = 2:length(names), MrOS.(names{i})(ii) = Hyps.(names{i})(contains(Hyps.ID,MrOS.ID{ii})); end
    end
end

MrOS = [MrOS(:,1),MrOS(:,sort(MrOS.Properties.VariableNames(2:end)))];
writetable(MrOS,'MrOSv1.csv')
MrOS2 = readtable('MrOSv2.csv');
MrOS2.INS = []; MrOS2.DEP = []; MrOS2.ANX = [];
MrOS2 = [MrOS2(:,1),MrOS2(:,sort(MrOS2.Properties.VariableNames(2:end)))];

for i = 1:size(MrOS,1)
    if sum(contains(MrOS2.ID,MrOS.ID{i})) == 1
        MrOS{i,2:end} = nanmean([MrOS{i,2:end};MrOS2{contains(MrOS2.ID,MrOS.ID{i}),2:end}]);
    end
end
writetable(MrOS,'MrOS.csv')