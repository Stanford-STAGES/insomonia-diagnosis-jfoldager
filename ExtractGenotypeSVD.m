startup
% Input: csv of subjects (rows) and principal components (columns)
% Output: PCA ordered by dbID
%% 
ALL = readtable('All4_MDS.csv');
ALL.None1 = [];
ALL.None2 = [];
pcas = ALL{:,4:end};
warning('off','all');
wsc = readtable('WSC_NEWEST.csv',...
    detectImportOptions('WSC_NEWEST.csv'));
mros = readtable('MrOS_IDRosettaStone_NSRR-HAdbGaP-MrOSID.xlsx');
ssc = readtable('SSC.csv');  ssc.ID = string(ssc.ID);
load ssc_gen.mat
mros = rmmissing(mros); mros.HA_ID = string(mros.HA_ID);
m = 1;
w = 1;
s = 1;
SSC = table();
WSC = table();
MROS = table();
for i = 1:size(ALL,1)
    if all(ismember(ALL.IID{i}, '0':'9')) && ...
            sum(strcmp(mros.HA_ID,ALL.IID{i})) == 1
        MROS.ID{m} = upper(mros.NSRR_ID{strcmp(mros.HA_ID,ALL.IID{i})});
        MROS.PCs(m,:) = ALL{i,4:end};
        m = m + 1;
    elseif sum(contains(lower(ALL.IID{i}), {'na','wisc'})) > 0
        nums = regexp(ALL.IID{i},'\d');
        gwas = ALL.IID{i}(nums(1):nums(end));
        if sum(strcmp(wsc.GWAS,gwas)) == 1
            WSC.ID{w} = upper(wsc.ID{strcmp(wsc.GWAS,gwas)});
            WSC.PCs(w,:) = ALL{i,4:end};
            w = w + 1;
        end
    elseif any(ismember(lower(ALL.FID{i}), 'a':'z')) && ...
            ~startsWith(ALL.IID{i},'NA')
        ids = split(ALL.FID{i},'_'); id = ids{end};
        if sum(strcmp(ssc.GWAS,id)) == 1
            SSC.ID{s} = upper(ssc.ID{strcmp(ssc.GWAS,id)});
            SSC.PCs(s,:) = ALL{i,4:end};
            s = s + 1;
        end
    end
end
MROS = sortrows(MROS); mros = MROS{:,2};
WSC = sortrows(WSC); wsc = WSC{:,2};
SSC = sortrows(SSC); ssc = SSC{:,2};
all = [WSC;MROS;SSC];
grps = [ones(size(mros,1),1);2*ones(size(wsc,1),1);3*ones(size(ssc,1),1)];
figure,gscatter(all(:,1),all(:,2),grps);
xlabel('PC 1'), ylabel('PC 2')
legend({'MrOS','WSC','SSC'})
figure,gscatter(all(:,2),all(:,3),grps);
xlabel('PC 2'), ylabel('PC 3')
legend({'MrOS','WSC','SSC'})
writetable(all,'PCA1to10.csv');