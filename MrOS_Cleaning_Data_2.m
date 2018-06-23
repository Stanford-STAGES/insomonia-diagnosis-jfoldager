startup
load('mros_gen')
load('mros_psg')
ID = readtable('MROS_DATASET_WITH_HA_ID.csv');
% 
oldIDs  = ID.ID_Unique;
fieldd  = fields(FEATURES);
X       = ones(length(fieldd), size(mros_gen.snps,1)).*(-1);    
idx     = 1;
genIDs  = cellfun(@str2double,mros_gen.id);
mrosIDs = cell(size(X)); mrosIDs = {''};
for i = 1:length(oldIDs)
    if ~isnan(ID.HA_ID(i)) && sum(genIDs == ID.HA_ID(i)) == 1
        X(idx,:) = mros_gen.data_expected(genIDs == ID.HA_ID(i),:);
        mrosIDs{idx} = oldIDs{i};
        idx = idx + 1;
    end
end