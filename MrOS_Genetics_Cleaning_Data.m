startup
extDiskGenPath  = 'mros/genetics/';
listing         = extractfield(dir(extDiskGenPath), 'name')';
genListing      = listing(contains(listing,'.gen'));
sampleListing   = listing(contains(listing,'.sample'));
%% Raw imputation
nColumnOffset   = 5;
% 500k sampleListing{1}
fileID          = fopen(strcat(extDiskGenPath,sampleListing{1}),'r');
C500               = textscan(fileID,'%s %s %s %s %s %s');
                fclose(fileID);
nSubjects(1)    = length(C500{1})-2;      
fileID          = fopen(strcat(extDiskGenPath,genListing{1}),'r');
C500               = textscan(fileID, repmat('%s ',1,3*nSubjects(1)+5));
                fclose(fileID);  
nSNPs(1)        = length(C500{1,2});
SNPs.S1         = C500{1,2};
spiller = nSubjects(1) == (length(C500)-nColumnOffset)/3;             

DATA500k = ones(nSubjects(1),length(C500{1})*3)*(-1);
idxx = 1;
for sub = 1:nSubjects(1)
    i = 1;
    idx = 1;
    while idx <= length(C500{1})
        DATA500k(sub,i:i+2) = [str2double(C500{idxx+5}(idx)), str2double(C500{idxx+6}(idx)), str2double(C500{idxx+7}(idx))];
        idx = idx + 1;
        i = i + 3;
    end
    idxx = idxx + 3;
end
%%
[~, idTable, ~] = xlsread(strcat(extDiskIDPath,idsListing{end}));
fileID500k      = fopen(strcat(extDiskGenPath,sampleListing{1}),'r');
sample500k      = textscan(fileID500k,'%s %s %s %s %s');
                fclose(fileID500k);
%% 
usedIDs         = repmat({''},length(edfListing),1);
data            = ones(sum(nSubjects),size(DATA500k,2))*(-1);
sub             = 1;
fun             = 1;
for i = 1:length(edfListing)
    curID       = split(edfListing{i},'_'); curID = curID{1};
    if sum(contains(usedIDs,curID)) == 1, continue, end
    if sum(contains(idTable(:,2),curID)) ~= 1, funny{fun}  = curID; fun = fun + 1; continue, end
    curOldID    = idTable{contains(idTable(:,2),curID),1};
    curImput    = idTable{contains(idTable(:,2),curID),3};
    dataPoint   = ones(1,size(DATA500k,2))*(-1);
    if contains(curImput,'affy6') && sum(contains(sample6k{1,2},curOldID)) == 1
        posID       = find(contains(sample6k{1,2},curOldID)) - 2;
        data(sub,:) = DATA6k(posID,:);
        usedIDs{sub}= curID;
        sub         = sub + 1;
    elseif contains(curImput,'affy500k') && sum(contains(sample500k{1,2},curOldID)) == 1
        posID = find(contains(sample500k{1,2},curOldID)) - 2;
        data(sub,:) = DATA500k(posID,:);
        usedIDs{sub}= curID;
        sub         = sub + 1;
    else
        funny{fun}  = curID; fun = fun + 1;
    end
end
usedIDs(strcmp('',usedIDs)) = [];
data( any(data==-1,2), : ) = []; 
%%


Y = repmat([0,1,2],size(data,1),size(data,2)/3);
YY = Y.*data;
YYY = zeros(size(YY,1),size(YY,2)/3);
for i = 1:size(YY,1)
    idx = 1;
    idxx = 1;
    while idx <= size(YY,2)-2
       YYY(i,idxx)=sum(YY(i,idx:idx+2));
       idx = idx + 3;
       idxx = idxx + 1;
    end
end

wsc_gen.data_expected = YYY;
save('wsc_gen','wsc_gen')

% [coeff,score,latent] = pca(zscore(wsc_gen.data_expected));
% scatter(score(:,1),score(:,2))
% stem(latent/sum(latent))
% biplot(coeff(:,1:2),'varlabels',wsc_gen.snps(:,1))




