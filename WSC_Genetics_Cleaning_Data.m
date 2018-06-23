startup
extDiskPSGPath  = 'Data/wsc/psg/';
extDiskMetaPath = 'Data/wsc/_meta/';
extDiskGenPath  = 'Data/wsc/genetics/';
extDiskIDPath   = 'Data/wsc/genetics/id_maps/';
listing         = extractfield(dir(extDiskGenPath), 'name')';
genListing      = listing(contains(listing,'.gen'));
sampleListing   = listing(contains(listing,'.sample'));
listing         = extractfield(dir(extDiskPSGPath), 'name')';
edfListing      = listing(contains(lower(listing),'.edf'));
listing         = extractfield(dir(extDiskIDPath), 'name')';
idsListing      = listing(contains(lower(listing),'.xlsx'));
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

% 6k sampleListing{2}
fileID          = fopen(strcat(extDiskGenPath,sampleListing{2}),'r');
C6               = textscan(fileID,'%s %s %s %s %s %s');
                fclose(fileID);
nSubjects(2)    = length(C6{1})-1;      
fileID          = fopen(strcat(extDiskGenPath,genListing{2}),'r');
C6               = textscan(fileID, repmat('%s ',1,3*nSubjects(2)+5));
                fclose(fileID);   
spiller = nSubjects(2) == (length(C6)-nColumnOffset)/3;      
nSNPs(2)        = length(C6{1,2});
SNPs.S2         = C6{1,2};
% if sum(diff(nSNPs)) ~= 0
contains(SNPs.S1,SNPs.S2);
discardedSPNs = [C500{1,1}(contains(SNPs.S1,SNPs.S2)==0)...
    C500{1,2}(contains(SNPs.S1,SNPs.S2)==0)];
for i = 1:length(C500)
    C500{i} = C500{i}(contains(SNPs.S1,SNPs.S2));
end

% end

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

DATA6k = ones(nSubjects(2),length(C6{1})*3)*(-1);
idxx = 1;
for sub = 1:nSubjects(2)
    i = 1;
    idx = 1;
    while idx <= length(C6{1})
        DATA6k(sub,i:i+2) = [str2double(C6{idxx+5}(idx)), str2double(C6{idxx+6}(idx)), str2double(C6{idxx+7}(idx))];
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
fileID6k        = fopen(strcat(extDiskGenPath,sampleListing{2}),'r');
sample6k        = textscan(fileID6k,'%s %s %s %s %s');
                fclose(fileID6k);
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
wsc_gen.id = usedIDs;
wsc_gen.snps = [C500{1},SNPs.S2];
wsc_gen.data = data;

Y = repmat([2,1,0],size(data,1),size(data,2)/3);
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
save(strcat('Data/WSC/wsc_gen_new'),'wsc_gen')

% [coeff,score,latent] = pca(zscore(wsc_gen.data_expected));
% scatter(score(:,1),score(:,2))
% stem(latent/sum(latent))
% biplot(coeff(:,1:2),'varlabels',wsc_gen.snps(:,1))




