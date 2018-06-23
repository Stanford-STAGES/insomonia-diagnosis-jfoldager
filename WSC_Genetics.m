startup
extDiskGenPath  = 'Data/wsc/genetics/';
listing         = extractfield(dir(extDiskGenPath), 'name')';
genListing      = listing(contains(listing,'.gen'));
sampleListing   = listing(contains(listing,'.sample'));

fileID          = fopen(strcat(extDiskGenPath,sampleListing{1}),'r');
C               = textscan(fileID,'%s %s %s %s %s %s');
                fclose(fileID);
nSubjects       = length(C{1})-1;      
fileID          = fopen(strcat(extDiskGenPath,genListing{2}),'r');
C               = textscan(fileID, repmat('%s ',1,3*nSubjects+5));
                fclose(fileID);   
                
DATA = ones((length(C)-5)/3-2,length(C{1})*3)*(-1);
idxx = 1;
for sub = 1:(length(C)-5)/3-2 
    i = 1;
    idx = 1;
    while idx <= length(C{1})
        DATA(sub,i:i+2) = [str2double(C{idxx+5}(idx)), str2double(C{idxx+6}(idx)), str2double(C{idxx+7}(idx))];
        idx = idx + 1;
        i = i + 3;
    end
    idxx = idxx + 3;
end
[~,score,latent] = pca(zscore(DATA));
scatter3(score(:,1),score(:,2),score(:,3))
stem(latent/sum(latent))

Y = repmat([2,1,0],size(DATA,1),size(DATA,2)/3);
YY = Y.*DATA;
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

[coeff,score,latent] = pca(zscore(YYY));
scatter(score(:,1),score(:,2))
stem(latent/sum(latent))
biplot(coeff(:,1:2),'varlabels',C{1,1})

YYYY = tsne(YYY);
scatter(YYYY(:,1),YYYY(:,2))