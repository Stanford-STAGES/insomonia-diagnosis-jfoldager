startup
extDiskPath     = 'mros/';
extDiskGenPath  = 'mros/genetics/';
listing         = extractfield(dir(extDiskGenPath), 'name')';
genListing      = listing(contains(listing,'.gen'));
sampleListing   = listing(contains(listing,'.sample'));

fileID          = fopen(strcat(extDiskGenPath,sampleListing{1}),'r');
Csample         = textscan(fileID,'%s %s %s %s %s %s');
                fclose(fileID);
nSubjects       = length(Csample{1})-2;      
fileID          = fopen(strcat(extDiskGenPath,genListing{1}),'r');
C               = textscan(fileID, repmat('%s ',1,3*nSubjects+5));
                fclose(fileID);   
                
DATA = ones(nSubjects,length(C{1})*3)*(-1);
idxx = 1;
for sub = 1:nSubjects
    i = 1;
    idx = 1;
    while idx <= length(C{1})
        DATA(sub,i:i+2) = [str2double(C{idxx+5}(idx)), str2double(C{idxx+6}(idx)), str2double(C{idxx+7}(idx))];
        idx = idx + 1;
        i = i + 3;
    end
    idxx = idxx + 3;
end

DATA(any(isnan(DATA), 2), :) = [];


Y = repmat([2,1,0],size(DATA,1),size(DATA,2)/3);
YY = Y.*DATA;
for i = 1:size(YY,1)%YYY = zeros(size(YY,1),size(YY,2)/3);

    idx = 1;
    idxx = 1;
    while idx <= size(YY,2)-2
       YYY(i,idxx)=sum(YY(i,idx:idx+2));
       idx = idx + 3;
       idxx = idxx + 1;
    end
end

%% Store in model
mros_gen.id             = Csample{1,1}(3:end);
mros_gen.snps           = [C{1,1},C{1,2}];
mros_gen.data           = DATA;
mros_gen.data_expected  = YYY;
%% Save variable
save(strcat(extDiskPath,'mros_gen'),'mros_gen')
save(strcat('Data/MrOS/mros_gen'),'mros_gen')
