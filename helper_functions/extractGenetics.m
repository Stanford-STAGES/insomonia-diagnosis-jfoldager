function out = extractGenetics(genpath,savepath)
% startup
extDiskPath     = savepath;
extDiskGenPath  = genpath;
listing         = extractfield(dir(extDiskGenPath), 'name')';
genListing      = listing(contains(listing,'.gen'));
sampleListing   = listing(contains(listing,'.sample'));
%% Get sample and imputation files
fileID          = fopen(strcat(extDiskGenPath,sampleListing{1}),'r');
Csample         = textscan(fileID,'%s %s %s %s %s %s');
                fclose(fileID);
nSubjects       = length(Csample{1})-2;      
fileID          = fopen(strcat(extDiskGenPath,genListing{1}),'r');
C               = textscan(fileID, repmat('%s ',1,3*nSubjects+5));
                fclose(fileID);   
%% Extract probabilities
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
% Remove empty rows if a given subject lacked one SNP information
DATA(any(isnan(DATA), 2), :) = [];
%% Calculate expected dosage from probabilities
multiplicationMatrix = repmat([2,1,0],size(DATA,1),size(DATA,2)/3);
multiplicationMatrix = multiplicationMatrix.*DATA;
expectedDosage       = ones(nSubjects,length(C{1}))*(-1);
for i = 1:size(multiplicationMatrix,1)
    idx = 1;
    idxx = 1;
    while idx <= size(multiplicationMatrix,2)-2
       expectedDosage(i,idxx)   = sum(multiplicationMatrix(i,idx:idx+2));
       idx                      = idx + 3; 
       idxx                     = idxx + 1;
    end
end
%% Store in model
gen.id             = Csample{1,1}(3:end);
gen.snps           = [C{1,1},C{1,2}];
gen.data           = DATA;
gen.data_expected  = expectedDosage;
%% Save variable
save(strcat(extDiskPath,'ssc_gen'),'gen')
save(strcat('Data/SSC/ssc_gen'),'gen')
end