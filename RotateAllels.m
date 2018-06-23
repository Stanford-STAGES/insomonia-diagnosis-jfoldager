startup
%% 
% Input: SNP data from all cohorts and table from Jansen showing effect
% allele directions.
% Output: correct effect direction to match Jansen in our used SNP data.
%%
T = readtable('AllData.csv');
Jansen      = readtable('SNPsJansen.csv');
nSubjects   = 2436;      
fileID      = fopen('mros\genetics\Imputation.csv','r');
dataLogan   = textscan(fileID, repmat('%s ',1,3*nSubjects+5));
            fclose(fileID);   
SNPs        = T(:,1:239);
SNPscopy    = SNPs;
snps        = dataLogan{1,2};
A1Logan     = dataLogan{1,4};
A2Logan     = dataLogan{1,5};
rsIDsLogan  = SNPs.Properties.VariableNames;
rsIDsJansen = Jansen.rsID;
shit        = {};
exclude     = [];
rotate      = [];
dircs       = [];
newestSNPs  = table();
for i = 1:length(rsIDsLogan)
    curSNP = rsIDsLogan{i};
    if sum(ismember(snps,curSNP)) == 1 && ...
           sum(ismember(rsIDsJansen,curSNP)) == 1 
       a1jansen = Jansen{ismember(rsIDsJansen,curSNP),5}; 
       a2jansen = Jansen{ismember(rsIDsJansen,curSNP),6};
       a1logan  = A1Logan{ismember(snps,curSNP)};
       a2logan  = A2Logan{ismember(snps,curSNP)};

       if strcmp(a1jansen,a1logan) && strcmp(a2jansen,a2logan)
            dircs    = [dircs; Jansen{ismember(rsIDsJansen,curSNP),7}>0];
            rotate   = [rotate;false];
       elseif strcmp(a1jansen,a2logan) && strcmp(a2jansen,a1logan)
            rotate   = [rotate;true];
            dircs    = [dircs; Jansen{ismember(rsIDsJansen,curSNP),7}>0];
       else
           shit = [shit;curSNP];
           exclude = [exclude;i];
       end
    else
        exclude = [exclude;i];
    end
end

%%
dataV2 = load('V2data.mat');
dataV2.X(:,exclude) = [];
for i = 1:size(dataV2.X,2)
    if rotate(i) == 1
        dataV2.X{:,i} = 2 - dataV2.X{:,i};
    end
end
snp2gene       = readtable('SNP2GeneNew.csv');
finalSNPsTable = table();
nearestGenes    = table();
for i = 1:size(dataV2.X,2)
    finalSNPsTable = [finalSNPsTable;...
        Jansen(strcmp(Jansen.rsID,dataV2.X.Properties.VariableNames{i}),[2:6,9:10])];
    nearestGenes.NearestGene{i,1} = snp2gene.NearestGene{...
        strcmp(snp2gene.rsID,dataV2.X.Properties.VariableNames{i})};
end
writetable([finalSNPsTable,nearestGenes],'FinalSNPs.csv');
% Xav2 = table2array(dataV2.X);
% MAFv2 = (sum(round(Xav2) == 0)+(sum(round(Xav2) == 1)/2))./...
%     (size(Xav2,1));
% dataV2.X(:,MAFv2<0.05) = [];
Xav2 = table2array(dataV2.X);

dataV1 = load('V1data.mat');
dataV1.X(:,exclude) = [];
for i = 1:size(dataV1.X,2)
    if rotate(i) == 1
        dataV1.X{:,i} = 2 - dataV1.X{:,i};
    end
end
% Xav1 = table2array(dataV1.X);
% MAFv1 = (sum(round(Xav1) == 0)+(sum(round(Xav1) == 1)/2))./...
%     (size(Xav1,1));
% dataV1.X(:,MAFv2<0.05) = [];
Xav1 = table2array(dataV1.X);

X = dataV1.X;
Y = dataV1.Y;
IDs = dataV1.IDs;
save('V1newdata.mat','X','Y','IDs');
X = dataV2.X;
Y = dataV2.Y;
IDs = dataV2.IDs;
save('V2newdata.mat','X','Y','IDs');

save('JansenDirections.mat','dircs','rotate');


