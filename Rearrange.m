startup
% Input: GenetypeV*.csv and PhenotypeV*.csv (V = version)
% Object: Rearrange, remove missing rows and categorize cohort by number
% Output: V*data.mat file containing X (SNPs), Y (Phenotypes) and IDs
%% Load
visit   = 1;
Xorg    = readtable(strcat('GenotypesV',num2str(visit),'.csv')); 
Yorg    = readtable(strcat('PhenotypesV',num2str(visit),'.csv')); 
Yorg.Cohorts(strcmp(Yorg.Cohort,'WSC'),1)   = 1;
Yorg.Cohorts(strcmp(Yorg.Cohort,'MrOS'),1)  = 2;
Yorg.Cohorts(strcmp(Yorg.Cohort,'SSC'),1)   = 3;
Yorg    = Yorg(:,2:end);
Xorg    = Xorg(:,3:end);
Z       = [Xorg,Yorg];
Z       = rmmissing(Z);
Z.AGE(Z.AGE>100) = Z.AGE(Z.AGE>100)/10000;
Z.TST = Z.TST/60;
Z.SOL = Z.SOL/60;
Z.REML = Z.REML/60;
Z(Z.BMI < 5,:) = [];
X   = Z(:,1:239);
Y   = Z(:,241:end);
IDs = Z(:,240);
save(strcat('V',num2str(visit),'data'),'X','Y','IDs');