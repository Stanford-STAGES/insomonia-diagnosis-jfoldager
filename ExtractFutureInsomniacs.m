startup
T1 = readtable('WSC.csv');
Ids = split(T1.ID,'_'); Ids = Ids(:,1);
T1.TST     = [];          
T1.WA     = [];          
T1.N1     = [];          
T1.N2     = [];          
T1.N3     = [];          
T1.REM     = [];         
T1.REML     = [];        
T1.SOL     = [];        
T1.SE     = [];        
T1.AW25     = [];        
T1.AW5     = [];        
T1.AWN125     = [];        
T1.AWN15     = [];
T2 = readtable('WSChypnogram.csv');T2(:,2:11) = [];
T  = [];
for i = 1:size(T1,1)
   if sum(strcmp(T2.ID,T1.ID{i})) == 1
       T = [T;T1(i,:),T2(strcmp(T2.ID,T1.ID{i}),2:end)];
   end
end
T(:,contains(T.Properties.VariableNames,'PC_')) = [];
[~,~,ic] = unique(Ids,'stable');
BeforeInsomnia  = table();
AfterInsomnia   = table();
Controlv1         = table();
Controlv2         = table();
for i = 1:length(ic)
    if sum(ic == i) > 1 && ...
            sum(T.INS(ic == i)) ~= sum(ic == i) && ...
                    sum(T.INS(ic == i)) > 0 
        t = T(ic == i,:);
        befInsIdxs      = diff(t.INS) == 1;
        BeforeInsomnia  = [BeforeInsomnia; t(befInsIdxs,:)];
        aftInsIdxs      = find(befInsIdxs==1) + 1;
        if ~isempty(t(aftInsIdxs,:))
        AfterInsomnia   = [AfterInsomnia; t(aftInsIdxs,:)];
        end
    elseif sum(T.INS(ic == i)) == 0 && sum(ic == i) > 1
        t       = T(ic == i,:);
        Controlv1 = [Controlv1; t(1,:)];
        Controlv2 = [Controlv2; t(2,:)];
    end
end
BeforeInsomnia.INS      = [];
AfterInsomnia.INS       = [];
Controlv1.INS           = []; 
Controlv2.INS           = []; 

% Insomnia
for i = 2:size(BeforeInsomnia,2)
    X = AfterInsomnia{:,i} - BeforeInsomnia{:,i};
    if sum(ismember({'ANX','DEP','SEX'},AfterInsomnia.Properties.VariableNames{i})) == 1
        [~,p(i)] = chi2gof(X);
    else
        [~,p(i)] = ttest(X); 
    end
end
pValsFDRIns    = mafdr(p);
pValsBONIns    = p<(0.05/size(BeforeInsomnia,2));
idxsFDRIns     = find(pValsFDRIns<(0.05)) + 1;
idxsBONIns     = find(pValsBONIns) + 1;

% Control
for i = 2:size(Controlv1,2)
    X = Controlv2{:,i} - Controlv1{:,i};
    if sum(ismember({'ANX','DEP','SEX'},Controlv2.Properties.VariableNames{i})) == 1
        [~,p(i)] = chi2gof(X);
    else
        [~,p(i)] = ttest(X); 
    end
end
pValsFDRCon    = mafdr(p);
pValsBONCon    = p<(0.05/size(Controlv1,2));
idxsFDRCon     = find(pValsFDRCon<(0.05)) + 1;
idxsBONCon     = find(pValsBONCon) + 1;

idxsBONIns(ismember(idxsBONIns,idxsBONCon)==0)
idxsFDRIns(ismember(idxsFDRIns,idxsFDRCon)==0)

diffMatCon  = table2array(Controlv2(:,2:end)) - table2array(Controlv1(:,2:end));
diffMatIns  = table2array(AfterInsomnia(:,2:end)) - table2array(BeforeInsomnia(:,2:end));

diffTab     = array2table([diffMatCon;diffMatIns],...
    'VariableNames',Controlv2.Properties.VariableNames(2:end));

res     = table();
res.INS = [zeros(size(diffMatCon,1),1);ones(size(diffMatIns,1),1)];
for i = 1:size(diffTab,2)
    mdls            = fitglm([diffTab(:,i),res],'Distribution','binomial');
    pVal(i)         = mdls.Coefficients.pValue(2);
    effect(i)       = mdls.Coefficients.Estimate(2);
end

%%
BeforeInsomnia          = rmmissing(BeforeInsomnia(:,2:end));
Controlv1               = rmmissing(Controlv1(:,2:end));
CATexp  = [BeforeInsomnia.DEP,BeforeInsomnia.ANX,BeforeInsomnia.SEX]; BeforeInsomnia.DEP = []; BeforeInsomnia.ANX = []; BeforeInsomnia.SEX = [];
CATcon  = [Controlv1.DEP,Controlv1.ANX,Controlv1.SEX]; Controlv1.DEP = []; Controlv1.ANX = []; Controlv1.SEX = [];
Aa      = table2array(BeforeInsomnia);
Ca      = table2array(Controlv1);
Z       = zscore([Aa;Ca]);
Az      = Z(1:size(Aa,1),:);
Cz      = Z(1+size(Aa,1):end,:);
catDist = pdist2(CATexp,CATcon,'jaccard');
% Check for correlations
a = corrcoef([Az;Cz]);
b = abs(a)>0.95;
figure, subplot(3,1,1),imagesc(a), colorbar
subplot(3,1,2),imagesc(b), colorbar
c = [Az;Cz];
e = array2table(c,'VariableNames',BeforeInsomnia.Properties.VariableNames);
i = 1;
while i <= size(c,2)
    if sum(b(:,i)) > 1
        d = find(b(:,i) == 1); d(d==i) = [];
        c(:, d) = [];
        e(:, d) = [];
        BeforeInsomnia(:, d) = [];
        Controlv1(:, d) = [];
        a = corrcoef(c);
        b = abs(a)>0.95;
        i = 0;
    end
    i = i + 1;
end
subplot(3,1,3),imagesc(abs(b)>0.95), colorbar
[~,idxss]       = sort(std(c));
c(:,idxss(1:5)) = [];
BeforeInsomnia(:,idxss(1:5)) = [];
Controlv1(:,idxss(1:5)) = [];
Az = c(1:size(Az,1),:);
Cz = c(1+size(Az,1):end,:);
expirmentset    = Az;
%% All info
dists           = pdist2(Az(:,1:6),Cz(:,1:6))/6 + catDist/3;
[mins,idxss]    = sort(dists');
usedIdx         = [];
controlset      = nan(size(expirmentset));
for i = 1:size(idxss,2)
    if sum(ismember(usedIdx,idxss(1,i))) == 0
        controlset(i,:) = Cz(idxss(1,i),:);
        usedIdx = [usedIdx;idxss(1,i)];
    else
        found = false;
        ii = 2;
        while ~found
            if sum(ismember(usedIdx,idxss(ii,i))) == 0
            controlset(i,:) = Cz(idxss(ii,i),:);
            usedIdx         = [usedIdx;idxss(ii,i)];
            found           = true;
            end
            ii = ii + 1;
        end
    end
end
%% Remove noise and make data
criterion = ((std(expirmentset) < 0.1)*1 + (std(controlset) < 0.1))*1 > 0;

AaT      = BeforeInsomnia(:,7:end);          AaT(:,criterion) = [];
CaT      = Controlv1(usedIdx,7:end);    CaT(:,criterion) = [];
names = AaT.Properties.VariableNames;

controlset      = controlset(:,7:end);
controlset(:,criterion) = [];

expirmentset    = expirmentset(:,7:end);
expirmentset(:,criterion) = [];


AzT = array2table(expirmentset,'VariableNames',names);
CzT = array2table(controlset,'VariableNames',names);
labels   = [zeros(size(controlset,1),1);...
    ones(size(expirmentset,1),1)];
DATAz = [AzT;CzT];
% Classifier app
DATAz.Labels  = labels;
%% CSV write
csvwrite('Xfuture.csv',[controlset;expirmentset]);
csvwrite('Yfuture.csv',labels);
writetable(DATAz,'Zfuture.csv');
fid = fopen('featurenames.csv', 'w') ;
fprintf(fid, '%s,', names{1,1:end-1}) ;
fprintf(fid, '%s\n', names{1,end}) ;
fclose(fid) ;
dlmwrite('featurenames.csv', names(2:end,:), '-append') ;
%% No depression or anxiety
Aznm            = Az(sum(CATexp(:,1:2),2) == 0,:);
AaT             = AaT(sum(CATexp(:,1:2),2) == 0,:);
expirmentset    = Aznm;
Cznm            = Cz(sum(CATcon(:,1:2),2) == 0,:);
dists           = dists(sum(CATexp(:,1:2),2) == 0,sum(CATcon(:,1:2),2) == 0);
[mins,idxss]    = sort(dists');
usedIdx         = [];
controlset      = nan(size(expirmentset));
for i = 1:size(idxss,2)
    if sum(ismember(usedIdx,idxss(1,i))) == 0
        controlset(i,:) = Cz(idxss(1,i),:);
        usedIdx = [usedIdx;idxss(1,i)];
    else
        found = false;
        ii = 2;
        while ~found
            if sum(ismember(usedIdx,idxss(ii,i))) == 0
            controlset(i,:) = Cz(idxss(ii,i),:);
            usedIdx         = [usedIdx;idxss(ii,i)];
            found           = true;
            end
            ii = ii + 1;
        end
    end
end
csvwrite('Xnomental.csv',[controlset;expirmentset]);
labels   = [zeros(size(controlset,1),1);...
    ones(size(expirmentset,1),1)];
csvwrite('Ynomental.csv',labels);
% DATA  = [AaT;CaT(usedIdx,:)];
% DATAz = [AzT;CzT(usedIdx,:)];
% % Classifier app
% DATA.Labels  = labels;
% DATAz.Labels  = labels;
%% Other plots
% Anxiety, depression and insomnia a.f.o age
data    = [T.AGE, T.ANX, T.DEP, T.INS];
sex     = [T.SEX];
% Males
data(sex == 0,:) = [];
[~,idxs] = sort(data);
data = data(idxs(:,1),:);
data(any(isnan(data), 2), :) = [];
data = [data(:,1), movmean(data(:,2:end),150)];
figure,
hold on
plot(data(:,1),data(:,2),'-r')
plot(data(:,1),data(:,3),'-g')
plot(data(:,1),data(:,4),'-b')
grid on, box off
xlabel('Age'),ylabel('Probability')
% Females
data    = [T.AGE, T.ANX, T.DEP, T.INS];
data(sex == 1,:) = [];
[~,idxs] = sort(data);
data = data(idxs(:,1),:);
data(any(isnan(data), 2), :) = [];
data = [data(:,1), movmean(data(:,2:end),150)];
hold on
plot(data(:,1),data(:,2),'--r')
plot(data(:,1),data(:,3),'--g')
plot(data(:,1),data(:,4),'--b')
legend('Anxiety','Depression','Insomnia',...
    'Location','best')
title('--- (males), -\,-\,- (females)')
print(strcat(plotPath3,'MentalHealth_vs_Age'),'-depsc');



