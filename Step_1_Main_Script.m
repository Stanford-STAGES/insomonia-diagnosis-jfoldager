startup
%% Load genetic data
load('wsc_gen_new')
% wsc     = load('wsc_psg');
load('mros_gen');
% mro     = load('mros_psg');
load('ssc_gen')
% ssc 	= load('ssc_psg');
rsTable = readtable('MrOS_rs429358.xlsx');
bIds = readtable('MrOS_IDRosettaStone_NSRR-HAdbGaP-MrOSID.xlsx');
bIds = rmmissing(bIds,1);
bIds.rs42 = nan(size(bIds,1),1);
for i = 1:size(bIds,1)
    if sum(contains(rsTable.ID,bIds.MrOS_ID{i})) == 1
        dose = rsTable.Dose_of_C(contains(rsTable.ID,bIds.MrOS_ID{i}));
        if dose == 0
        bIds.rs42(i) = 2;
        elseif dose == 2
            bIds.rs42(i) = 0;
        else 
            bIds.rs42(i) = 1;
        end
        
    end
end
bIds = rmmissing(bIds,1);
%% Load demographics
MrOS = readtable('MrOSV2.csv');
WSC = readtable('WSCv1.csv'); WSC.ANX = []; WSC.INS = []; WSC.DEP = [];
SSC = readtable('SSC.csv'); SSC.ANX = []; SSC.INS = []; SSC.DEP = [];
%% Get common SNPs
sscSNPs     = gen.snps(:,2);
wscSNPs     = wsc_gen.snps(:,2);
mroSNPs     = mros_gen.snps(:,2);
allSNPs     = mintersect(sscSNPs,wscSNPs, mroSNPs);
sscSNPsIdxs = ismember(sscSNPs,allSNPs);
sscSNPs     = sscSNPs(sscSNPsIdxs);
wscSNPsIdxs = ismember(wscSNPs,allSNPs);
wscSNPs     = wscSNPs(wscSNPsIdxs);
mroSNPsIdxs = ismember(mroSNPs,allSNPs);
mroSNPs     = mroSNPs(mroSNPsIdxs);
allChrs     = wsc_gen.snps(wscSNPsIdxs);
allSNPs     = mroSNPs; % get the right order
%% Minor allele frequency
w               = repmat(wscSNPsIdxs',3,1); w = w(:);
s               = repmat(sscSNPsIdxs',3,1); s = s(:);
m               = repmat(mroSNPsIdxs',3,1); m = m(:);
XX              = [wsc_gen.data(:,w);gen.data(:,s);mros_gen.data(:,m)];
XXX             = sum(XX);
XXX             = XXX.*2;
XXXX            = reshape(XXX,3,[]);
allDosages      = sum(XXXX);
minAllelFreq    = ((XXXX(2,:)./2)+XXXX(3,:))./allDosages;
minAllelFreqIdx = find(minAllelFreq < 0.05);
minAllelFreq    = minAllelFreq(minAllelFreqIdx);
minAllelFreqSNP = allSNPs(minAllelFreqIdx);
%% WSC eeg correction
allSubjects = wsc_gen.id;
%% WSC
Xwsc = wsc_gen.data_expected(:,wscSNPsIdxs);
Y = [];
WSCsubjs = [];
for i = 1:length(allSubjects)
    if sum(contains(WSC.ID,allSubjects{i})) == 1
       Y = [Y; WSC(contains(WSC.ID,allSubjects{i}),:)];
    end
end
Ynonan = rmmissing(Y,1);
Xwsc(~ismember(allSubjects,Y.ID),:) = [];
cohortID = array2table([repmat("WSC",length(allSubjects(ismember(allSubjects,Y.ID))),1)...
    ,allSubjects(ismember(allSubjects,Y.ID))]...
    ,'VariableNames',{'Cohort','ID'});
% Write tables
XtableWSC  = array2table(Xwsc,'VariableNames',allSNPs);
YtableWSC  = Y; 
YtableWSC.CohortID = ones(size(YtableWSC,1),1);
Xall       = XtableWSC;
Yall       = YtableWSC(:,2:end);
XallIDs    = [YtableWSC.ID, XtableWSC]; XallIDs.Properties.VariableNames{1} = 'ID';
YallIDs    = YtableWSC;

% [rhoWSC,pvalWSC] = corrcoef(table2array(Ynonan(:,2:end)));
% figure,imagesc(rhoWSC)
% ax = gca; 
% ax.YTick = 1:length(rhoWSC); ax.YTickLabel = Yall.Properties.VariableNames;
% ax.XTick = 1:length(rhoWSC); ax.XTickLabel = Yall.Properties.VariableNames;
% xtickangle(ax,90), colorbar, axis xy
% print('Plots/Poster/CorrWSC','-depsc')
%% MrOS
ID      = readtable('Data/MrOS_IDRosettaStone_NSRR-HAdbGaP-MrOSID.xlsx');
idx     = 1;
idxx    = 1;
genIDs  = cellfun(@str2double,mros_gen.id);
mros_gen.data_expected = mros_gen.data_expected(:,mroSNPsIdxs);
T = table();
for i = 1:size(bIds,1)
    if sum(ismember(genIDs,bIds.HA_ID(i))) == 1
        T = [T; array2table(upper(bIds.NSRR_ID(i)),'VariableNames',{'ID'}),...
            array2table(mros_gen.data_expected(ismember(genIDs,bIds.HA_ID(i)),:),...
            'VariableNames',mroSNPs)];
    end
end
T.rs429358 = bIds.rs42;
names = MrOS.Properties.VariableNames;
for i = 2:length(names), T.(names{i}) = nan(size(T,1),1);end
for ii = 1:size(T,1)
    if sum(contains(MrOS.ID,T.ID{ii})) == 1
for i = 2:length(names), T.(names{i})(ii) = MrOS.(names{i})(contains(MrOS.ID,T.ID{ii})); end
    end
end
T           = rmmissing(T,1);
XtableMRO   = T(:,2:240);
YtableMRO   = T(:,241:end);
YtableMRO   = YtableMRO(:,sort(YtableMRO.Properties.VariableNames));
YtableMRO.CohortID = ones(size(YtableMRO,1),1) + 1;
Xall        = [Xall;XtableMRO];
Yall        = [Yall;YtableMRO];
XallIDs     = [XallIDs;T(:,1),XtableMRO];
YallIDs     = [YallIDs;T(:,1),YtableMRO];

% Ynonan      = YtableMRO; Ynonan.SEX = [];  Ynonan.CohortID = [];
% [rhoMrOS,pvalMrOS] = corrcoef(table2array(Ynonan));
% figure,imagesc(rhoMrOS)
% ax = gca; 
% names = Yall.Properties.VariableNames;
% names(contains(names,'SEX')) = [];
% ax.YTick = 1:length(rhoMrOS); ax.YTickLabel = names;
% ax.XTick = 1:length(rhoMrOS); ax.XTickLabel = names;
% xtickangle(ax,90), colorbar, axis xy
% print('Plots/Poster/CorrMrOS','-depsc')
%% SSC
T = table();
Xssc = [];
Yssc = table();
a = 0;
for i = 1:size(SSC,1)
    if sum(contains(gen.id,SSC.GWAS{i})) == 1
        Xssc = [Xssc;gen.data_expected(contains(gen.id,SSC.GWAS{i}),sscSNPsIdxs)];
        Yssc = [Yssc;SSC(i,:)];
    end
end
Yssc.CohortID = ones(size(Yssc,1),1) + 2;
Yssc.GWAS = [];
Yssc.ID = string(Yssc.ID);
Yssc = [Yssc(:,1),Yssc(:,sort(Yssc.Properties.VariableNames(2:end)))];
YallIDs     = [YallIDs;Yssc];
Xssc        = array2table(Xssc,'VariableNames',sscSNPs);
Xall        = [Xall;Xssc];
XallIDs     = [XallIDs;Yssc(:,1),Xssc];

Yssc.ID = [];
Yall        = [Yall;Yssc];

X           = table2array(Xall);
Y           = table2array(Yall);

YallIDs.Cohort(YallIDs.CohortID == 1)   = "WSC";
YallIDs.Cohort(YallIDs.CohortID == 2)   = "MrOS";
YallIDs.Cohort(YallIDs.CohortID == 3)   = "SSC";
XallIDs.Cohort(YallIDs.CohortID == 1)   = "WSC";
XallIDs.Cohort(YallIDs.CohortID == 2)   = "MrOS";
XallIDs.Cohort(YallIDs.CohortID == 3)   = "SSC";
YallIDs.CohortID                        = [];

YallIDs = [YallIDs(:,end), YallIDs(:,1:end-1)];
XallIDs = [XallIDs(:,end), XallIDs(:,1:end-1)];

writetable(YallIDs,'PhenotypesV1.csv');
writetable(XallIDs,'GenotypesV1.csv');

% Ynonan      = rmmissing(Yssc,1); Ynonan.CohortID = [];
% [rhoSSC,pvalSSC] = corrcoef(table2array(Ynonan));
% figure,imagesc(rhoSSC)
% ax = gca; 
% names = Yall.Properties.VariableNames;
% names(contains(names,'CohortID')) = [];
% ax.YTick = 1:length(rhoSSC); ax.YTickLabel = names;
% ax.XTick = 1:length(rhoSSC); ax.XTickLabel = names;
% xtickangle(ax,90), colorbar, axis xy
% print('Plots/Poster/CorrSSC','-depsc')

% Ynonan      = rmmissing(Yall,1);%  Ynonan.CohortID = [];
% [rhoY,pvalY] = corrcoef(table2array(Ynonan));
% figure,imagesc(rhoY)
% ax = gca; 
% names = Yall.Properties.VariableNames;
% ax.YTick = 1:length(rhoY); ax.YTickLabel = names;
% ax.XTick = 1:length(rhoY); ax.XTickLabel = names;
% xtickangle(ax,90), colorbar, axis xy
% print('Plots/Poster/CorrAll','-depsc')

% Ynonan      = rmmissing(Yall,1);%  Ynonan.CohortID = [];
% [rhoY,pvalY] = corrcoef(table2array(Ynonan));
% figure,imagesc(abs(rhoY)>0.9)
% ax = gca; 
% names = Yall.Properties.VariableNames;
% ax.YTick = 1:length(rhoY); ax.YTickLabel = names;
% ax.XTick = 1:length(rhoY); ax.XTickLabel = names;
% xtickangle(ax,90), c = colorbar; c.Ticks = [0,1]; colormap(gray(2)), axis xy
% print('Plots/Poster/CorrIsAbove90','-depsc')




