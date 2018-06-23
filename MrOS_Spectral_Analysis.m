startup
% Band summary
TBand = readtable('mros-visit2-eeg-band-summary-0.3.0.csv');
TBand = TBand(:,sort(TBand.Properties.VariableNames));
TBand = [TBand(:,contains(TBand.Properties.VariableNames,'nsrrid')) ...
    TBand(:,~contains(TBand.Properties.VariableNames,'nsrrid'))];
fBands = TBand(1,contains(TBand.Properties.VariableNames,'band'));
% Spectral summary
TSpectral = readtable('mros-visit2-eeg-spectral-summary-0.3.0.csv');
TSpectral = TSpectral(:,sort(TSpectral.Properties.VariableNames));
TSpectral = [TSpectral(:,contains(TSpectral.Properties.VariableNames,'nsrrid')) ...
    TSpectral(:,~contains(TSpectral.Properties.VariableNames,'nsrrid'))];
% Everything else
TData = readtable('mros-visit2-dataset-0.3.0.csv');
TData = TData(:,sort(TData.Properties.VariableNames));
TData = [TData(:,contains(TData.Properties.VariableNames,'nsrrid')) ...
    TData(:,~contains(TData.Properties.VariableNames,'nsrrid'))];
clc % cleans off warnings
% GEN = readtable('IMPUTED_EXAMPLE_JONATHAN.csv');
%% Convention
% '': Missing
% 1 : No clinically significant insomnia        0–7 
% 2 : Subthreshold insomnia                     8–14
% 3 : Clinical insomnia (moderate severity)     15–21
% 4 : Clinical insomnia (severe)                22–28
%% Extracting relevant info
% Get all that gave insomnia answer(s)
pop             = TData(str2double(TData.slisiscr) >= 0,:); 
% Get ID's
popIDs          = pop.nsrrid; 
% Grouping the insomniacs
ISI             = str2double(pop.slisiscr);
pop.ins         = ISI > 0;
pop.ISI         = ISI;
% Sort and save patients based on severity index
insomniaSorted = sortrows(pop,size(pop,2),'descend');
insomniaSorteds = insomniaSorted(:,contains(insomniaSorted.Properties.VariableNames,...
    {'nsrrid','inssev','slfalslp','slstyslp',...
    'slwkerly', 'slsatpat', 'slprnotc', 'slprworr', 'slprintr'}));%...
    %'poxusual'}));
insomniaSorteds.slfalslp = str2double (insomniaSorteds.slfalslp); pop.ins_slfalslp =  str2double (pop.slfalslp);
insomniaSorteds.slstyslp = str2double (insomniaSorteds.slstyslp); pop.ins_slstyslp = str2double (pop.slstyslp);
insomniaSorteds.slwkerly = str2double (insomniaSorteds.slwkerly); pop.ins_slwkerly = str2double (pop.slwkerly);
insomniaSorteds.slsatpat = str2double (insomniaSorteds.slsatpat); pop.ins_slsatpat = str2double (pop.slsatpat);
insomniaSorteds.slprnotc = str2double (insomniaSorteds.slprnotc); pop.ins_slprnotc = str2double (pop.slprnotc);
insomniaSorteds.slprworr = str2double (insomniaSorteds.slprworr); pop.ins_slprworr = str2double (pop.slprworr);
insomniaSorteds.slprintr = str2double (insomniaSorteds.slprintr); pop.ins_slprintr = str2double (pop.slprintr);
% insomniaSorteds.poxusual = str2double (insomniaSorteds.poxusual);
% save('Data/MrOS/insomniaseverity.mat','insomniaSorteds')
% Extract patients from "Spectrum Summary" and "Band Summary"
popBand         = TBand(ismember(TBand.nsrrid, popIDs),:);
popBand         = popBand(:,~contains(popBand.Properties.VariableNames,'band'));
popSpectrum     = TSpectral(ismember(TSpectral.nsrrid, popIDs),:);
popSpectrum     = popSpectrum(:,~contains(popSpectrum.Properties.VariableNames,'band'));
% Get insomnia ID's and their corresponding severity
insomniaIDs             = pop.nsrrid(pop.ins);
% Identify and assign insomniacs and their severity in other datasets
popBand.insomnia        = ismember(popBand.nsrrid, insomniaIDs);
for i=1:length(insomniaIDs),popBand.ISI(ismember(popBand.nsrrid, insomniaIDs{i}),1)=pop(contains(pop.nsrrid,insomniaIDs{i}),:).ISI; for ii = 1:7,popBand.ins_detailed(ismember(popBand.nsrrid, insomniaIDs{i}),ii) = pop{contains(pop.nsrrid,insomniaIDs{i}),end-7+ii};end, end
popSpectrum.insomnia    = ismember(popSpectrum.nsrrid, insomniaIDs);
for i=1:length(insomniaIDs),popSpectrum.ISI(ismember(popSpectrum.nsrrid, insomniaIDs{i}),1)=pop(contains(pop.nsrrid,insomniaIDs{i}),:).ISI; for ii = 1:7,popSpectrum.ins_detailed(ismember(popSpectrum.nsrrid, insomniaIDs{i}),ii) = pop{contains(pop.nsrrid,insomniaIDs{i}),end-7+ii};end, end
band.C3 = popBand(contains(popBand.signal,'C3'),:);
band.C4 = popBand(contains(popBand.signal,'C4'),:);
spec.C3 = popSpectrum(contains(popBand.signal,'C3'),:);
spec.C4 = popSpectrum(contains(popBand.signal,'C4'),:);
% % Transform tables to arrays
% band.C3Insomnias= table2array(band.C3(:,end-2)); band.C3InsomniaSeverity= table2array(band.C3(:,end));
% band.C4Insomnias= table2array(band.C4(:,end-1)); band.C4InsomniaSeverity= table2array(band.C4(:,end));
% spec.C3Insomnias= table2array(spec.C3(:,end-1)); spec.C3InsomniaSeverity= table2array(spec.C3(:,end));
% spec.C4Insomnias= table2array(spec.C4(:,end-1)); spec.C4InsomniaSeverity= table2array(spec.C4(:,end));
% Remove all un-usable datafields
band.C3(:,{'visit','signal'})      = [];
band.C4(:,{'visit','signal'})      = [];
spec.C3(:,{'visit','signal'})      = [];
spec.C4(:,{'visit','signal'})      = [];
% Transform tables to arrays
band.C3Mat      = table2array(band.C3(:,2:end-3)); band.C3Names = band.C3(:,2:end-3).Properties.VariableNames;
band.C4Mat      = table2array(band.C4(:,2:end-3)); band.C4Names = band.C4(:,2:end-3).Properties.VariableNames;
spec.C3Mat      = table2array(spec.C3(:,2:end-3)); spec.C3Names = spec.C3(:,2:end-3).Properties.VariableNames;
spec.C4Mat      = table2array(spec.C4(:,2:end-3)); spec.C4Names = spec.C4(:,2:end-3).Properties.VariableNames;
k               = 4;
band.C3Mat      = transpose(knnimpute(band.C3Mat', k,'Distance','euclidean'));
band.C4Mat      = transpose(knnimpute(band.C4Mat', k,'Distance','euclidean'));
spec.C3Mat      = transpose(knnimpute(spec.C3Mat', k,'Distance','euclidean'));
spec.C4Mat      = transpose(knnimpute(spec.C4Mat', k,'Distance','euclidean'));
% Split data into control and insomnia
% band.C3MatIns   = band.C3Mat(band.C3Insomnias,:);
% band.C3MatCon   = band.C3Mat(~band.C3Insomnias,:);
% band.C4MatIns   = band.C4Mat(band.C4Insomnias,:);
% band.C4MatCon   = band.C4Mat(~band.C4Insomnias,:);
% spec.C3MatIns   = spec.C3Mat(spec.C3Insomnias,:);
% spec.C3MatCon   = spec.C3Mat(~spec.C3Insomnias,:);
% spec.C4MatIns   = spec.C4Mat(spec.C4Insomnias,:);
% spec.C4MatCon   = spec.C4Mat(~spec.C4Insomnias,:);
% Replace missing values with average of k NN
% Use as many controls as insomniacs
% if 1 == 2
%     cBandPat = randperm(length(band.C3MatCon),length(band.C3MatIns));
%     cSpecPat = randperm(length(spec.C3MatCon),length(spec.C3MatIns));
%     band.C3MatCon = band.C3MatCon(cBandPat,:);
%     spec.C3MatCon = spec.C3MatCon(cSpecPat,:);
%     band.C3All    = [band.C3MatCon;band.C3MatIns];
%     spec.C3All    = [spec.C3MatCon;spec.C3MatIns];
% else
%     band.C3All    = [band.C3MatCon;band.C3MatIns];
%     spec.C3All    = [spec.C3MatCon;spec.C3MatIns];   
% end
%% Anomaly detection: Normal distribution test plot (QQ-plot)
% Extreme is obs. 433 and 95
band.C3(433,:)                  = [];
band.C3(95,:)                   = [];
band.C3Mat(433,:)               = [];
band.C3Mat(95,:)                = [];
spec.C3(433,:)                  = [];
spec.C3(95,:)                   = [];
spec.C3Mat(433,:)               = [];
spec.C3Mat(95,:)                = [];
% figure('units','normalized','outerposition',[0.1 0.05 0.6 0.9])
% for i = 1:size(band.C3Mat,2)
%     subplot(8,6,i)
%     qqplot(band.C3Mat(:,i));
%     xlabel('')
%     ylabel('')
%     title(band.C3Names(i))
%     axis tight
%     axis off
%     box off
% end
% print(strcat(plotPath2,'QQ_C3_Band.eps'),'-depsc')
% figure('units','normalized','outerposition',[0.1 0.05 0.6 0.9])
% for i = 1:size(band.C3Mat,2)
%     subplot(8,6,i)
%     boxplot(band.C3Mat(:,i));
%     xlabel('')
%     ylabel('')
%     title(band.C3Names(i))
%     axis tight
%     axis off
%     box off
% end
% print(strcat(plotPath2,'Boxplot_C3_Band.eps'),'-depsc')
%% PCA
[~,score,~] = pca(zscore(spec.C3Mat));
figure('units','normalized','outerposition',[0.15 0.05 0.7 0.9])
    subplot('Position',[0.4 0.4 0.25 0.55])
    gscatter(score(:,1),score(:,2),spec.C3.ISI);
    title('ISI')
    set(gca,'xtick',[])
    subplot('Position',[0.1 0.1 0.25 0.25])
    gscatter(score(:,1),score(:,2),spec.C3.ins_detailed(:,1));
    title('slfalslp')
    subplot('Position',[0.4 0.1 0.25 0.25])
    gscatter(score(:,1),score(:,2),spec.C3.ins_detailed(:,2));
    title('slstyslp')
    xlabel('Principal Component 1')
    subplot('Position',[0.7 0.1 0.25 0.25])
    gscatter(score(:,1),score(:,2),spec.C3.ins_detailed(:,3));
    title('slwkerly')
    subplot('Position',[0.1 0.4 0.25 0.25])
    gscatter(score(:,1),score(:,2),spec.C3.ins_detailed(:,4));
    title('slsatpat')
    set(gca,'xtick',[])
    ylabel('Principal Component 2')
    subplot('Position',[0.7 0.4 0.25 0.25])
    gscatter(score(:,1),score(:,2),spec.C3.ins_detailed(:,5));
    title('slprnotc')
    set(gca,'xtick',[])
    subplot('Position',[0.1 0.7 0.25 0.25])
    gscatter(score(:,1),score(:,2),spec.C3.ins_detailed(:,6));
    title('slprworr')
    set(gca,'xtick',[])
    subplot('Position',[0.7 0.7 0.25 0.25])
    gscatter(score(:,1),score(:,2),spec.C3.ins_detailed(:,7));
    title('slprintr')
    set(gca,'xtick',[])
% print(strcat(plotPath2,'Scatter_PCA_C3_AllQuestions.eps'),'-depsc')    
%% Outliers
fig = figure(2);
for i = 1:size(band.C3Mat,2)
    subplot(8,6,i)
    qqplot(band.C3Mat(:,i));
    xlabel('')
    ylabel('')
    title(band.C3.Properties.VariableNames(i))
    axis tight
    axis off
    box off
end

figure(3);
data = [str2double(pop.hwbmi), ...
    str2double(pop.girace), ...
    str2double(pop.vs2age1), ...
    str2double(pop.epcar), ...
    str2double(pop.epeat), ...
    str2double(pop.epeds), ...
    str2double(pop.eppub), ...
    str2double(pop.epread), ...
    str2double(pop.eprest), ...
    str2double(pop.eptalk), ...
    str2double(pop.eptraf), ...
    str2double(pop.eptv), ...
    str2double(pop.slfalslp), ...
    str2double(pop.sllivlev), ...
    str2double(pop.slprintr), ...
    str2double(pop.slprnotc), ...
    str2double(pop.slprworr), ...
    str2double(pop.poxbeer), ...
    str2double(pop.poxcig), ...
    str2double(pop.poxcigar), ...
    str2double(pop.poxcoff), ...
    str2double(pop.poxqual1), ...
    str2double(pop.poxqual2), ...
    str2double(pop.poxqual3), ...
    str2double(pop.poxslpmn)];
data = transpose(knnimpute(data', k,'Distance','euclidean'));
data = zscore(data);
[coeff,score,~] = pca(data);
gscatter(score(:,1),score(:,2),pop.ins + pop.ISI > 2)
legend({'Control','Insomnia'})
xlabel('Principal Component 1')
ylabel('Principal Component 2')













