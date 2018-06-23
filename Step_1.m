startup
% Cohort folders MUST:
%       - be lower case
%       - not contain '.'
dataPath    = 'Data\test\';
dataDir     = dir(dataPath);
cohorts     = extractfield(dataDir, 'name')';
cohorts     = cohorts(cell2mat(extractfield(dataDir,'isdir')));
cohorts(contains(cohorts,'.')) = [];
coh         = 2;
sub         = 1;
N           = length(dir(strcat(dataPath,'**\*.edf')));
errorSubj   = zeros(size(cohorts));
shouldPlot  = false;
segLength   = 5;
nFeatures   = 88;
OUTPUT      = nan(N,nFeatures);
%% 
subOverall=1;
% for coh = 1:length(cohorts)   
curCohort   = cohorts{coh};
curDataPath = strcat(dataPath,curCohort,'\');
listing     = extractfield(dir(curDataPath), 'name')';
edflisting  = listing(endsWith(lower(listing),'.edf'));
nSubjPrCohort(coh)  = length(edflisting); 
if strcmp(curCohort,'ssc')
    stalisting      = listing(endsWith(lower(listing),'.sta'));
    evtlisting      = listing(endsWith(lower(listing),'.evts'));
    sscData         = load(strcat(dataPath,curCohort,'\SSC_demographics'));
    sscData         = sscData.t;
    ID              = string(num2cell(sscData.PatID));  
    AGE             = sscData.Age;
    BMI             = (sscData.Weight)./((sscData.Height).^2);
    SEX             = strcmpi(string(sscData.Gender),'male'); 
elseif strcmp(curCohort,'wsc')
    stalisting      = listing(endsWith(lower(listing),'.sta'));
    evtlisting      = listing(endsWith(lower(listing),'.sco'));
    csvlisting      = listing(endsWith(lower(listing),'.csv'));
    wscData         = load(strcat(dataPath,curCohort,'\WSC_demographics'));
    wscData         = wscData.t;
    ID              = strcat(wscData.PatID,'_',num2str(wscData.StudyNumber));   
    AGE             = wscData.Age;
    BMI             = wscData.BMI;
    SEX             = strcmpi(string(wscData.Gender),'male'); 
elseif strcmp(curCohort,'mros')
    evtlisting      = listing(endsWith(lower(listing),'.xml'));
    mrosTimes       = readtable(strcat(dataPath,curCohort,'\Times-visit1'));
    mrosData        = readtable(strcat(dataPath,curCohort,'\mros-visit1-dataset-0.3.0.csv'));
    ID              = mrosData.nsrrid;   
    AHI             = str2double(mrosData.pordi4pa);   
    PLMI            = mrosData.poavgplm;                
    BMI             = str2double(mrosData.hwbmi);  
    SEX             = mrosData.gender == 2;
    AGE             = mrosData.vsage1;
    load(strcat(dataPath,curCohort,'\mros_gen'));
end
%% Load PSG and event data
% tic,for sub = 1:nSubjPrCohort(coh)
curSubj         = sub; 
subjFile        = edflisting{curSubj}; 
subjID          = split(subjFile,'.'); 
if strcmp(curCohort,'ssc'),id = split(subjID{1},'_'); id = string(id{2});
elseif strcmp(curCohort,'wsc'),id = split(subjID{1},' '); id = id{1};
elseif strcmp(curCohort,'mros'),id = split(subjID{1},'-'); id = upper(id{3}); 
    idNum = str2double(strrep(id,'AA',''));
%     if sum(idNum == str2double(genetics.mros_gen.id)) ~= 1, 
%     errorSubj(coh)              = errorSubj(coh) + 1;
%     errorSubjs{errorSubj(coh)}  = subjID;
%     continue; end
end
fprintf('Running %s from %s\n',id,upper(curCohort))
idIdx           = find(strcmp(ID,id)); 
subjEnd         = subjID{2}; 
subjID          = upper(strrep(strrep(strrep(strrep(subjID{1},' ','_'),'A','WSC_'),'-','_'),'_visit1','_1')); 
fullEDFpath     = strcat(curDataPath,subjFile);
fullSTApath     = strcat(curDataPath,strrep(subjFile,subjEnd,'STA'));
fullSCOpath     = strcat(curDataPath,strrep(subjFile,subjEnd,'SCO'));
fullEVTpath     = strcat(curDataPath,strrep(subjFile,subjEnd,'EVTS'));
fullCSVpath     = strcat(curDataPath,strrep(subjFile,strcat('.',subjEnd),'_export.csv'));
fullXMLpath     = strcat(curDataPath,strrep(subjFile,strcat('.',subjEnd),'-nsrr.xml'));
% try
if exist(fullEDFpath, 'file') == 2  && ...
        exist(fullSTApath, 'file') == 2 && ...
            exist(fullSCOpath, 'file') == 2
            [header, data]  = loadEDF(fullEDFpath);
            events          = loadEVT(fullSTApath,fullSCOpath,header.T0);
            events          = concatstructs(events,CLASS_codec.parseSCOfile(fullSCOpath));
            [ahi,plmi,plmai]= analyze_wsc_events(events,header.samplerate(1));
            for i = 1:length(header.label),fieldd=header.label{i};header.label{i}=fieldd(~isspace(fieldd));end
elseif exist(fullEDFpath, 'file') == 2  && ...
        exist(fullCSVpath, 'file') == 2
            [header, data]  = loadEDF(fullEDFpath);
            [events,ahi,plmi,plmai]          = csv2events(fullCSVpath);
            for i = 1:length(header.label),fieldd=header.label{i};header.label{i}=fieldd(~isspace(fieldd));end
elseif exist(fullEDFpath, 'file') == 2  && ...
            exist(fullSTApath, 'file') == 2 && ...
                exist(fullEVTpath, 'file') == 2
            [header, data]  = loadEDF(fullEDFpath);
            events          = loadEVT(fullSTApath,fullEVTpath,header.T0);
            [SCOStruct, ~]  = CLASS_codec.parseSSCevtsFile(fullEVTpath);
            events          = concatstructs(events,SCOStruct);
            [ahi,plmi,plmai]= analyze_ssc_events(events);
            for i = 1:length(header.label),fieldd=header.label{i};header.label{i}=fieldd(~isspace(fieldd));end
elseif exist(fullEDFpath, 'file') == 2  && ...
            exist(fullXMLpath, 'file') == 2
            [header, data]  = loadEDF(fullEDFpath);
            events          = xml2events(fullXMLpath,mrosTimes.Start2LightsOff...
                (strcmp(mrosTimes.nsrrid,id)));
            ahi=AHI(idIdx);plmi=PLMI(idIdx);
            for i = 1:length(header.label),fieldd=header.label{i};header.label{i}=fieldd(~isspace(fieldd));end
else
    errorSubj(coh)              = errorSubj(coh) + 1;
    errorSubjs{errorSubj(coh)}  = subjID;
%     continue;
end
% catch
%     errorSubj(coh)              = errorSubj(coh) + 1;
%     errorSubjs{errorSubj(coh)}  = subjID;
% end
% if sum(contains(header.label,'C3')) == 0, continue; end
fs              = header.samplerate(contains(header.label,'C3')); fs = fs(1);
maxFs           = max(header.samplerate);
epochLength     = 30;
winLength       = 1;
%% Filtering
if sum(contains(header.label,'C3-A2')) == 1
    eegData = filtering(data{contains(header.label,'C3-A2')},fs,'eeg');
elseif sum(contains(header.label,'C3-M2')) == 1
    eegData = filtering(data{contains(header.label,'C3-M2')},fs,'eeg');
elseif sum(strcmp(header.label,'C3')) == 1 && sum(strcmp(header.label,'A2')) == 1 
    eegData = filtering(data{strcmp(header.label,'C3')} ...
        - data{strcmp(header.label,'A2')},fs,'eeg');
elseif sum(strcmp(header.label,'C3')) == 1
    eegData = filtering(data{strcmp(header.label,'C3')},fs,'eeg');
else 
    eegData = data(contains(header.label,'C3'));
    eegData = filtering(eegData{1},fs,'eeg');
end
data{end+1} = eegData;
header.label{end+1} = 'C3Ref';
header.samplerate(end+1) = fs;
eegDatas.C3Ref = eegData;
%% Hypnogram
data{end+1}                 = events.hypnogram;
header.label{end+1}         = 'Hypnogram';
header.samplerate(end+1)    = 1/epochLength;
%% Artefact removal and prepare ECG signal
PSG                         = artefact_removal(data,header,events,epochLength,lower(curCohort));
ecg                         = data(contains(header.label,{'ECG','EKG'}));
ecgFs                       = header.samplerate(contains(header.label,{'ECG','EKG'}));
if length(ecg) > 1, ecg = ecg(1); ecgFs = ecgFs(1); end
ecg                         = ecg{1};
ecgRes                      = reshape(ecg(1:end-mod(length(ecg) ...
                                ,ecgFs*epochLength)),ecgFs*epochLength,[]);
ecgRes                      = ecgRes(:,PSG.Org_LightsOffEpoch:PSG.Org_LastSleepEpoch);
ecg                         = ecgRes(:);    
%% Parameters
times                   = PSG.times;
hyp                     = PSG.Hypnogram;
timess.Pmodified        = floor(size(times,1).*PPeriods);
timess.Pmodified(:,1)   = timess.Pmodified(:,1) + 1;
for i = 1:size(CPeriods,1)
    logicall = (times(:,4:end) >= CPeriods(i,1:3)).*(times(:,4:end) <= CPeriods(i,4:6));
    logicall(any(logicall==0, 2), :)    = 0;
    indices                             = find(logicall(:,1));
    try
        timess.Cmodified(i,:)                      = [indices(1) indices(end)];
    catch
        timess.Cmodified(i,:)                      = [NaN NaN];
    end
end
test = eegData(:);
test = reshape(test(1:end-mod(length(test) ...
                                ,fs*30)),fs*30,[]);
win         = hamming(1*fs);
[pxx,f]     = pwelch(test,win,[],[],fs,'psd');

a = csvread('matlab_spectrum.csv');
b = csvread('python_spectrum.csv');
c = a - b;
figure
subplot(3,1,1)
imagesc(abs(c)), colorbar,axis xy, title('Absolute difference');
subplot(3,1,2)
imagesc(log(abs(c))), colorbar, axis xy, title('Log absolute difference');
subplot(3,1,3)
histogram(abs(sort(c(c<100))),1000), title('Histogram of absolute differences');

figure,plot(f,10*log10(pxx(:,1)))

figure
subplot(2,1,1)
imagesc([1 size(pxx,2)],[min(f) max(f)],10*log10(pxx)) , title('MATLAB pwelch output'), ylabel('Hz');
axis xy, colorbar
subplot(2,1,2)
imagesc([1 size(b,2)],[min(f) max(f)],b) , title('Scipy welch output'), ylabel('Hz');
axis xy, colorbar

caxis([0 50])
arousals = pxx(f>4,:)';
a = arousals(2:end,:)./arousals(1:end-1,:);
%% Feature extraction
SEQ_SPECTRAL    = eeg_welch(eegDatas,fs,timess,freqs,30,winLength);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% I.D.D %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
qs              = [0.25 0.50 0.75];
if ~isempty(idIdx)
% Demographics
IDD_DEMOG       = [coh, AGE(idIdx),SEX(idIdx),BMI(idIdx),ahi,plmi];
% Hypnogram
IDD_HYPNO       = struct2cell(hypnogram_analysis(events.hypnogram));
% EEG Spectral
IDD_SPECTRAL    = eeg_welch_summary(SEQ_SPECTRAL,qs);
OUTPUT(subOverall,:) = [IDD_DEMOG,IDD_HYPNO{1:22},IDD_SPECTRAL(1:2:end)];
subOverall = subOverall + 1;
end

