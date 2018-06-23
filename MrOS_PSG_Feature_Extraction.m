% Author: Jonathan Foldager
%% Initial/startup code
startup      
visit           = 1;
refElectrode    = {'A','M'}; refElectrode = {strcat(refElectrode{visit},'1'), strcat(refElectrode{visit},'2')};
extDiskPSGPath  = strcat('edfs/visit',num2str(visit),'/');
extDiskHypPath  = strcat('mros/polysomnography/annotations-events-nsrr/visit',num2str(visit),'/');
extDiskMetaPath = strcat('mros/polysomnography/events/visit',num2str(visit));
listing         = extractfield(dir(extDiskPSGPath), 'name')';
listing         = listing(contains(listing,'.edf'));
nSubjects       = length(listing);
nSubject        = 1;
T               = readtable(strcat('Times-visit',num2str(visit))); 
nCurrentSubjects= nSubjects; 
FEATURES        = struct();
douches         = {};
tic,for nSubject  = 1:nCurrentSubjects
try
curSubj         = nSubject;
if visit == 2 && true
    load('insomniaseverity.mat'); 
    worstInsomniacs = insomniaSorteds.nsrrid(1:nCurrentSubjects); 
    leastInsomniacs = insomniaSorteds.nsrrid(end-nCurrentSubjects-1:end);
    subjID = strcat('mros-visit',num2str(visit),'-', worstInsomniacs{curSubj});
elseif visit == 2 && false
    load('insomniaseverity.mat'); 
    subjID = split(listing(curSubj),'.'); subjID = subjID{1};
elseif visit == 1
    subjID = split(listing(curSubj),'.'); subjID = subjID{1};
end
% index = find(strcmp(lower(insomniaSorteds.nsrrid), nsrrid));
nsrrid = split(subjID,'-'); nsrrid = nsrrid{end};
fullPathAndFile = strcat(extDiskPSGPath,subjID,'.edf');
if exist(fullPathAndFile, 'file') == 2
    [header, data]  = loadEDF(fullPathAndFile);
    for i = 1:length(header.labels),fieldd=header.labels{i};header.labels{i}=fieldd(~isspace(fieldd));end
end
lightsOffSec    = T.Start2LightsOff(contains(lower(T.nsrrid), lower(nsrrid)));
winLength       = 1; epochLength = 30;
fs              = header.samplerate(contains(header.labels,'C3'));
%% Filtering
C3Ref    = filtering(data{contains(header.labels,'C3')} ...
    - data{contains(header.labels,refElectrode(2))},fs,'eeg');
C3       = filtering(data{contains(header.labels,'C3')},fs,'eeg');
C4Ref       = filtering(data{contains(header.labels,'C4')} ...
    - data{contains(header.labels,refElectrode(1))},fs,'eeg');
C4       = filtering(data{contains(header.labels,'C4')},fs,'eeg');

data{contains(header.labels,'C3')} = C3;
data{contains(header.labels,'C4')} = C4;
data{end+1} = C3Ref; header.labels{end+1} = 'C3Ref'; header.samplerate(end+1) = fs;
data{end+1} = C4Ref; header.labels{end+1} = 'C4Ref'; header.samplerate(end+1) = fs;
%% Hypnogram
if exist(strcat(extDiskMetaPath,'_',nsrrid,'.mat'), 'file') == 2
    load(strcat(extDiskMetaPath,'_',nsrrid));
else
    events = xml2events(strcat(extDiskHypPath,...
        'mros-visit',num2str(visit),'-',nsrrid,'-nsrr.xml'),epochLength,lightsOffSec);
end
data{end+1}                 = events.hyp;
header.labels{end+1}        = 'Hypnogram';
header.samplerate(end+1)    = 1/30;
features.insomnia_hypnogram = insomnia_hypnogram_analysis(events.sleepHyp,epochLength);
% features.hypnogram          = hypnogram_analysis(events.sleepHyp);
%% Artefact removal
PSG                         = artefact_removal(data,header,events.lightsoffSec,epochLength);
% interval                    = 1:PSG.New_LastSleepEpoch;
%% ECG
% ecg   = data{contains(header.labels,'ECGR')} - data{contains(header.labels,'ECGL')};
% ecgFs   = header.samplerate(contains(header.labels,'ECGL'));
% ecg = reshape(ecg(1:end-mod(length(ecg) ...
%             ,ecgFs*epochLength)),ecgFs*epochLength,[]); 
% ecg = ecg(:,PSG.Org_LightsOffEpoch:end); ecg = ecg(:);
%% Features
% features.bandpower  = power_spectral_analysis(PSG.C3Ref(:,interval),fs,winLength,epochLength);
% features.coherence  = coherence(PSG.C4Ref(:,interval),PSG.C3Ref(:,interval),fs,winLength,epochLength);
% features.hpC3       = hjort_parameters(PSG.C3Ref(:,interval),epochLength,fs,250);
% features.hpC4       = hjort_parameters(PSG.C4Ref(:,interval),epochLength,fs,250);
% features.br         = brain_rate(PSG.C3Ref(:,interval),fs,winLength,epochLength);
% features.bsi        = brain_symmetry_index(PSG.C4Ref(:,interval),PSG.C3Ref(:,interval),fs,winLength,epochLength );
% features.statistics = temporal_statistics(PSG.C3Ref(:,interval),fs,epochLength);
% features.hrv        = heart_analysis(ecg,ecgFs,PSG.Hypnogram,true,PSG.ArtefactEpochs,epochLength);
features.insomnia_eeg= insomnia_eeg_analysis(PSG,fs,winLength);
FEATURES.(nsrrid) = features;
fprintf('%i/%i done.',curSubj,nCurrentSubjects);
catch
    fprintf('%s is funny.',nsrrid);
    douches{nSubject} = nsrrid;
end
toc,end
save(strcat('mros/features/'),'FEATURES')