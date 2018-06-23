% Author: Jonathan Foldager
%% Initial/startup code
startup      
% Set path and get list of all edf, sco, sta files in directory
extDiskPath     = strcat('Data/mros/psg/edfs/visit1/');
extEvtDiskPath  = strcat('Data/mros/psg/annotations-events-nsrr/visit1/');
times           = readtable('Times-visit1.csv');
listing         = extractfield(dir(extDiskPath), 'name')';
edflisting      = listing(contains(lower(listing),'.edf'));
slidePath       = '\figs\';
visit           = 1;
refElectrode    = {'A','M'}; refElectrode = {strcat(refElectrode{visit},'1'), strcat(refElectrode{visit},'2')};
%% Load genetic data to see which subjects to include
% takeAll = false;
% try load('mros_gen'); catch, 
    takeAll = true; 
%     fprintf('Genetics file does not exist!'); end
%% 
nSubjects       = length(edflisting);
nSubject        = 1;
nCurrentSubjects= nSubjects; 
FEATURES        = struct();
funny           = 1;
tic,for nSubject= 2746:nCurrentSubjects
curSubj         = nSubject;
if takeAll, subjFile = edflisting{nSubject};  else, subjFiles = edflisting(contains(edflisting,wsc_gen.id{nSubject})); subjFile = subjFiles{1}; end % Take first
subjID          = split(subjFile,'.'); subjID = strrep(subjID{1},' ','_'); id = split(subjID,'-'); id=id{end};
fullEDFPath     = strcat(extDiskPath,subjFile);
fullEvtPath     = strcat(extEvtDiskPath,strrep(subjFile,'.edf','-nsrr.xml'));
winLength       = 1; epochLength = 30;
%% Load PSG and event data
if exist(fullEDFPath, 'file') == 2  && ...
        exist(fullEvtPath, 'file') == 2  
            try
            [header, data]  = loadEDF(fullEDFPath);
            events          = xml2events(fullEvtPath,epochLength,times.Start2LightsOff(contains(lower(times.nsrrid),id)));
            for i = 1:length(header.label),fieldd=header.label{i};header.label{i}=fieldd(~isspace(fieldd));end
            catch
                funnyIDs{funny} = id;
                funny = funny + 1;
                continue;
            end
else
    continue;
end
fs              = header.samplerate(contains(header.label,'C3'));
%% Filtering
C3Ref = filtering(data{contains(header.label,'C3')} - data{contains(header.label,refElectrode{2})},fs,'eeg');
C4Ref = filtering(data{contains(header.label,'C4')} - data{contains(header.label,refElectrode{1})},fs,'eeg');
data{contains(header.label,'C3')} = filtering(data{contains(header.label,'C3')},fs,'eeg');
data{contains(header.label,'C4')} = filtering(data{contains(header.label,'C4')},fs,'eeg');
data{end+1}                 = data{contains(header.label,'C3')}; % C3Ref;
data{end+1}                 = data{contains(header.label,'C4')}; % C4Ref;
header.label{end+1}         = 'C3Ref'; 
header.label{end+1}         = 'C4Ref'; 
header.samplerate(end+1)    = fs; 
header.samplerate(end+1)    = fs; 
%% Hypnogram
data{end+1}                 = events.hypnogram;
header.label{end+1}         = 'Hypnogram';
header.samplerate(end+1)    = 1/epochLength;
%% Artefact removal
PSG                         = artefact_removal(data,header,events.lightsoffEpoch,epochLength,'mros');
%% Features
features                    = insomnia_hypnogram_analysis(events.hypnogram,epochLength);
features                    = concatstructs(features,insomnia_eeg_analysis(PSG,fs,winLength));
%%
usedIDs{curSubj} = id;
FEATURES.(id) = features;
save(strcat('Data/mros/mros_psg'),'FEATURES')
fprintf('%i/%i done.',curSubj,nCurrentSubjects);
toc,end