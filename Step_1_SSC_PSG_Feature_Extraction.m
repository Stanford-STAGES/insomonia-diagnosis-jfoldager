% Author: Jonathan Foldager
%% Initial/startup code
startup      
% Set path and get list of all edf, sco, sta files in directory
extPSGDiskPath     = strcat('Data/ssc/psg/');
extGenDiskPath     = strcat('Data/ssc/genetics/');
% extDiskMetaPath = strcat('F:/DTU/Stanford/Data/ssc/_meta/');
listing         = extractfield(dir(extPSGDiskPath), 'name')';
edflisting      = listing(contains(lower(listing),'.edf'));
stalisting      = listing(contains(lower(listing),'.sta'));
listing         = extractfield(dir(extGenDiskPath), 'name')';
idConvTab       = strcat('Data/ssc/genetics/',listing(contains(lower(listing),'.csv')));
idConvTab       = readtable(idConvTab{1});
%% Load genetic data to see which subjects to include
try load('ssc_gen');catch, disp('No genetics file found'); end
%% Main loop
nSubjects       = length(edflisting);
nCurrentSubjects= nSubjects; 
FEATURES        = struct();
funny           = 0;
tic,for nSubject= 1:nCurrentSubjects
curSubj         = nSubject;
subjFile        = edflisting{curSubj}; 
subjID          = split(subjFile,'.'); subjID = strrep(subjID{1},' ','_'); id = split(subjID,'_'); id=str2num(id{2});
fullEDFPath     = strcat(extPSGDiskPath,subjFile);
fullSTAPath     = strcat(extPSGDiskPath,strrep(subjFile,'EDF','STA'));
%% Load PSG and event data
if exist(fullEDFPath, 'file') == 2  && ...
        exist(fullSTAPath, 'file') == 2  && ...  
                sum(idConvTab.DbID == id) == 1 && ...  
                    sum(strcmp(gen.id,idConvTab.GWASLocation(idConvTab.DbID == id))) == 1
            gwasID      = idConvTab.GWASLocation(idConvTab.DbID == id);
            gwasIDidx   = find(contains(gen.id,gwasID));
            try
            [header, data]  = loadEDF(fullEDFPath);
            events          = loadEVT(fullSTAPath,header.T0);
            for i = 1:length(header.label),fieldd=header.label{i};header.label{i}=fieldd(~isspace(fieldd));end
            catch
                funny = funny + 1;
                funnyIDs{funny} = id;
                continue;
            end
else
    funny = funny + 1;
    funnyIDs{funny} = id;
    continue;
end
fs              = header.samplerate(contains(header.label,'C3')); fs = fs(1);
if fs <= 70, continue; end
epochLength     = 30;
winLength       = 1;
%% Filtering
if sum(contains(header.label,'C3-A2')) == 1
    eegData = filtering(data{contains(header.label,'C3-A2')},fs,'eeg');
elseif sum(contains(header.label,'C3-M2')) == 1
    eegData = filtering(data{contains(header.label,'C3-M2')},fs,'eeg');
elseif sum(strcmp(header.label,'C3')) == 1
    eegData = filtering(data{strcmp(header.label,'C3')},fs,'eeg');
else 
    eegData = data(contains(header.label,'C3'));
    eegData = filtering(eegData{1},fs,'eeg');
end
    data{end+1} = eegData;
    header.label{end+1} = 'C3Ref';
    header.samplerate(end+1) = fs;
%% Hypnogram
data{end+1}                 = events.hypnogram;
header.label{end+1}         = 'Hypnogram';
header.samplerate(end+1)    = 1/epochLength;
%% Artefact removal
PSG                         = artefact_removal(data,header,events.lightsoffepoch,epochLength,'ssc');
%% Features
features                    = insomnia_hypnogram_analysis(events.hypnogram,epochLength);
features                    = concatstructs(features,insomnia_eeg_analysis(PSG,fs,winLength));
%%
usedIDs{curSubj} = gwasID{1};
FEATURES.(strcat('ssc_',gwasID{1})) = features;
save(strcat('Data/ssc/ssc_psg'),'FEATURES')
save(strcat('Data/SSC/ssc_psg'),'FEATURES')
fprintf('%i/%i done.',curSubj,nCurrentSubjects);
toc,end
