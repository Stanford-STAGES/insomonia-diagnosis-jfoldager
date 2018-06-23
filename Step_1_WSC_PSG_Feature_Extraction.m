% Author: Jonathan Foldager
%% Initial/startup code
startup      
% Set path and get list of all edf, sco, sta files in directory
extDiskPath     = strcat('Data/wsc/psg/');
extDiskMetaPath = strcat('Data/wsc/_meta/');
listing         = extractfield(dir(extDiskPath), 'name')';
% edflisting      = listing(contains(lower(listing),'.edf'));
edflisting      = listing(contains(listing,'.EDF'));
scolisting      = listing(contains(lower(listing),'.sco'));
stalisting      = listing(contains(lower(listing),'.sta'));
stalisting      = stalisting(~contains(lower(stalisting),'.sta2'));
slidePath       = 'figs\';
%% Load genetic data to see which subjects to include
try load('wsc_gen'); catch, fprintf('Genetics file does not exist!'); end
%% 
nSubjects       = length(wsc_gen.id);
nSubject        = 1;
nCurrentSubjects= nSubjects; 
FEATURES        = struct();
% tic,for nSubject= 1:nCurrentSubjects
curSubj         = nSubject;
subjFiles       = edflisting(contains(edflisting,wsc_gen.id{nSubject})); subjFile = subjFiles{1}; % Take first
subjID          = split(subjFile,'.'); subjID = strrep(subjID{1},' ','_'); id = split(subjID,'_'); id=id{1};
fullEDFPath     = strcat(extDiskPath,subjFile);
fullSCOPath     = strcat(extDiskPath,strrep(subjFile,'EDF','SCO'));
fullSTAPath     = strcat(extDiskPath,strrep(subjFile,'EDF','STA'));
%% Load PSG and event data
if exist(fullEDFPath, 'file') == 2 && ...           % edf must be present         exist(fullSCOPath, 'file') == 2 && ...      % sco must be present
            exist(fullSTAPath, 'file') == 2 && ...  % sta must be present
                sum(contains(wsc_gen.id,id)) == 1  && ...  % genetics must be present
                    sum(contains(IDs,edflisting{curSubj})) == 1 &&  ...  % Alex
                        sum(contains(usedIDs,id)) == 0
%     [header, data]  = loadEDF(fullEDFPath);
%     events          = loadEVT(fullSTAPath,header.T0);
    events.hypnogram    = class(ids.train == curSubj);
%     for i = 1:length(header.label),fieldd=header.label{i};header.label{i}=fieldd(~isspace(fieldd));end
else
%     continue;
end
winLength       = 1; epochLength = 30;
% fs              = header.samplerate(contains(header.label,'C3'));
%% Filtering
% str = {'F3','FZ','CZ','C3','PZ','O1','F4','C4'};
% idx = 1;
% for i = 1:length(data)
%     if contains(header.label{i},str)
%         data{i} = filtering(data{i},fs,'eeg');
%         eegChannels{idx} = header.label{i}; idx = idx + 1;
%         if contains(header.label{i},'C3'), header.label{i} = 'C3Ref'; end
%     end
% end
%% Hypnogram
% data{end+1}                 = events.hypnogram;
% header.label{end+1}         = 'Hypnogram';
% header.samplerate(end+1)    = 1/epochLength;
%% Artefact removal
% PSG                         = artefact_removal(data,header,events.lightsoffepoch,epochLength,'wsc');
%% Features
features                    = insomnia_hypnogram_analysis(events.hypnogram,epochLength);
% features.insomnia_eeg       = insomnia_eeg_analysis(PSG,fs,winLength);
% allFields                   = fields(events);
% for i = 1:length(allFields)
%     if contains(allFields{i},'h_')
%         features.(allFields{i}) = events.(allFields{i});
%     end
% end
%%
usedIDs{curSubj} = id;
FEATURES.(id) = features;
fprintf('%i/%i done.',curSubj,nCurrentSubjects);
% toc,end
% save(strcat('Data/wsc/'),'FEATURES_ALEX')