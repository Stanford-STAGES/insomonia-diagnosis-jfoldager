startup
% This script calculates temporal and spectral features from an .EDF
% in specified segLength and stores the array along with demographics
% in a .h5 file.
% Input: .EDF, .csv containing demographics and .EVT file
% Ouput: .h5 and .csv file containing features summaries and a.f.o. time
%%
% Cohort folders MUST:
%       - be lower case
%       - not contain '.'
dataPath    = 'Data\';
dataDir     = dir(dataPath);
cohorts     = extractfield(dataDir, 'name')';
cohorts     = cohorts(cell2mat(extractfield(dataDir,'isdir')));
cohorts((contains(cohorts,'_')+contains(cohorts,'.'))>0) = [];
coh         = 1;
sub         = 1;
N           = length(dir(strcat(dataPath,'**\*.edf')));
ids         = cell(N,2);                
errorSubj   = zeros(size(cohorts));
shouldPlot  = false;
segLength   = 30;
features    = table();
%% Main loop
pth = 'wsc\psgs\';
for coh = 1:length(cohorts)   
curCohort           = cohorts{coh};
curDataPath         = strcat(dataPath,curCohort,'\');
listing             = extractfield(dir(curDataPath), 'name')';
listing             = extractfield(dir(pth), 'name')';
edflisting          = listing(endsWith(lower(listing),'.edf'));
nSubjPrCohort(coh)  = length(edflisting); 
if strcmp(curCohort,'ssc')
    stalisting      = listing(endsWith(lower(listing),'.sta'));
    sscData         = readtable('ssc\Events.csv');
    tic,for sub = 1:nSubjPrCohort(coh) 
    try
    subjEDFFile     = edflisting{sub}; 
    subjSTAFile     = strrep(subjEDFFile,'.EDF','.STA'); 
    subjID          = split(subjEDFFile,'_'); 
    id              = subjID{2};
    [header, data]  = loadEDF(strcat(pth,subjEDFFile));
    events          = loadEVT(strcat(pth,subjSTAFile));
    for i = 1:length(header.label),fieldd=header.label{i};header.label{i}=fieldd(~isspace(fieldd));end
    EEGfs           = header.samplerate(contains(header.label,'C3')); EEGfs = EEGfs(1);
    epochLength     = 30;
    winLength       = 1;
%     EEG
    if sum(contains(header.label,'C3-A2')) == 1
        EEG = filtering(data{contains(header.label,'C3-A2')},EEGfs,'eeg');
    elseif sum(contains(header.label,'C3-M2')) == 1
        EEG = filtering(data{contains(header.label,'C3-M2')},EEGfs,'eeg');
    elseif sum(strcmp(header.label,'C3')) == 1 && sum(strcmp(header.label,'A2')) == 1 
        EEG = filtering(data{strcmp(header.label,'C3')} ...
            - data{strcmp(header.label,'A2')},EEGfs,'eeg');
    elseif sum(strcmp(header.label,'C3-')) == 1
        EEG = filtering(data{strcmp(header.label,'C3-')},EEGfs,'eeg');
    else 
        EEG = data(contains(header.label,'C3'));
        EEG = filtering(EEG{1},EEGfs,'eeg');
    end
%     ECG
    ECG                         = data(contains(header.label,{'ECG','EKG'}));
    ECGfs                       = header.samplerate(contains(header.label,{'ECG','EKG'}));
    if length(ECG) > 1, ECG = ECG(1); ECGfs = ECGfs(1); end
    ECG                         = ECG{1};
    ecgRes                      = reshape(ECG(1:end-mod(length(ECG) ...
                                    ,ECGfs*epochLength)),ECGfs*epochLength,[]);
    ecgRes                      = ecgRes(:,events.LOF:events.LON);
    ECG                         = ecgRes(:);  
    Features
    EEGTEMPORAL                 = temporal_statistics(EEG,EEGfs,segLength,events);
    HYPNOGRAM                   = hypnogram_analysis(events.hypnogram);
    EEGSPECTRAL                 = spectral_statistics(EEG,EEGfs,segLength,events);
    EEGARMODEL                  = estimate_ARmodel(EEG,EEGfs,segLength,5);
    ECGHRV                      = heart_analysis(ECG,ECGfs,segLength,events);
    sscData.ARI(strcmp(string(sscData.ID),id)) = sscData.ARI(strcmp(string(sscData.ID),id))/(HYPNOGRAM.TST/60);
    sscData.PLMI(strcmp(string(sscData.ID),id)) = sscData.PLMI(strcmp(string(sscData.ID),id))/(HYPNOGRAM.TST/60);
    sscData.AHI(strcmp(string(sscData.ID),id)) = sscData.AHI(strcmp(string(sscData.ID),id))/(HYPNOGRAM.TST/60);
%     Save
    features = [features; sscData(strcmp(string(sscData.ID),id),:), ...
        HYPNOGRAM, EEGTEMPORAL, EEGSPECTRAL, ECGHRV];
    writetable(features,'_features\SSC\SSC.csv');  
    writetable(EEGARMODEL(events.LOF:events.LON,:),strcat('_features\SSC\AR_',id,'.csv'));  
    try
    h5create(strcat('_features\SSC\AR_',id,'.h5'),'/values',size(table2array(EEGARMODEL)))
    h5write(strcat('_features\SSC\AR_',id,'.h5'),'/values',table2array(EEGARMODEL));
    h5create(strcat('_features\SSC\AR_',id,'_chunk.h5'),'/values',size(table2array(EEGARMODEL)),'ChunkSize',[6 5])
    h5write(strcat('_features\SSC\AR_',id,'_chunk.h5'),'/values',table2array(EEGARMODEL));
    catch
    end
    fprintf('%i/%i done.',sub,nSubjPrCohort(coh));
    catch
        errorSubj(coh)              = errorSubj(coh) + 1;
        errorSubjs{errorSubj(coh)}  = id;
        fprintf('%s is funny.',id);
    end
    toc,end
elseif strcmp(curCohort,'wsc')
    wscData         = readtable('WSCall.csv');           
    wscData.TST     = [];          
    wscData.WA     = [];          
    wscData.N1     = [];          
    wscData.N2     = [];          
    wscData.N3     = [];          
    wscData.REM     = [];         
    wscData.REML     = [];        
    wscData.SOL     = [];        
    wscData.SE     = [];        
    wscData.AW25     = [];        
    wscData.AW5     = [];        
    wscData.AWN125     = [];        
    wscData.AWN15     = [];
    stalisting      = listing(endsWith(lower(listing),'.sta'));
    tic,for sub = 1:nSubjPrCohort(coh) 
    try
    subjEDFFile     = edflisting{sub}; 
    subjSTAFile     = strrep(subjEDFFile,'.EDF','.STA'); 
    subjID          = split(subjEDFFile,' '); 
    id              = subjID{1};
    [header, data]  = loadEDF(strcat('wsc\psgs\',subjEDFFile));
    events          = loadEVT(strcat('wsc\psgs\',subjSTAFile));
    for i = 1:length(header.label),fieldd=header.label{i};header.label{i}=fieldd(~isspace(fieldd));end
    EEGfs           = header.samplerate(contains(header.label,'C3')); EEGfs = EEGfs(1);
    epochLength     = 30;
    winLength       = 1;
    % EEG
    if sum(contains(header.label,'C3-A2')) == 1
        EEG = filtering(data{contains(header.label,'C3-A2')},EEGfs,'eeg');
    elseif sum(contains(header.label,'C3-M2')) == 1
        EEG = filtering(data{contains(header.label,'C3-M2')},EEGfs,'eeg');
    elseif sum(strcmp(header.label,'C3')) == 1 && sum(strcmp(header.label,'A2')) == 1 
        EEG = filtering(data{strcmp(header.label,'C3')} ...
            - data{strcmp(header.label,'A2')},EEGfs,'eeg');
    elseif sum(strcmp(header.label,'C3-')) == 1
        EEG = filtering(data{strcmp(header.label,'C3-')},EEGfs,'eeg');
    else 
        EEG = data(contains(header.label,'C3'));
        EEG = filtering(EEG{1},EEGfs,'eeg');
    end
    %% ECG
    ECG                         = data(contains(header.label,{'ECG','EKG'}));
    ECGfs                       = header.samplerate(contains(header.label,{'ECG','EKG'}));
    if length(ECG) > 1, ECG = ECG(1); ECGfs = ECGfs(1); end
    ECG                         = ECG{1};
    ecgRes                      = reshape(ECG(1:end-mod(length(ECG) ...
                                    ,ECGfs*epochLength)),ECGfs*epochLength,[]);
    ecgRes                      = ecgRes(:,events.LOF:end);
    ECG                         = ecgRes(:);  
    %% Features
    EEGTEMPORAL                 = temporal_statistics(EEG,EEGfs,segLength,events);
    EEGSPECTRAL                 = spectral_statistics(EEG,EEGfs,segLength,events);
    EEGARMODEL                  = estimate_ARmodel(EEG,EEGfs,segLength,5);
    ECGHRV                      = heart_analysis(ECG,ECGfs,segLength,events);
    HYPNNOGRAM                    = hypnogram_analysis(events.hypnogram);
    %% Save
    features = [features; wscData(strcmp(wscData.ID,id),:), HYPNNOGRAM];
    writetable(features,'_features\WSChypnogram.csv');  
    writetable(EEGARMODEL,strcat('_features\AR_',id,'.csv'));  
    try
    h5create(strcat('_features\AR_',id,'.h5'),'/values',size(table2array(EEGARMODEL)))
    h5write(strcat('_features\AR_',id,'.h5'),'/values',table2array(EEGARMODEL));
    h5create(strcat('_features\AR_',id,'_chunk.h5'),'/values',size(table2array(EEGARMODEL)),'ChunkSize',[6 5])
    h5write(strcat('_features\AR_',id,'_chunk.h5'),'/values',table2array(EEGARMODEL));
    catch
    end
    catch
        errorSubj(coh)              = errorSubj(coh) + 1;
        errorSubjs{errorSubj(coh)}  = id;
        fprintf('%s is funny.',id);
    end
    fprintf('%i/%i done.',sub,nSubjPrCohort(coh));
    toc,end
elseif strcmp(curCohort,'mros')
    evtlisting      = listing(endsWith(lower(listing),'.xml'));
    mrosTimes       = readtable(strcat(curDataPath,'Times-visit1.csv'));
    mrosData        = readtable(strcat(curDataPath,'MrOSv1.csv'));
    tic,for sub = 1:nSubjPrCohort(coh) 
    try
        %% ID checking
        subjFile        = edflisting{sub}; 
        subjID          = split(subjFile,'.'); 
        id = split(subjID{1},'-'); id = upper(id{3});
        if sum(strcmp(mrosData.ID,id)) == 1
        idIdx           = find(strcmp(mrosTimes.nsrrid,id));
        subjEnd         = subjID{2}; 
        fullEDFpath     = strcat(curDataPath,subjFile);
        fullXMLpath     = strcat(curDataPath,strrep(subjFile,strcat('.',subjEnd),'-nsrr.xml'));
        %% Load data
        [header, data]  = loadEDF(fullEDFpath);
        events          = loadHypnogramMrOS(fullXMLpath,mrosTimes.Start2LightsOff(idIdx));
        for i = 1:length(header.label),fieldd=header.label{i};header.label{i}=fieldd(~isspace(fieldd));end
        EEGfs           = header.samplerate(contains(header.label,'C3')); EEGfs = EEGfs(1);
        epochLength     = 30;
        winLength       = 1;
        %% EEG
        if sum(contains(header.label,'C3-A2')) == 1
            EEG = filtering(data{contains(header.label,'C3-A2')},EEGfs,'eeg');
        elseif sum(contains(header.label,'C3-M2')) == 1
            EEG = filtering(data{contains(header.label,'C3-M2')},EEGfs,'eeg');
        elseif sum(strcmp(header.label,'C3')) == 1 && sum(strcmp(header.label,'A2')) == 1 
            EEG = filtering(data{strcmp(header.label,'C3')} ...
                - data{strcmp(header.label,'A2')},EEGfs,'eeg');
        elseif sum(strcmp(header.label,'C3')) == 1
            EEG = filtering(data{strcmp(header.label,'C3')},EEGfs,'eeg');
        else 
            EEG = data(contains(header.label,'C3'));
            EEG = filtering(EEG{1},EEGfs,'eeg');
        end
        %% ECG
        ECG                         = data(contains(header.label,{'ECG','EKG'}));
        ECGfs                       = header.samplerate(contains(header.label,{'ECG','EKG'}));
        if length(ECG) > 1, ECG = ECG(1); ECGfs = ECGfs(1); end
        ECG                         = ECG{1};
        ecgRes                      = reshape(ECG(1:end-mod(length(ECG) ...
                                        ,ECGfs*epochLength)),ECGfs*epochLength,[]);
        ecgRes                      = ecgRes(:,events.LOF:end);
        ECG                         = ecgRes(:);    
        %% Features
        EEGTEMPORAL                 = temporal_statistics(EEG,EEGfs,segLength,events);
        EEGSPECTRAL                 = spectral_statistics(EEG,EEGfs,segLength,events);
        ECGHRV                      = heart_analysis(ECG,ECGfs,segLength,events);
        %% Save
        features = [features; mrosData(strcmp(mrosData.ID,id),:), EEGTEMPORAL, EEGSPECTRAL, ECGHRV];
        writetable(features,'MrOSar.csv');
        end
        fprintf('%i/%i done.',sub,nSubjPrCohort(coh));
    catch
        errorSubj(coh)              = errorSubj(coh) + 1;
        errorSubjs{errorSubj(coh)}  = id;
        fprintf('%s is funny.',id);
    end
    toc,end
end
