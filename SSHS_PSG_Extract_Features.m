startup
%% Main loop
pth             = 'Data\shhs\';
listing         = extractfield(dir(pth), 'name')';
edflisting      = listing(endsWith(lower(listing),'.edf'));
evtlisting      = listing(endsWith(lower(listing),'.xml'));
shhsData        = readtable('SHHSfutureAll.csv');
shhsData.ID     = string(shhsData.ID);
segLength       = 5;
features        = table();
tic,for i = 1:length(edflisting)
try
    %% ID checking
    subjFile        = edflisting{i}; 
    subjID          = split(subjFile,'.'); 
    id = split(subjID{1},'-'); id = upper(id{2});
    if sum(strcmp(shhsData.ID,id)) == 1
    fullEDFpath     = strcat(pth,subjFile);
    fullXMLpath     = strcat(pth,strrep(subjFile,'.edf','-nsrr.xml'));
    %% Load data
    [header, data]  = loadEDF(fullEDFpath);
    lights          = data{contains(lower(header.label),'ligh')};
    lightsFs        = header.samplerate(contains(lower(header.label),'ligh'));
    LOFsecs         = find(diff(lights) == 1);
    if length(LOFsecs) > 1, LOFsec = LOFsecs((diff(LOFsecs)/lightsFs)>(60*10))+1; 
        if isempty(LOFsec), LOFsec = LOFsecs(end); else, LOFsec = LOFsec(1); end
    else LOFsec = 1;
    end
    events          = loadHypnogramSHHS(fullXMLpath,LOFsec);
    for ii = 1:length(header.label),fieldd=header.label{ii};header.label{ii}=fieldd(~isspace(fieldd));end
    EEGfs           = header.samplerate(contains(header.label,{'C3','EEG'})); EEGfs = EEGfs(1);
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
        EEG = data(contains(header.label,'EEG'));
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
    EEGTEMPORAL                 = temporal_statistics(EEG,EEGfs,segLength,events,false);
    EEGSPECTRAL                 = spectral_statistics(EEG,EEGfs,segLength,events,false);
    EEGARMODEL                  = estimate_ARmodel(EEG,EEGfs,segLength,5);
    try
    h5create(strcat('F:\DTU\Stanford\Data\_features\SHHS\',id,'.h5'),'/demographics',size(table2array(shhsData(strcmp(shhsData.ID,id),2:end))));
    h5write(strcat('F:\DTU\Stanford\Data\_features\SHHS\',id,'.h5'),'/demographics',shhsData{strcmp(shhsData.ID,id),2:end});
    h5create(strcat('F:\DTU\Stanford\Data\_features\SHHS\',id,'.h5'),'/temporal',size(table2array(EEGTEMPORAL)),'ChunkSize',[6 13])
    h5write(strcat('F:\DTU\Stanford\Data\_features\SHHS\',id,'.h5'),'/temporal',table2array(EEGTEMPORAL));
    h5create(strcat('F:\DTU\Stanford\Data\_features\SHHS\',id,'.h5'),'/spectral',size(table2array(EEGSPECTRAL)),'ChunkSize',[6 8])
    h5write(strcat('F:\DTU\Stanford\Data\_features\SHHS\',id,'.h5'),'/spectral',table2array(EEGSPECTRAL));
    h5create(strcat('F:\DTU\Stanford\Data\_features\SHHS\',id,'.h5'),'/armodel',size(table2array(EEGARMODEL)),'ChunkSize',[6 5])
    h5write(strcat('F:\DTU\Stanford\Data\_features\SHHS\',id,'.h5'),'/armodel',table2array(EEGARMODEL));
    catch
    end
    HYPNOGRAM                   = hypnogram_analysis(events.hypnogram);
    ECGHRV                      = heart_analysis(ECG,ECGfs,30,events);
    EEGTEMPORAL                 = temporal_statistics(EEG,EEGfs,30,events,true);
    EEGSPECTRAL                 = spectral_statistics(EEG,EEGfs,30,events,true);
    %% Save
    features = [features; shhsData(strcmp(shhsData.ID,id),:),HYPNOGRAM, EEGTEMPORAL, EEGSPECTRAL, ECGHRV];
    writetable(features,'SHHSfeatures.csv');
    end
    fprintf('%i/%i done.',i,length(edflisting));
catch
    fprintf('%s is funny.',id);
end
toc,end