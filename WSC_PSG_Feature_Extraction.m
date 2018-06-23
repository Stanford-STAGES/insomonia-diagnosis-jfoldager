startup
%% Main loop
pth             = 'Data\wsc\psgs\';
listing         = extractfield(dir(pth), 'name')';
edflisting      = listing(endsWith(lower(listing),'.edf'));
evtlisting      = listing(endsWith(lower(listing),'.sta'));
demoData        = readtable('WSCall.csv');
demoData.ID     = string(demoData.ID);
demoData        = [demoData(:,1:5),demoData(:,10:13),demoData(:,17),demoData(:,21)];
segLength       = 5;
sizes           = zeros(length(edflisting),3);
savepth         = 'Data\_features\WSC\5sec\';
load FeatureNames.mat
features        = array2table(nan(length(edflisting),549),'VariableNames',names);
features.ID     = cell(size(features,1),1);
tic,for i = 1:length(edflisting)
try
    %% ID checking
    subjFile        = edflisting{i}; 
    subjID          = split(subjFile,' '); 
    id = subjID{1}; 
    if sum(strcmp(demoData.ID,id)) == 1
    fullEDFpath     = strcat(pth,subjFile);
    fullSTApath     = strcat(pth,strrep(subjFile,'.EDF','.STA'));
    %% Load data
    [header, data]  = loadEDF(fullEDFpath);
    events          = loadEVT(fullSTApath);
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
    elseif sum(strcmp(header.label,'C3-x')) == 1
        EEG = filtering(data{strcmp(header.label,'C3-x')},EEGfs,'eeg');
    elseif sum(contains(header.label,'C3')) == 1
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
    s = [size(EEGTEMPORAL,1);size(EEGSPECTRAL,1);size(EEGARMODEL,1)];
    for ii = 1:3
        sizes(i,ii) = s(ii);
    end
    try
    h5create(strcat(savepth,id,'.h5'),'/demographics',size(table2array(demoData(strcmp(demoData.ID,id),2:end))));
    h5write(strcat(savepth,id,'.h5'),'/demographics',demoData{strcmp(demoData.ID,id),2:end});
    h5create(strcat(savepth,id,'.h5'),'/temporal',size(table2array(EEGTEMPORAL)),'ChunkSize',[6 13])
    h5write(strcat(savepth,id,'.h5'),'/temporal',table2array(EEGTEMPORAL));
    h5create(strcat(savepth,id,'.h5'),'/spectral',size(table2array(EEGSPECTRAL)),'ChunkSize',[6 8])
    h5write(strcat(savepth,id,'.h5'),'/spectral',table2array(EEGSPECTRAL));
    h5create(strcat(savepth,id,'.h5'),'/armodel',size(table2array(EEGARMODEL)),'ChunkSize',[6 5])
    h5write(strcat(savepth,id,'.h5'),'/armodel',table2array(EEGARMODEL));
    catch
    end
    HYPNOGRAM                   = hypnogram_analysis(events.hypnogram);
    ECGHRV                      = heart_analysis(ECG,ECGfs,30,events);
    EEGTEMPORAL                 = temporal_statistics(EEG,EEGfs,30,events,true);
    EEGSPECTRAL                 = spectral_statistics(EEG,EEGfs,30,events,true);
    %% Save
    features(i,2:end) = [demoData(strcmp(demoData.ID,id),2:end),HYPNOGRAM, EEGTEMPORAL, EEGSPECTRAL, ECGHRV];
    features.ID{i} = id;
    writetable(features,'WSCfeatures.csv');
    end
    fprintf('%i/%i done.',i,length(edflisting));
catch
    fprintf('%s is funny.',id);
end
toc,end