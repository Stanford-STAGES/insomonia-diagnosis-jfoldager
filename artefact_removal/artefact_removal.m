function PSG = artefact_removal(data,header,events,epochLength,cohort)
    lightsOffEpoch = events.LOF;
    lightsOnEpoch = events.LON;
    if contains(cohort,'mros')
        C3Ref = data{contains(header.label,'C3Ref')}; 
        fs= header.samplerate(contains(header.label,'C3'));
        C3Ref = reshape(C3Ref(1:end-mod(length(C3Ref) ...
                ,fs(1)*epochLength)),fs(1)*epochLength,[]);   
        C3Ref = C3Ref(:,lightsOffEpoch:lightsOnEpoch); 
        try
            STAT = data{contains(header.label,'STAT')}; 
            fs= header.samplerate(contains(header.label,'STAT'));
            STAT = reshape(STAT(1:end-mod(length(STAT) ...
                    ,fs*epochLength)),fs*epochLength,[]);
            STAT = STAT(:,lightsOffEpoch:lightsOnEpoch);
        catch
            STAT = zeros(size(C3Ref));
        end
    %%   Requirements for it to be an artefact:
        status    =   sum(STAT)...
                    + sum(abs(C3Ref)>=250); 
        crSig = C3Ref;
    elseif contains(cohort,'wsc')
    %%   Requirements for it to be an artefact:  
        str = {'F3','FZ','CZ','C3','PZ','O1'};
        first = true;
        for i = 1:length(data)
            if contains(header.label{i},str)
                fs      = header.samplerate(i);
                crSig   = reshape(data{i}(1:end-mod(length(data{i}) ...
                ,fs*epochLength)),fs*epochLength,[]); crSig = crSig(:,lightsOffEpoch:lightsOnEpoch); 
                if first
                    status = zeros(1,size(crSig,2));
                    first = false;
                end
                status  = status + sum(abs(crSig)>=250)>0;              
            end
        end
        if length(status) == length(data{end})
            status = status + (data{end} == 7)';
        elseif length(status) > length(data{end})
            data{end} = [0;data{end}];
        elseif length(status) < length(data{end})
            data{end} = data{end};
        end
    elseif contains(cohort,'ssc')
    %%   Requirements for it to be an artefact:  
        fs = header.samplerate(contains(header.label,'C3Ref'));
        crSig = reshape(data{contains(header.label,'C3Ref')}(1:end-mod(length(data{contains(header.label,'C3Ref')}) ...
            ,fs*epochLength)),fs*epochLength,[]); 
        status = sum(crSig(:,lightsOffEpoch:lightsOnEpoch) > 250) ; 
        
        status = status + (data{end} == 7)';
    end
    if length(data{end})  > length(status), data{end} = data{end}(1:length(status)); end
    if length(status)  > length(data{end}), data{end} = [data{end}; 0]; end
    t = linspace(0,length(crSig(:)),length(crSig(:)))*1/fs(1);
    t = reshape(t(1:end-mod(length(t) ...
            ,fs(1)*epochLength)),fs(1)*epochLength,[]); 
    isArtefact    = status > 0;
    t = t(:,isArtefact == 0);
    times = events.times;
%%  Artefact Removal
    PSG = struct();
    for i = 1:length(data)
        signal      = data{i}; fs = header.samplerate(i);
        signalRes   = reshape(signal(1:end-mod(length(signal),fs*epochLength))...
            ,fs*epochLength,[]);
        fieldd = header.label{i};
        if sum(strncmpi(cellstr(('0':'9')')',fieldd,1)) == 1, fieldd = strcat('eeg_',fieldd);end
        if ~contains(header.label(i),'Hypnogram'), signalRes   = signalRes(:,lightsOffEpoch:lightsOnEpoch);end
        fieldd(~ismember(fieldd,['A':'Z' 'a':'z' '0':'9'])) = '_';
        PSG.(strrep(fieldd(~isspace(fieldd)),' ',''))  = signalRes(:,isArtefact == 0);
    end
    PSG.t                   = t;
    PSG.times               = times(isArtefact == 0,:);
    PSG.Org_ArtefactEpochs  = isArtefact;
    org_LastSleepEpoch      = find(((data{contains(header.label,'Hypnogram')}> 0).*...
        data{contains(header.label,'Hypnogram')}< 7) ==1); 
    new_LastSleepEpoch      = find(((PSG.Hypnogram > 0).*(PSG.Hypnogram < 7))==1); 
    PSG.Org_LastSleepEpoch  = org_LastSleepEpoch(end);
    PSG.New_LastSleepEpoch  = new_LastSleepEpoch(end);
    PSG.Org_LightsOffEpoch  = lightsOffEpoch;
    PSG.New_Artefacts = isArtefact(PSG.Org_LightsOffEpoch:PSG.Org_LastSleepEpoch);
end