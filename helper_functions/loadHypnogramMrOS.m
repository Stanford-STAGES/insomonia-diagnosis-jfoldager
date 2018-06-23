function events = loadHypnogramMrOS(file,start2lof)
            
DOMnode = xml2struct(file);
stru    = DOMnode.PSGAnnotation.ScoredEvents.ScoredEvent;
stru    = cell2struct(stru,{'structs'});
events.hypnogram = nan(2000,1); % pre-allocating 2000 epochs
idx     = 1;
nAHIs   = 0;
nPLMs   = 0;
T       = table();
start2lof = floor(start2lof/30);
for i = 1:length(stru)
    if  contains(stru(i).structs.EventType.Text,'Stage')
        stageStr = split(stru(i).structs.EventConcept.Text,'|');
        stageDou = str2double(stageStr(2));
        stageDur = stru(i).structs.Duration.Text;
        nEpochs = floor(str2double(stageDur)/30);            
        events.hypnogram(idx:idx+nEpochs-1) = repmat(stageDou,nEpochs,1);
        idx = idx+nEpochs;
    end
end
events.hypnogram    = events.hypnogram(~isnan(events.hypnogram));
events.hypnogram    = events.hypnogram(start2lof:end);
events.LOF          = start2lof;
events.LON          = length(events.hypnogram);
sLastSleep          = find(events.hypnogram ~= 0); events.MW = sLastSleep(end);
end