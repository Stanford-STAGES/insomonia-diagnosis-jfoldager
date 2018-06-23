function events = loadHypnogramSHHS(file,start2lof)
            
DOMnode = xml2struct(file);
stru    = DOMnode.PSGAnnotation.ScoredEvents.ScoredEvent;
stru    = cell2struct(stru,{'structs'});
events.hypnogram = nan(2000,1); % pre-allocating 2000 epochs
idx     = 1;
T                   = table();
T.TST               = nan(1,1);
T.SOL               = nan(1,1);
T.SE                = nan(1,1);
T.AW25              = nan(1,1);
T.AW5               = nan(1,1);
T.AWN125            = nan(1,1);
T.AWN15             = nan(1,1);
T.WA                = nan(1,1);
T.N1                = nan(1,1);
T.N2                = nan(1,1);
T.N3                = nan(1,1);
T.REM               = nan(1,1);
T.REML              = nan(1,1);
% T.AHI               = nan(1,1);
T.PLMI              = nan(1,1);
nAHIs   = 0;
nPLMs   = 0;
start2lof = ceil(start2lof/30);
for i = 1:length(stru)
    if contains(stru(i).structs.EventConcept.Text,'Apnea') || ...
        contains(stru(i).structs.EventConcept.Text,'Hypopnea')
    nAHIs = nAHIs + 1;
    end
    if contains(stru(i).structs.EventConcept.Text,'PLM') && ...
            str2double(stru(i).structs.Start.Text)>=start2LOF(nSubject)
        nPLMs = nPLMs + 1;
    end
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
SOP                 = find((events.hypnogram>0).*(events.hypnogram<7) == 1); SOP=SOP(1);
sleephyp            = events.hypnogram(SOP:end);
T.TST       = sum(events.hypnogram>0)*30/60;
T.SOL       = (SOP-1)*30/60;
T.SE        = sum((sleephyp > 0).*(sleephyp<7) == 1)/length(sleephyp);
T.WA        = sum(sleephyp == 0)/length(sleephyp);
T.N1        = sum(sleephyp == 1)/length(sleephyp);
T.N2        = sum(sleephyp == 2)/length(sleephyp);
T.N3        = sum(sleephyp == 3)/length(sleephyp);
T.REM       = sum(sleephyp == 5)/length(sleephyp);
try a  = find(sleephyp == 5); T.REML = (a(1)-1)*30; catch, end
a = diff(find(diff(diff(cumsum((sleephyp == 0) + (sleephyp == 7)))) ~= 0));
a = a(1:2:end);
T.AW25 = sum(a>=(2.5*60)/30)/(T.TST/60);
T.AW5  = sum(a>=(5*60)/30)/(T.TST/60);
a = diff(find(diff(diff(cumsum((sleephyp == 0) + (sleephyp == 1) + (sleephyp == 7)))) ~= 0));
a = a(2:2:end);
T.AWN125  = sum(a>=(2.5*60)/30)/(T.TST/(60));
T.AWN15    = sum(a>=(5*60)/30)/(T.TST/(60));     

try REML = find(events.hypnogram==5); T.REML = (REML(1)-1)*30/60; catch, REML = nan; end
% T.AHI = nAHIs/(T.TST/60);
T.PLMI = nPLMs/(T.TST/60);
events.T = T;
end