startup
pth = 'mros-visit1-dataset-0.3.0.csv';
TData = readtable(pth, ...
    detectImportOptions(pth));
tStudyStart = datetime(string(TData.poststtp),'Format','HH:mm:ss');
tStudyStart = datevec(tStudyStart(~isnat(tStudyStart),:)); 
tLightsOff  = datetime(string(TData.postlotp),'Format','HH:mm:ss');
tLightsOff  = datevec(tLightsOff(~isnat(tLightsOff),:)); 
tSleepOnset = datetime(string(TData.postontp),'Format','HH:mm:ss');
tSleepOnset = datevec(tSleepOnset(~isnat(tSleepOnset),:)); 
tStudyStop  = datetime(string(TData.postendp),'Format','HH:mm:ss');
tStudyStop  = datevec(tStudyStop(~isnat(tStudyStop),:));
start2LOF   = timediff(tStudyStart,tLightsOff);
LOF2SO      = timediff(tLightsOff,tSleepOnset);
LOF2Stop    = timediff(tLightsOff,tStudyStop);
start2Stop  = timediff(tStudyStart,tStudyStop);
lightsoffEpochs = floor(start2LOF/30) + 1;

extDiskHypPath  = strcat('visit1\');
listing         = extractfield(dir(extDiskHypPath), 'name')';
listing         = listing(contains(listing,'.xml'));
nSubjects       = length(listing);
nSubject        = 1;

T                   = table();
T.ID                = cell(nSubjects,1);
T.TST               = nan(nSubjects,1);
T.SOL               = nan(nSubjects,1);
T.SE                = nan(nSubjects,1);
T.AW25              = nan(nSubjects,1);
T.AW5               = nan(nSubjects,1);
T.AWN125            = nan(nSubjects,1);
T.AWN15             = nan(nSubjects,1);
T.WA                = nan(nSubjects,1);
T.N1                = nan(nSubjects,1);
T.N2                = nan(nSubjects,1);
T.N3                = nan(nSubjects,1);
T.REM               = nan(nSubjects,1);
T.REML              = nan(nSubjects,1);
T.AHI               = nan(nSubjects,1);
T.PLMI              = nan(nSubjects,1);
tic,for nSubject = 1063:nSubjects
    id = split(listing{nSubject},'-'); id = upper(id{3});
file = strcat(extDiskHypPath,listing{nSubject});
try
DOMnode = xml2struct(file);
stru    = DOMnode.PSGAnnotation.ScoredEvents.ScoredEvent;
stru    = cell2struct(stru,{'structs'});
hyp     = nan(2000,1); % pre-allocating 2000 epochs
idx     = 1;
T.ID{nSubject}= id;
curSubj = nSubject;
nAHIs = 0;
nPLMs = 0;
for i = 1:length(stru)
    if contains(stru(i).structs.EventConcept.Text,'Apnea') || ...
            contains(stru(i).structs.EventConcept.Text,'Hypopnea')
        nAHIs = nAHIs + 1;
    end
    if contains(stru(i).structs.EventConcept.Text,'PLM') && ...
            str2double(stru(i).structs.Start.Text)>=start2LOF(nSubject)
        nPLMs = nPLMs + 1;
    end
    if contains(stru(i).structs.EventType.Text,'Stage')
        stageStr = split(stru(i).structs.EventConcept.Text,'|');
        stageDou = str2double(stageStr(2));
        stageDur = stru(i).structs.Duration.Text;
        nEpochs = floor(str2double(stageDur)/30);            
        hyp(idx:idx+nEpochs-1) = repmat(stageDou,nEpochs,1);
        idx = idx+nEpochs;
    end
    eventIdx(i) = ~contains(stru(i).structs.EventType.Text,'Stage');
end
hyp = hyp(~isnan(hyp));
hyp = hyp(lightsoffEpochs(nSubject):end);
SOP            = find((hyp>0).*(hyp<7) == 1); SOP=SOP(1);
sleephyp        = hyp(SOP:end);
T.TST(curSubj)       = sum(hyp>0)*30;
T.SOL(curSubj)       = (SOP-1)*30;
T.SE(curSubj)        = sum((sleephyp > 0).*(sleephyp<7) == 1)/length(sleephyp);
T.WA(curSubj)        = sum(sleephyp == 0)/length(sleephyp);
T.N1(curSubj)        = sum(sleephyp == 1)/length(sleephyp);
T.N2(curSubj)        = sum(sleephyp == 2)/length(sleephyp);
T.N3(curSubj)        = sum(sleephyp == 3)/length(sleephyp);
T.REM(curSubj)       = sum(sleephyp == 5)/length(sleephyp);
try a  = find(sleephyp == 5); T.REML(curSubj) = (a(1)-1)*30; catch, end
a = diff(find(diff(diff(cumsum((sleephyp == 0) + (sleephyp == 7)))) ~= 0));
a = a(1:2:end);
T.AW25(curSubj) = sum(a>=(2.5*60)/30)/(T.TST(curSubj)/60/60);
T.AW5(curSubj)  = sum(a>=(5*60)/30)/(T.TST(curSubj)/60/60);
a = diff(find(diff(diff(cumsum((sleephyp == 0) + (sleephyp == 1) + (sleephyp == 7)))) ~= 0));
a = a(2:2:end);
T.AWN125(curSubj)  = sum(a>=(2.5*60)/30)/(T.TST(curSubj)/60/60);
T.AWN15(curSubj)    = sum(a>=(5*60)/30)/(T.TST(curSubj)/60/60);     

try REML = find(hyp==5); T.REML(curSubj) = (REML(1)-1)*30; catch, REML = nan; end
T.AHI(curSubj) = nAHIs/(T.TST(curSubj)/60/60);
T.PLMI(curSubj) = nPLMs/(T.TST(curSubj)/60/60);
fprintf('%i out of %i is done.\n',nSubject,nSubjects);
catch
    fprintf('%s SUCKS.\n',id);
end
toc,end
writetable(T,'MrOSHypnogramsV1.csv');