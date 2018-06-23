function out = xml2events(xmlpath,lightsoffSec)
    epochLength = 30;
    DOMnode = xml2struct(xmlpath);
    stru    = DOMnode.PSGAnnotation.ScoredEvents.ScoredEvent;
    stru    = cell2struct(stru,{'structs'});
    hyp     = nan(2000,1); % pre-allocating 2000 epochs
    idx     = 1;
    lightsoffEpoch = floor(lightsoffSec/epochLength) + 1;
    out.LOF = lightsoffEpoch;
    eventIdx= false(length(stru),1);
    for i = 1:length(stru)
        if contains(stru(i).structs.EventType.Text,'Stage')
            stageStr = split(stru(i).structs.EventConcept.Text,'|');
            stageDou = str2double(stageStr(2));
            stageDur = stru(i).structs.Duration.Text;
            nEpochs = floor(str2double(stageDur)/epochLength);            
            hyp(idx:idx+nEpochs-1) = repmat(stageDou,nEpochs,1);
            idx = idx+nEpochs;
        end
        if contains(stru(i).structs.EventConcept.Text,'Recording Start Time')
            TRS = split(stru(i).structs.ClockTime.Text,' ');
            TRS = replace(TRS{2},'.',':');
        end
        eventIdx(i) = ~contains(stru(i).structs.EventType.Text,'Stage');
    end
    hyp = hyp(~isnan(hyp));
    out.times     = datevec(datetime(TRS) + seconds((0:length(hyp)-1)*30));
    hyp = hyp(lightsoffEpoch:end);
    out.LON      = length(hyp);
    out.times = out.times(lightsoffEpoch:end,:);
%     hypIdxs = find(diff(hyp)~=0); lastSleepEpoch = hypIdxs(end);
    out.hypnogram= hyp;
    out.TST              = sum(out.hypnogram>0)*30;
    out.TSTN1            = sum(out.hypnogram>1)*30;
%     events.sleepHyp = hyp(1:lastSleepEpoch);
%     lastSleepSec    = lastSleepEpoch*epochLength;
    %% Events
%     EVENTS  = struct();
%     evtCount= zeros(5,1);
%     for i = 1:length(stru)
%         curTime = str2double(stru(i).structs.Start.Text);
%         if eventIdx(i) && curTime >= lightsoffSec && curTime < lastSleepSec
%             curEvent = stru(i);
%             type = lower(strrep(curEvent.structs.EventConcept.Text,'|',' '));
%             ss   = hyp(ceil((str2double(curEvent.structs.Start.Text) - lightsoffSec)/epochLength)+1);
%             if contains(type,{'plm','lm','limb'})
%                 evtCount(1) = evtCount(1) + 1;
%                 EVENTS.PLM(evtCount(1),1) = ss;
%                 EVENTS.PLM(evtCount(1),2) = str2double(curEvent.structs.Duration.Text);
%                 EVENTS.PLM(evtCount(1),3) = str2double(curEvent.structs.Start.Text);
%             elseif contains(type,{'arousal'})
%                 evtCount(2) = evtCount(2) + 1;
%                 EVENTS.ARO(evtCount(2),1) = ss;
%                 EVENTS.ARO(evtCount(2),2) = str2double(curEvent.structs.Duration.Text);
%                 EVENTS.ARO(evtCount(2),3) = str2double(curEvent.structs.Start.Text);
%             elseif contains(type,{'hypopnea'})
%                 evtCount(3) = evtCount(3) + 1;
%                 EVENTS.HYP(evtCount(3),1) = ss;
%                 EVENTS.HYP(evtCount(3),2) = str2double(curEvent.structs.Duration.Text);
%                 EVENTS.HYP(evtCount(3),3) = str2double(curEvent.structs.Start.Text);
%             elseif contains(type,{'apnea'})
%                 APNtypes = {'obstructive','central','mixed'};
%                 evtCount(4) = evtCount(4) + 1;
%                 EVENTS.APN(evtCount(4),1) = ss;
%                 EVENTS.APN(evtCount(4),2) = str2double(curEvent.structs.Duration.Text);
%                 EVENTS.APN(evtCount(4),3) = str2double(curEvent.structs.Start.Text);
%             elseif contains(type,{'spo2', 'desaturation'})
%                 evtCount(5) = evtCount(5) + 1;
%                 EVENTS.SP2(evtCount(5),1) = ss;
%                 EVENTS.SP2(evtCount(5),2) = str2double(curEvent.structs.Duration.Text);
%                 EVENTS.SP2(evtCount(5),3) = str2double(curEvent.structs.Start.Text);
%             end
%         end
%     end
%     
%     s2 = EVENTS.ARO(:,3);
%     e2 = sum(EVENTS.ARO(:,2:3),2);
%     for plm = 1:length(EVENTS.PLM)
%         s1 = EVENTS.PLM(plm,3);
%         e1 = EVENTS.PLM(plm,3) + EVENTS.PLM(plm,2);
%         curAsscApn = find((s2 - 0.5 <= e1).*(e2 >= s1 - 0.5));
%         EVENTS.PLM(plm,4) = sum((s2 - 0.5 <= e1).*(e2 >= s1 - 0.5));
%     end
%     
%     fi = fields(EVENTS);
%     events.Events = nan(5,25);
%     qs = [0.25 0.50 0.75];
%     sss= [0 1 2 3 5];
%     for i = 1:length(fi)
%         ii = 1;
%         ss = 1;
%         while ii+4 <= size(events.Events,2)
%             events.Events(i,ii)      = length(EVENTS.(fi{i})(EVENTS.(fi{i})(:,1) == sss(ss),2));
%             events.Events(i,ii+1)      = length(EVENTS.(fi{i})(EVENTS.(fi{i})(:,1) == sss(ss),2)) / ...
%                 sum(hyp == sss(ss));
%             events.Events(i,(ii+2):(ii+4))    = quantile(EVENTS.(fi{i})(EVENTS.(fi{i})(:,1) == sss(ss),2),qs); 
%             ii = ii + 5;
%             ss = ss + 1;
%         end
%     end
%     events.Events(isnan(events.Events)) = 0;
end