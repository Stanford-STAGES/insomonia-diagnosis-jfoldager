function [out,ahi,plmi,plmai] = csv2events(fullpath)
    plmai = [];
    epochLength = 30;
    fid = fopen(fullpath);
    nSemicolons = textscan(fid,'%s',1); nSemicolons = count(nSemicolons{1},';');
    characters = repmat('%s ',1,nSemicolons + 1);
    C = textscan(fid,characters,'Delimiter',';','TreatAsEmpty',{'NA','na'});
    fclose(fid);
%% Times
    startTime       = C{end-1}(1);
    if etime(datevec(C{end-1}(2)),datevec(startTime)) < 0,startTime = C{end-1}(2);end
    lof                 = find(contains(lower(C{end}),'lights out')); lof = lof(1);
    lightsOffTime   = C{end-1}(contains(lower(C{end}),'lights out'));
    T0                  = split(startTime,'.');T0=T0{1};
    TLOF                = split(lightsOffTime,'.');TLOF=TLOF{1};
    out.LOF  = ceil(timediff(datevec(datetime(T0)),datevec(datetime(TLOF)))/epochLength);
    T0                  = datevec(TLOF);  
    lightsOnTime    = C{end-1}(contains(lower(C{end}),'lights on'));
    lastStageTime = C{end-1}(contains(lower(C{end}),'stage')); 
    lastStageTime = datevec(lastStageTime{end});
    endTime = C{end-1}(contains(lower(C{end}),'paused')); 
    endTime = endTime(end);
    if etime(datevec(endTime),lastStageTime) < 0        
        endTime         = C{end-1}(end);
    else
        try
            endTime     = C{end-1}(contains(lower(C{end}),'paused')); endTime = endTime(end);
        catch
            endTime     = C{end-1}(end);
        end
    end    
    timediffs           = [timediff(datevec(startTime{:}),datevec(lightsOffTime{:}))...
                            timediff(datevec(startTime{:}),datevec(lightsOnTime{:}))...
                            timediff(datevec(startTime{:}),datevec(endTime{:}))];
    setPointEpoch       = ceil(timediffs(1)/epochLength) + 1;
    lon                 = ceil(timediffs(2)/epochLength) + 1;
    endTimeEpoch        = ceil(timediffs(3)/epochLength) + 1;
%% Hypnogram
    criterion           = (contains(lower(C{end}),'stage')).*(~contains(lower(C{end}),'no stage')) == 1;
    criterion           = criterion(lof:end);
    evts                = C{end}(lof:end);
    ts                  = C{end-1}(lof:end);
    allStages           = evts(criterion); 
    allStages           = strtrim(split(allStages,'-')); allStages = allStages(:,2);
    allStageTimes       = ts(criterion);
    startTimes          = repmat(datevec(startTime{:}),length(allStageTimes),1);
    timediffs           = timediff(startTimes,datevec(allStageTimes));
    allSleepEpochs      = floor(timediffs/epochLength)+1;
    allSleepEpochs      = [allSleepEpochs; endTimeEpoch];
    sleepStageDiff      = diff(allSleepEpochs);
    sleepStageNames     = {'W','N1','N2','N3','N3','R'};
    hypnogramConcat     = zeros(size(allStages));
    for i = 1:length(sleepStageNames)
        hypnogramConcat(contains(allStages,sleepStageNames{i})) = i - 1; 
    end
    hyp = nan(2000,1);
    idx = 1;
    hypnogramConcat = [0;hypnogramConcat];
    sleepStageDiff  = [allSleepEpochs(1)-setPointEpoch;sleepStageDiff];
    for i = 1:length(sleepStageDiff)
        hyp(idx:idx+sleepStageDiff(i)-1) = repmat(hypnogramConcat(i),sleepStageDiff(i),1);
        idx = idx + sleepStageDiff(i);
    end
    hyp(isnan(hyp))     = [];
    out.hypnogram       = hyp;
    out.LON             = lon;
    out.times           = datevec(datetime(T0)+seconds((0:size(out.hypnogram,1)-1)*30));
    out.TST              = sum(out.hypnogram>0)*30;
    out.TSTN1            = sum(out.hypnogram>1)*30;  
%% Events
    allEvents       = C{end}(~(contains(lower(C{end}),{'stage','paused'})));
    allEventTimes   = C{end-1}(~contains(lower(C{end}),{'stage','paused'}));
    startTimes      = repmat(datevec(startTime{:}),length(allEventTimes),1);
    timediffs       = timediff(startTimes,datevec(allEventTimes));
    allEventEpochs  = floor(timediffs/epochLength);  
    criterion       = (allEventEpochs>=lof).*(allEventEpochs<=lon) > 0;
    events          = allEvents(criterion);
    timediffs       = timediffs(criterion);
    allEventEpochs  = allEventEpochs(criterion);
%     AHI
    evtTypes    = {'hypopnea','apnea'};
    TSThrs      = (out.TST/60/60);
    ahi         = sum(contains(lower(events),evtTypes))/TSThrs;
%     PLMI
    criterion   = (contains(lower(events),'lm')).*(~contains(lower(events),'arousal')) > 0;
    sa  = timediffs(criterion);
    eventsLM    = events(criterion);
    
    a  = split(eventsLM,'Dur:');
    b  = strrep(split(a(:,2),'sec'),' ','');
    
    for i = 1:length(eventsLM)
        sa(i,2) = sa(i,1) + str2double(b{i,1});
    end
    sa = sa'; sa = sa(:); sa = diff(sa); sa = sa(2:2:end);
    plmCount    = 0;
    curCount    = 0;
    idx         = 1;
    lastSize    = 0;
    for i = 1:length(sa)
        if sa(i) <= 90
            disp(sa(i))
            curCount = curCount + 1;
            if curCount >= 3
               plmCount(idx) = curCount + 1;
            end
        else 
            curCount = 0;
            if length(plmCount) > lastSize
                idx = idx + 1; 
                lastSize = length(plmCount); 
            end
        end
    end
    plmi = sum(plmCount(:))/TSThrs;
%     PLMAI
%     criterion   = (contains(lower(events),'lm')).*(contains(lower(events),'arousal')) > 0;
%     sa  = timediffs(criterion);
%     eventsLM    = events(criterion);
%     
%     a  = split(eventsLM,'Dur:');
%     b  = strrep(split(a(:,2),'sec'),' ','');
%     
%     for i = 1:length(eventsLM)
%         sa(i,2) = sa(i,1) + str2double(b{i,1});
%     end
%     sa = sa'; sa = sa(:); sa = diff(sa); sa = sa(2:2:end);
%     plmCount    = 0;
%     curCount    = 0;
%     idx         = 1;
%     lastSize    = 0;
%     for i = 1:length(sa)
%         if sa(i) <= 90
%             disp(sa(i))
%             curCount = curCount + 1;
%             if curCount >= 3
%                plmCount(idx) = curCount + 1;
%             end
%         else 
%             curCount = 0;
%             if length(plmCount) > lastSize
%                 idx = idx + 1; 
%                 lastSize = length(plmCount); 
%             end
%         end
%     end
%     plmai = sum(plmCount(:))/TSThrs;
%     eventStages     = zeros(size(allEventEpochs));
%     
%     sleepStageNames     = {'W','N1','N2','N3','N3','R'};
%     for i = 1:length(allEventEpochs)
%         idxs = find(allEventEpochs(i)>=allSleepEpochs); 
%         if isempty(idxs) || idxs(end) > length(hypnogramConcat), continue; end
%         idx = idxs(end);
%         eventStages(i)  = hypnogramConcat(idx);
%         ss              = lower(sleepStageNames(eventStages(i)+1)); 
%         curEvent        = strtrim(split(events{i},'-'));   
%         try
%             if length(curEvent) == 4
%                 if contains(lower(curEvent(1)),'respiratory')
%                     supType = lower(strrep(curEvent(1),' ','_')); 
%                     subType = lower(strrep(curEvent(3),' ','_'));
%                     type = strcat(supType,'_',subType);
%                     eventAppearance.(type{:})(i,1) = sscanf(curEvent{2},'Dur: %d');
%                     eventAppearance.(type{:})(i,2) = sscanf(curEvent{4},'Desat %d');
%                     type = strcat(type,'_',ss);                    
%                     eventAppearance.(type{:})(i,1) = sscanf(curEvent{2},'Dur: %d');
%                     eventAppearance.(type{:})(i,2) = sscanf(curEvent{4},'Desat %d');
%                 elseif contains(lower(curEvent(1)),'desaturation')
%                     type = lower(curEvent(1));
%                     eventAppearance.(type{:})(i,1) = sscanf(curEvent{2},'Dur: %d');
%                     eventAppearance.(type{:})(i,2) = sscanf(curEvent{3},'Min %d');
%                     eventAppearance.(type{:})(i,3) = sscanf(curEvent{4},'Drop %d');
%                     type = strcat(type,'_',ss);   
%                     eventAppearance.(type{:})(i,1) = sscanf(curEvent{2},'Dur: %d');
%                     eventAppearance.(type{:})(i,2) = sscanf(curEvent{3},'Min %d');
%                     eventAppearance.(type{:})(i,3) = sscanf(curEvent{4},'Drop %d');
%                 end
%             elseif length(curEvent) == 3
%                 supType = lower(strrep(curEvent(1),' ','_')); 
%                 subType = lower(strrep(curEvent(3),' ','_'));
%                 type = strcat(supType,'_',subType);
%                 eventAppearance.(type{:})(i,1) = sscanf(curEvent{2},'Dur: %d');
%                 type = strcat(type,'_',ss);
%                 eventAppearance.(type{:})(i,1) = sscanf(curEvent{2},'Dur: %d');
%             end
%         catch
% %             fprintf('Could not get %s \n subject %s',events{i},fullpath)
%         end
%     end
%     if exist('eventAppearance','var')
%         allFields = fields(eventAppearance);
%         for i = 1:length(allFields)
%             data    = eventAppearance.(allFields{i});
%             data( ~any(data,2), : ) = [];  
%             eventAppearance.(allFields{i}) = data;
%             out.(strcat(allFields{i},'_n')) = size(data,1);
%             if size(data,1) > 1
%                 q       = quantile(data,[0.25,0.5,0.75]);
%                 if size(data,2) < 3, q = reshape(q,length(q),[]); end
%                 for ii = 1:size(q,2)
%                     out.(strcat(allFields{i},'_q1'))(ii) = q(1,ii);
%                     out.(strcat(allFields{i},'_q2'))(ii) = q(2,ii);
%                     out.(strcat(allFields{i},'_q3'))(ii) = q(3,ii);
%                 end
%             end
%         end  
%     else
% %         fprintf('NO DATA for subject %s',fullpath)
%     end
end