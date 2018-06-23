function [ahi,plmi,plmai] = analyze_ssc_events(events)
    plmai = [];
% Calculates AHI, PLMI, PLMAI
    TSThrs      = (events.TST/60/60);
    lof         = find(contains(lower(events.label),{'lights off'}));
    lon         = find(contains(lower(events.label),{'lights on'}));
    evts        = events.label(lof + 1: lon - 1);
%     AHI
    evtTypes    = {'hypopnea','apnea'};
    TSThrs      = (events.TST/60/60);
    ahi         = sum(contains(lower(evts),evtTypes))/TSThrs;
%     PLMI
    ep = events.epoch(contains(lower(events.category),'plm'));
    sa = events.startStopSamples(contains(lower(events.category),'plm'),:);
    sa = sa'; sa = diff(sa(:)); sa = sa(2:2:end)/events.samplerate;
    plmCount    = 0;
	curCount    = 0;
    idx         = 1;
    lastSize    = 0;
    for i = 1:length(sa)
        if sa(i) <= 90
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
%     ep = events.epoch(contains(lower(events.label),'lma'));
%     hy = events.hypnogram(ep-events.LOF+1);
%     sa = events.start_stop_matrix(contains(lower(events.label),'lma'),:);
%     sa = sa'; sa = diff(sa(:)); sa = sa(2:2:end)/fs;
%     plmCount    = 0;
%     idx         = 1;
%     lastSize    = 0;
%     for i = 1:length(sa)
%         if sa(i) <= 90
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
end

