function [ahi,plmi,plmai] = analyze_wsc_events(events,fs)
    plmai = [];
% Calculates AHI, PLMI, PLMAI
    TSThrs      = (events.TST/60/60);
%     AHI
    evtTypes    = {'hypopnea','apnea'};
    sss         = [0 1 2 3 5];
    ahi         = nan(5,length(evtTypes));
    for i = 1:length(evtTypes)
        ep = events.epoch(contains(lower(events.label),evtTypes{i}));
        hy = events.hypnogram(ep-events.LOF+1);
        for s = 1:5
           ahi(s,i) = sum(hy == sss(s));
        end
    end
    ahi = sum(ahi(:))/TSThrs;
%     PLMI
    ep = events.epoch(contains(lower(events.label),'lm'));
    sa = events.start_stop_matrix(contains(lower(events.label),'lm'),:);
    sa = sa'; sa = diff(sa(:)); sa = sa(2:2:end)/fs;
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
% %     PLMAI
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

