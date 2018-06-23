function out = loadEVT(STApath)%,EVTpath,T0)
%   Load stages
    fileID  = fopen(STApath);
        STA     = textscan(fileID, '%f%f%f', 'Delimiter', '\t', 'EmptyValue', -1);
    fclose(fileID);
%%  Hypnogram and hypnogram analysis
%   Extract hypnogram and timestamps
    hypnogram           = STA{2}; hypnogram(hypnogram==4) = 3;
    lightsOffEpoch      = find(diff(hypnogram == 7) == -1); 
    if isempty(lightsOffEpoch), lightsOffEpoch = 1; else, lightsOffEpoch = lightsOffEpoch(1)+1; end
    if hypnogram(end,1) == 7
        lightsOnEpoch   = find(diff(hypnogram == 7) == 1); 
        lightsOnEpoch   = lightsOnEpoch(end);
    else
        lightsOnEpoch   = length(hypnogram);
    end
    out.hypnogram        = hypnogram(lightsOffEpoch:lightsOnEpoch);
    out.LOF   = lightsOffEpoch;
    out.LON   = lightsOnEpoch;
    sLastSleep    = find(out.hypnogram ~= 0); out.MW = sLastSleep(end);
%     times     = datevec(datetime(T0) + seconds((0:size(out.hypnogram,1)-1)*30));
%     times     = datevec(datetime(times(lightsOffEpoch,:)) + seconds((0:size(out.hypnogram,1)-1)*30));
%     out.times            = times;
%     out.TST              = sum(out.hypnogram>0)*30;
%     out.TSTN1            = sum(out.hypnogram>1)*30;
%%  Events
%   Load events
%     fileID = fopen(EVTpath);
%     if endsWith(EVTpath,'SCO')
%         SCO = textscan(fileID, '%f%f%f%f%s%s%s%f%f', 'Delimiter', '\t', 'EmptyValue', nan);
%     elseif endsWith(EVTpath,'EVTS')
%         textlines = regexp( fileread(EVTpath), '\r?\n', 'split');
% %         EVT = textscan(fileID, '%f%f%s%s%s%s', 'Delimiter', ',', 'EmptyValue', nan, ...
% %         'Headerlines', 2);
%         fs      = 256;
%         lof     = find(contains(lower(textlines),{'lights off'}));
%         lon     = find(contains(lower(textlines),{'lights on'}));
%         lofT    = strsplit(textlines{lof},',');
%         lofT    = strsplit(lofT{3},':');
%         T0(1,4:6) = [str2double(lofT{1}),str2double(lofT{2}), ...
%             roundn(str2double(lofT{3}),1)];
%         times                   = datevec(datetime(T0) ...
%             + seconds((0:size(out.hypnogram,1)-1)*30));
%         out.times            = times;
%         events.lightsofftime    = times(1,:);
%         events.lightsontime     = times(end,:);
%         evts = EVT{1,5}(lof+1:lon-1);
%         stgs = ~contains(EVT{1,6}(lof+1:lon-1),'stage');
%         samp = (EVT{1,1}(lof+1:lon-1)/(fs*30))-lightsOffEpoch;
%         evts = evts(stgs);
%         samp = samp(stgs);
%         EVENTS  = struct();
%         evtCount= zeros(5,1);
%         for e = 1:length(evts)
%             type = evts{e};
%             if contains(type,{'plm','lm','limb'})
%                 evtCount(1) = evtCount(1) + 1;
%                 EVENTS.PLM(evtCount(1),1) = ss;
%                 EVENTS.PLM(evtCount(1),2) = str2double(curEvent.structs.Duration.Text);
%             elseif contains(type,{'arousal'})
%                 evtCount(2) = evtCount(2) + 1;
%                 EVENTS.ARO(evtCount(2),1) = ss;
%                 EVENTS.ARO(evtCount(2),2) = str2double(curEvent.structs.Duration.Text);
%             elseif contains(type,{'hypopnea'})
%                 evtCount(3) = evtCount(3) + 1;
%                 EVENTS.HYP(evtCount(3),1) = ss;
%                 EVENTS.HYP(evtCount(3),2) = str2double(curEvent.structs.Duration.Text);
%             elseif contains(type,{'apnea'})
%                 evtCount(4) = evtCount(4) + 1;
%                 EVENTS.APN(evtCount(4),1) = ss;
%                 EVENTS.APN(evtCount(4),2) = str2double(curEvent.structs.Duration.Text);
%             elseif contains(type,{'spo2', 'desaturation'})
%                 evtCount(5) = evtCount(5) + 1;
%                 EVENTS.SP2(evtCount(5),1) = ss;
%                 EVENTS.SP2(evtCount(5),2) = str2double(curEvent.structs.Duration.Text);
%             end
%         end
%     end
%     fclose(fileID);
    
%     sleepStages         = {'WA','N1','N2','N3','','RE','','WA'};
%     includedEventRows   = find(EVT{1,1}((EVT{1,1} >= lightsOffEpoch).*(EVT{1,1} <= lightsOnEpoch) == 1));
%     includedEventEpochs = EVT{1,1}((EVT{1,1} >= lightsOffEpoch).*(EVT{1,1} <= lightsOnEpoch) == 1);
%     idx = 1;
%     for e = 1:length(includedEventRows)
%         curEvtRow = includedEventRows(e);
%         sleepStage= sleepStages(hypnogram(includedEventEpochs(curEvtRow))+1);
%         curEvtName= lower(strrep(strcat(SCO{1,5}(curEvtRow),'_',sleepStage),' ',''));
%         idxx = 1;
%         for i = 3:size(SCO,2)
%             if i == 5 || i == 7 || SCO{1,i}(e) == (-1)
%                 continue; 
%             else 
%                 events.(curEvtName{:})(idx,idxx) = SCO{1,i}(e);
%                 idxx = idxx + 1;
%             end
%         end
%         idx = idx + 1;
%     end
%     allFields = fields(events);
%         for i = 1:length(allFields)
%             if startsWith(allFields{i},'e_')
%                 data = events.(allFields{i});
%                 data( ~any(data,2), : ) = [];  
%                 events.(allFields{i}) = data;
%                 events.(strcat(allFields{i},'_n')) = size(data,1);
%                 if size(data,1) > 2
%                     q       = quantile(data,[0.25,0.5,0.75]);
%                     if size(data,2) < 3, q = reshape(q,length(q),[]); end
%                     for ii = 1:size(q,2)
%                         events.(strcat(allFields{i},'_q1'))(ii) = q(1,ii);
%                         events.(strcat(allFields{i},'_q2'))(ii) = q(2,ii);
%                         events.(strcat(allFields{i},'_q3'))(ii) = q(3,ii);
%                     end
%                 end
%             end
%         end  
end


















