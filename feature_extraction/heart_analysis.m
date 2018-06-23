function out = heart_analysis(ecg,fs,segmentLength,events)%,hyp,plotFlag,artefacts,epochLength,slidePath,timess,segmentLength)
% Raw, unfiltered ECG signal (1D)
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ANALYSIS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[~,qrsIdx,ECG,~]    = pan_tompkin(ecg,fs,false,segmentLength);
winHalf = 10;
minDist = floor(0.2*fs);
for i = 1:length(qrsIdx)
    [~,idx] = max(abs(ecg(qrsIdx(i)-winHalf:qrsIdx(i)+winHalf)));
    qrsIdx(i) = qrsIdx(i) + idx - winHalf - 1;
    if i > i && abs(qrsIdx(i-1)-qrsIdx(i)) <= minDist
        [~,idx] = max(abs(ecg(qrsIdx(i-1:i))));
        qrsIdx(i+idx-1) = [];
    end
end
% artFree             = artefacts==0;
% out.ECG_Reshaped    = ECG;
% out.ArtefactFreeECG = artFree;
locations           = zeros(size(ECG));
locations(qrsIdx)   = qrsIdx;
qrsIdxForAnalysis   = locations;%(:,artefacts==0);
% qrsIdxForAnalysis   = locations(:,artFree);
poinCare           = zeros(length(qrsIdxForAnalysis(qrsIdxForAnalysis>0))-1,2);
grp                 = zeros(length(qrsIdxForAnalysis(qrsIdxForAnalysis>0))-1,1);
artefactFree        = 1:size(qrsIdxForAnalysis,2);%find(artefacts == 0);
% segmentLength       = 5; segmentLength = epochLength/segmentLength;
% delay = 0;
% numComponents   = 1;
% if plotFlag
%     tECG= linspace(0,length(ecg),length(ecg))/fs;
%     figure
%         hold on
%         plot(tECG,abs(ecg))
%         plot(tECG(qrsIdx+delay),abs(ecg(qrsIdx+delay)),'o')
%         hold off
%         title('QRS detection')
%         xlabel('Time (sec)')
%         ylabel('Abs. mV')
%         print('-depsc2', '-loose', strcat(slidePath,'QRS.eps'));
%     figure
%         subplot(2,1,1)
%             hold on
%             plot(tECG,abs(ecg))
%             plot(tECG(qrsIdx+delay),abs(ecg(qrsIdx+delay)),'o')
%             hold off
%             title('QRS detection (raw ECG signal)')
%             ylabel('Abs. V')
%         subplot(2,1,2)
%             tECG= linspace(0,length(ecg),length(ecg))/fs;
%             tECGRes = reshape(tECG(1:end-mod(length(tECG) ...
%             ,fs*epochLength)),fs*epochLength,[]); tECGRes = tECGRes(:,artFree); tECGRes = tECGRes(:);
%             ecgRes = reshape(ecg(1:end-mod(length(ecg) ...
%             ,fs*epochLength)),fs*epochLength,[]); ecgRes = ecgRes(:,artFree); ecgRes = ecgRes(:);
%             hold on
%             newQRSIdx = qrsIdxForAnalysis(qrsIdxForAnalysis>0);
%             plot(tECGRes,abs(ecgRes))
%             plot(tECG(newQRSIdx+delay),abs(ecg(newQRSIdx+delay)),'o')
%             hold off
%             title('QRS detection (EEG and ECG artefact epochs excluded)')
%             xlabel('Time (sec)')
%             ylabel('Abs. V')
%         subplot(3,1,3)
%             tECGRes = reshape(tECG(1:end-mod(length(tECG) ...
%             ,fs*epochLength)),fs*epochLength,[]); tECGRes = tECGRes(:,artFree); tECGRes = tECGRes(:);
%             ecgRes = reshape(ecg(1:end-mod(length(ecg) ...
%             ,fs*epochLength)),fs*epochLength,[]); ecgRes = ecgRes(:,artFree); ecgRes = ecgRes(:);
%             hold on
%             newQRSIdx = qrsIdxForAnalysis2(qrsIdxForAnalysis2>0);
%             plot(tECGRes,abs(ecgRes))
%             plot(tECG(newQRSIdx+delay),abs(ecg(newQRSIdx+delay)),'o')
%             hold off
%             title('QRS detection (ECG artefact epochs excluded)')
%             xlabel('Time (sec)')
%             ylabel('Abs. V')
% end

% idx = 1;
% i = 1;
% Standard deviation of the average NN intervals calculated over 5-minute periods
% while i <= length(artefactFree)-segmentLength+1
%     if sum(diff(artefactFree(i:i+segmentLength-1))) == segmentLength - 1
%         q = locations(:,artefactFree(i):artefactFree(i)+segmentLength-1);
%         q = diff(q(q>0)./fs);
%         out.ANN5min(idx) = nanmean(q);        
%         idx = idx + 1;
%         i = i + segmentLength;
%     else 
%         i = i + 1;
%     end
% end
% out.SDANN5min =  nanstd(out.ANN5min);
idxPC= 1;
secSplit    = segmentLength; % if 30 then defualt
cur         = qrsIdxForAnalysis(:,1); 
curSecSplitsec     = reshape(cur(1:end-mod(length(cur) ...
            ,fs*secSplit)),fs*secSplit,[]);
out = table();
out.MNN     = nan(size(qrsIdxForAnalysis,2),size(curSecSplitsec,2));
out.MPULS   = nan(size(qrsIdxForAnalysis,2),size(curSecSplitsec,2));
out.RMSSD   = nan(size(qrsIdxForAnalysis,2),size(curSecSplitsec,2));
out.NN50    = nan(size(qrsIdxForAnalysis,2),size(curSecSplitsec,2));
out.PNN50   = nan(size(qrsIdxForAnalysis,2),size(curSecSplitsec,2));
out.SD2C    = nan(size(qrsIdxForAnalysis,2),size(curSecSplitsec,2));
out.SD1C    = nan(size(qrsIdxForAnalysis,2),size(curSecSplitsec,2));
out.PC_ang  = nan(size(qrsIdxForAnalysis,2),size(curSecSplitsec,2));
for ii = 1:size(qrsIdxForAnalysis,2)
    cur         = qrsIdxForAnalysis(:,ii); 
    curSecSplitsec     = reshape(cur(1:end-mod(length(cur) ...
                ,fs*secSplit)),fs*secSplit,[]);
    for i = 1:size(curSecSplitsec,2)
        curcur      = curSecSplitsec(:,i); 
        qrsIdxs     = curcur(curcur>0);
        qrsTime     = qrsIdxs./fs;
        qrsTimeDiff = diff(qrsTime);
        if ~isempty(qrsTimeDiff) && length(qrsTimeDiff) > 2
            instPulse           = 60./qrsTimeDiff;
            % Mean of the NN intervals
            out.MNN(ii,i)       = mean(qrsTimeDiff);  
            % Mean of the instant pulse
            out.MPULS(ii,i)       = mean(instPulse);  
            % Square root of the mean squared differences of successive NN intervals
            out.RMSSD(ii,i)     = rms(diff(qrsTimeDiff));  
            % Number of interval differences of successive NN intervals greater than 50 ms
            out.NN50(ii,i)      = sum(abs(diff(qrsTimeDiff))>0.05); 
            % Dividing NN50 by the total number of NN intervals
            out.PNN50(ii,i)     = out.NN50(ii,i)/length(qrsTimeDiff);  
            [a,b,phi]           = finding_length([qrsTimeDiff(1:end-1),qrsTimeDiff(2:end)]);
            out.SD2C(ii,i)      = a;
            out.SD1C(ii,i)      = b;
            out.PC_ang(ii,i)    = phi;
            try
                hypnogram(ii,i)      = events.hypnogram(ii);
            catch
            end
%             out.modC(ii,i)      = find((i>=timess.Cmodified(:,1)).*(i<=timess.Cmodified(:,2)));
%             out.modP(ii,i)      = find((i>=timess.Pmodified(:,1)).*(i<=timess.Pmodified(:,2)));
        end
    end
    qrsIdxs     = cur(cur>0);
    qrsTime     = qrsIdxs./fs;
    qrsTimeDiff = diff(qrsTime);
    % Poincare 
    if ~isempty(qrsTimeDiff)
        poinCare(idxPC:idxPC+length(qrsTimeDiff(1:end-1))-1,1)  = qrsTimeDiff(1:end-1);
        grp(idxPC:idxPC+length(qrsTimeDiff(1:end-1))-1,1)       = ones(size(qrsTimeDiff(1:end-1))).*i;
        poinCare(idxPC:idxPC+length(qrsTimeDiff(1:end-1))-1,2)  = qrsTimeDiff(2:end);
        idxPC           = idxPC+length(qrsTimeDiff(1:end-1)); 
    end
end
%     out.MNN         = out.MNN(:);
%     out.MPu         = out.MPu(:);
%     out.RMSSD       = out.RMSSD(:);
%     out.NN50        = out.NN50(:);
%     out.pNN50       = out.pNN50(:);
%     out.SD2C        = out.SD2C(:);
%     out.SD1C        = out.SD1C(:);
%     out.PC_ang      = out.PC_ang(:);
    
    out = calculateStatistics(out,hypnogram,[]);
    
    
%     out.HRV         = [out.MNN, ...
%                        out.RMSSD, ...
%                        out.NN50, ...
%                        out.pNN50, ...
%                        out.SD2C, ...
%                        out.SD1C, ...
%                        out.PC_ang, ...
%                         ];
    %% Assign sleep stage to each datapoint
%     grp             = grp(grp>0);
%     sleepPoinCare   = ones(size(grp))*(-1);
%     for i = 1:length(grp), sleepPoinCare(i) = hyp(grp(i)); end
%     poinCare       = [poinCare(poinCare(:,1)>0,1),poinCare(poinCare(:,2)>0,2)];
%     out.poinCare   = poinCare;    
%     out.sleepPoinCare   = sleepPoinCare;    
    %% Fit gaussian mixture model
%     gm              = fitgmdist(poinCare,numComponents);
%     y               = pdf(gm,poinCare);
%     poinCare        = poinCare(y>0.05,:);
%     grp             = grp(y>0.05);
%     sleepPoinCare   = sleepPoinCare(y>0.05);
%     if plotFlag 
%         figure
%             s = scatter(poinCare(:,1),poinCare(:,2),1,grp(grp>0));
%             phi = colorbar; colormap(jet(length(mean(hyp,1))))
%             s.LineWidth = 2;
%             title('Poincare plot')
%             xlabel('Time (sec)')
%             ylabel('Time (sec)')
%             xlim([0 2])
%             ylim([0 2])
%             
%         figure
%              gscatter(poinCare(:,1),poinCare(:,2),sleepPoinCare,'brygkm','xoxoxo');
%              
%         figure
%             xlimitttt = [0.5 1.5];
%             ylimitttt = xlimitttt;
%             subplot(2,3,1)
%             crit = sleepPoinCare == 0;
%             s = scatter(poinCare(crit,1),poinCare(crit,2),1,grp(crit)); s.LineWidth = 2;
%             xlim(xlimitttt)
%             ylim(ylimitttt)
%             colormap(jet(50))
%             phi = colorbar; %c.Ticks = 0:floor(max(grp)/4):max(grp);
%             title('Awake')
%             subplot(2,3,2)
%             crit = sleepPoinCare == 1;
%             s = scatter(poinCare(crit,1),poinCare(crit,2),1,grp(crit)); s.LineWidth = 2;
%             xlim(xlimitttt)
%             ylim(ylimitttt)
%             colormap(jet(50))
%             phi = colorbar; %c.Ticks = 0:floor(max(grp)/4):max(grp);
%             title('N1')
%             subplot(2,3,3)
%             crit = sleepPoinCare == 2;
%             s = scatter(poinCare(crit,1),poinCare(crit,2),1,grp(crit)); s.LineWidth = 2;
%             xlim(xlimitttt)
%             ylim(ylimitttt)
%             colormap(jet(50))
%             phi = colorbar; %c.Ticks = 0:floor(max(grp)/4):max(grp);
%             title('N2')
%             ylabel('Time (sec)')
%             subplot(2,3,4)
%             crit = ((sleepPoinCare == 3) + (sleepPoinCare == 4)) > 0;
%             s = scatter(poinCare(crit,1),poinCare(crit,2),1,grp(crit)); s.LineWidth = 2;
%             xlim(xlimitttt)
%             ylim(ylimitttt)
%             colormap(jet(50))
%             phi = colorbar; %c.Ticks = 0:floor(max(grp)/4):max(grp);
%             title('N3')
%             subplot(2,3,5)
%             crit = sleepPoinCare == 5;
%             s = scatter(poinCare(crit,1),poinCare(crit,2),1,grp(crit)); s.LineWidth = 2;
%             xlim(xlimitttt)
%             ylim(ylimitttt)
%             colormap(jet(50))
%             phi = colorbar; %c.Ticks = 0:floor(max(grp)/4):max(grp);
%             title('REM')
%             xlabel('Time (sec)')
%              
%         figure;
%             xlimitttt = [0.5 1.5];
%             ylimitttt = xlimitttt;
%             subplot(2,2,1)
%                 hold on
%                 plot(poinCare(:,1),poinCare(:,2),'*')
%                 gm = fitgmdist(poinCare,numComponents);
%                 gmPDF = @(x,y)pdf(gm,[x y]);
%                 h = fcontour(gmPDF);
%                 xlim(xlimitttt)
%                 ylim(ylimitttt)
%                 ylabel('Time (sec)')
%                 title('Poincare Plot')
%                 axis image
%                 hold off
%             subplot(2,2,3)
%                 hold on
%                 y = pdf(gm,poinCare);
%                 filterGM = poinCare(y>0.05,:);
%                 plot(filterGM(:,1),filterGM(:,2),'*')
%                 gm = fitgmdist(filterGM,numComponents);
%                 gmPDF = @(x,y)pdf(gm,[x y]);
%                 h = fcontour(gmPDF);
%                 xlim(xlimitttt)
%                 ylim(ylimitttt)
%                 ylabel('Time (sec)')
%                 xlabel('Time (sec)')
%                 axis image
%                 hold off
%             subplot(2,2,[2 4])
%                 histogram2(filterGM(:,1),filterGM(:,2),100)
%                 ylabel('Time (sec)')
%                 xlabel('Time (sec)')
%                 
%         figure
%             subplot(2,1,2)
%             hold on 
%             plot(out.MNN./max(out.MNN))
%             plot(out.PC_ang./max(out.PC_ang))
%             plot(out.SD2C./max(out.SD2C))
%             plot(out.SD1C./max(out.SD1C))
%             legend({'mu','angle', 'a', 'b'})
%             hold off
%             subplot(2,1,1)
%             plot(hyp)
%     end     
end