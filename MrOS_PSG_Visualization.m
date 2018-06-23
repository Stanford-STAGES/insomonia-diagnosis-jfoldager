%% Load and sort data
% This scripts is a generic approach to load the data properly.
% 
% 
% Author: Jonathan Foldager
%% Initial/startup code
startup                                                 % 
load('insomniaseverity.mat');
listing     = extractfield(dir('Data/MrOS/polysomnography/'), 'name')';
curSubj     = 1;
subjs       = listing(contains(listing,'.edf'));
subjID      = split(subjs(curSubj),'.'); subjID = subjID{1};
% EDF
if exist(strcat(subjID,'.EDF'), 'file') == 2
    filename    = strcat(subjID,'.EDF'); [data, header] = readEDF(filename);
end
%% Sampling frequency and time
fs      = header.samplerate(contains(header.labels,'C3'));
T       = 1/fs;
L       = length(data{contains(header.labels,'C3')});
tSec       = linspace(0,L,L)*T;
tMin       = tSec/60;
tHou       = tMin/60;
%% EEG
C3M2    = data{contains(header.labels,'C3')} - data{contains(header.labels,'M2')};
C4M1    = data{contains(header.labels,'C4')} - data{contains(header.labels,'M1')};
%% Filter(s)
% d = designfilt('bandpassfir', ...       % Response type
%        'StopbandFrequency1',0.1, ...    % Frequency constraints
%        'PassbandFrequency1',0.6, ...
%        'PassbandFrequency2',35, ...
%        'StopbandFrequency2',50, ...
%        'SampleRate',fs);               % Sample rate
% fvtool(d)
% C3M2Filt= filter(d,C3M2);
% C4M1Filt= filter(d,C4M1);
Fc1 = 0.3; % First Cutoff Frequency
Fc2 = 35; % Second Cutoff Frequency
[b,a] = butter(8/2,[Fc1/(fs/2),Fc2/(fs/2)],'bandpass');
C3M2 = filtfilt(b,a,C3M2); % zero-phase filtering
C4M1 = filtfilt(b,a,C4M1); % zero-phase filtering
%% Hypnogram
DOMnode = xml2struct(strcat(subjID,'-nsrr.xml'));
stru    = DOMnode.PSGAnnotation.ScoredEvents.ScoredEvent;
stru    = cell2struct(stru,{'structs'});
ss      = zeros(length(stru),1);
ssStart = zeros(length(stru),1);
ssDur   = zeros(length(stru),1);
hyp     = [];
for i = 1:length(stru)
    if contains(stru(i).structs.EventType.Text,'Stage')
        str = split(stru(i).structs.EventConcept.Text,'|');
        ss(i) = str2double(str(2));
        str = stru(i).structs.Start.Text;
        ssStart(i) = str2double(str);
        str = stru(i).structs.Duration.Text;
        ssDur(i) = str2double(str);
        hyp = [hyp; repmat(ss(i),fs*ssDur(i),1)];
    end
end
if length(hyp) < length(C3M2), hyp = [hyp; zeros(length(C3M2)-length(hyp),1)]; end
%% PSG example
figure('units','normalized','outerposition',[0.1 0.05 0.6 0.9])
hold on
newFs = 32;
y = decimate(hyp,fs/newFs,1);
plot(y./max(y),'black');
for i = 1:length(header.labels)
    sig = data{i};
    if contains(header.labels{i},{'Position','STAT','Airflow','Position', 'CannulaFlow'})
        if header.samplerate(i) > newFs
            y = decimate(sig,header.samplerate(i)/newFs,1);
        else
            y = resample(sig,newFs,header.samplerate(i),0);
        end        
    else
        if header.samplerate(i) > newFs
            y = decimate(sig,header.samplerate(i)/newFs);
        else
            y = resample(sig,newFs,header.samplerate(i));
        end        
    end  
    if length(y) ~= length(data{14})
        sprintf('ERROR!! At i = %i',i)
    else
        plot((y./abs(max(y)))-i*3,'black');
        sprintf("%s is done.",header.labels{i})
    end
end
hold off
axis tight
yticks(-i*3:3:0)
xticks([length(y)*0.25 length(y)*0.5 length(y)*0.75 length(y)])
set(gca,'YTickLabel',[flip(header.labels); 'Hypnogram']);
set(gca,'XTickLabel',[round(tHou(floor(length(tHou)*0.25)),1)...
    round(tHou(floor(length(tHou)*0.5)),1)...
    round(tHou(floor(length(tHou)*0.75)),1)...
    round(tHou(end),1)]);
xlabel('Time (hrs)')
% print(strcat(plotPath2,'PSG_Example.eps'),'-depsc')
%% Plot - before and after filtering
seconds = [60*60 60*61]; minutes = seconds/60; hours = minutes/60;
figure(1);
    subplot(2,1,1)
        plot(tSec,C3M2)
        xlim(seconds);
    subplot(2,1,2)
        plot(tSec,C3M2)
        xlim(seconds);
winSize = fs/4;
figure(2);
    subplot(2,1,1)
        spectrogram(C3M2,kaiser(winSize),[],[],fs,'yaxis');
        xlim(hours);
    subplot(2,1,2)
        spectrogram(C3M2,kaiser(winSize),[],[],fs,'yaxis');
        xlim(hours);
%% Plot 
seconds = [400 420]; minutes = seconds/60; hours = minutes/60;
winSize = fs/4;
% C3M2
figure(2);
    subplot(3,1,1)
        plot(tHou,hyp)
        yticklabels({'Awake','NREM1','NREM2','NREM3','','REM'})
        title('Hypnogram')
        axis tight
%         xlim(seconds);
    subplot(3,1,2)
        plot(tHou,C3M2)
        ylabel('$\mu V$')
        axis tight
        title('C3M2')
%         xlim(seconds);
    subplot(3,1,3)
        spectrogram(C3M2,kaiser(winSize),[],[],fs,'yaxis');
        colorbar off
        title('Spectrogram')
        ylabel('Hz')
        axis tight
%         xlim(hours);
% C4M1
figure('units','normalized','outerposition',[0.1 0.05 0.6 0.9])
    subplot(3,1,1)
        h = plot(tHou,hyp);
        title('Hypnogram')
        axis tight
        set(gca,'YTickLabel',{'Awake','NREM1','NREM2','NREM3','','REM'})
%         xlim(seconds);
    subplot(3,1,2)
        plot(tHou,C4M1)
        ylabel('$\mu V$')
        axis tight
        title('C4M1')
%         xlim(seconds);
    subplot(3,1,3)
        [F, fF, tF] = spectrogram(C4M1,kaiser(winSize),[],[],fs,'yaxis');
        imagesc(tF/60/60,fF,20*(log10(abs(F))))
        c = colorbar('southoutside');
        c.Label.String = 'Power/Hz (dB/Hz)';
        c.TickLabelInterpreter = 'latex';
        axis xy
        title('Spectrogram')
        ylabel('Hz')
        xlabel('Time (hrs)')
        ylim([0 50]);
%         xlim(hours);
print(strcat(plotPath2,'MrOS_Ins_EEG_Example.eps'),'-depsc')



%%
if false
k = 6;
figure
    subplot(k,1,1)
    hold on
    plot(features.bandpower.slow_waves_vs_total)
    plot(features.bandpower.delta_vs_total)
    plot(features.bandpower.theta_vs_total)
    plot(features.bandpower.alpha_vs_total)
    plot(features.bandpower.beta_vs_total)
%     plot([lightsOffEpo lightsOffEpo],[0 1],'r','LineWidth',3)
    title('Normalized spectral energy')
    legend({'slow waves','delta','theta','alpha','beta'})
    subplot(k,1,2)
    hold on
    plot(mean(PSG.Hypnogram(:,lightsOffEpo:end)),'b')
    set(gca,'Ydir','reverse')
%     plot([lightsOffEpo lightsOffEpo],[min(mean(PSG.Hypnogram)) max(mean(PSG.Hypnogram))],'r','LineWidth',3)
    ylim([min(mean(PSG.Hypnogram)) max(mean(PSG.Hypnogram))])
    title('Hypnogram')
    subplot(k,1,3)
    hold on
    plot(features.br)
%     plot([lightsOffEpo lightsOffEpo],[min(features.br) max(features.br)],'LineWidth',3)
    ylim([min(features.br) max(features.br)])
    title('Brain rate')
    subplot(k,1,4)
    hold on
    plot(features.hp.mob);
%   plot([lightsOffEpo lightsOffEpo],[min(features.hp.mob) max(features.hp.mob)],'LineWidth',3);
    ylim([min(features.hp.mob) max(features.hp.mob)])
    title('Hjorth (mobility)')
    subplot(k,1,5)
    hold on
    plot(features.hp.com)
%     plot([lightsOffEpo lightsOffEpo],[min(features.hp.com) max(features.hp.com)],'LineWidth',3)
    ylim([min(features.hp.com) max(features.hp.com)])
    title('Hjorth (complexity)')
    subplot(k,1,6)
    hold on
    plot(mean(PSG.HR(:,interval)))
%     plot([lightsOffEpo lightsOffEpo],[min(mean(PSG.HR)) max(mean(PSG.HR))],'LineWidth',3)
    ylim([min(mean(PSG.HR(:,interval))) max(mean(PSG.HR(:,interval)))])
    title('Mean heart rate')
    xlabel('Epochs')
end
