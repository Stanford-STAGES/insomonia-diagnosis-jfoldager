function out = estimate_wavelet(psg,fs,timess,shouldPlot)
%% Inputs
% eeg               : EEG signal
% fs                : Sampling frequency of EEG signal
%% Outputs
% out.mdl           : nWindows x mdlOrder matrix of AR model parameters
% out.ma_mdl        : out.mdl of moving averaged AR model parameters
% out.ma_mdl_pca    : pca of out.mdl of moving averaged AR model parameters
% out.modHyp        : nWindows x 1 hypnogram corresponding to ss in out.mdl
% out.modP          : nWindows x 1 integers corresponding to P in out.mdl
% out.modC          : nWindows x 1 integers corresponding to C in out.mdl
%% Function
% fprintf('Begin estimate_wavelet...\n')
[wt,f]      = cwt(psg.C3Ref(:,1),fs);
wavelet     = zeros(numel(wt),size(psg.C3Ref,2));
for i = 1:size(psg.C3Ref,2)
    [wt,~] = cwt(psg.C3Ref(:,i),fs);
    wavelet(:,i) = abs(wt(:)).^2;
    out.modC(i)      = find((i>=timess.Cmodified(:,1)).*(i<=timess.Cmodified(:,2)));
    out.modP(i)      = find((i>=timess.Pmodified(:,1)).*(i<=timess.Pmodified(:,2)));
    if shouldPlot 
        subplot(2,1,1)
        plot(psg.C3Ref(:,i));
        subplot(2,1,2)
        imagesc(t,f,log10(abs(wt).^2)), axis xy, ylim([0 20])
        cwt(psg.C3Ref(:,i),fs);
        title(strcat('Clock: ',num2str(events.times(psg.Org_LightsOffEpoch+i-1,4:6)),' SS: ',...
            num2str(psg.Hypnogram(psg.Org_LightsOffEpoch+i-1))));
        pause(1)
    end
%     fprintf('estimate_wavelet: %0.f percent done\n', 100*i/size(psg.C3Ref,2))
end   
% fprintf('estimate_wavelet: Beginning PCA.....\n')
[coe,sco,lat]= pca(zscore(wavelet(1:2:end,:)'));
out.pca.coe  = coe;
out.pca.sco  = sco;
out.pca.lat  = lat;
out.wavelet  = wavelet;
% fprintf('estimate_wavelet: Finished!\n')
end