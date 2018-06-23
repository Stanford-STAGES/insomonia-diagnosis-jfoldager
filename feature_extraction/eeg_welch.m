function out = eeg_welch(psg,fs,timess,freqs,segLength,winLength)
%     hyp         = repmat(psg.Hypnogram,30/segLength,1);  hyp = hyp(:);
%     for i = 1:length(psg.Hypnogram)
%         modC(i) = find((i>=timess.Cmodified(:,1)).*(i<=timess.Cmodified(:,2)));
%         modP(i) = find((i>=timess.Pmodified(:,1)).*(i<=timess.Pmodified(:,2)));
%     end
%     modC         = repmat(modC,30/segLength,1);  modC = modC(:);
%     modP         = repmat(modP,30/segLength,1);  modP = modP(:);
    %%
    sig         = psg.C3Ref(:);     
    sig         = reshape(sig(1:end-mod(length(sig) ...
                    ,fs*segLength)),fs*segLength,[]);
    win         = hamming(winLength*fs);
    [pxx,f]     = pwelch(sig,win,[],[],fs);
    %% 
    fieldss         = fields(freqs);
    out.One2TwentyEight = pwelch(sig,win,[],1:28,fs)';
    out.ALL         = [];
    totalPow        = bandpower(pxx,f, freqs.TOTAL,'psd');
    for i = 1:length(fieldss)-1
        out.(fieldss{i}) = (bandpower(pxx,f,freqs.(fieldss{i}),'psd')./totalPow)';
        out.ALL = [out.ALL,out.(fieldss{i})];
    end
    out.ALL = [out.ALL,out.One2TwentyEight,hyp,modC,modP];
    out.modH = hyp;
    out.modC = modC;
    out.modP = modP;
end