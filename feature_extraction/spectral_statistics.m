function features = spectral_statistics(sig,fs,segLength,events,takeavg)
    sig = reshape(sig(1:end-mod(length(sig) ...
                ,fs*segLength)),fs*segLength,[]);
    sig          = sig(:,events.LOF:end);
    win = hamming(1*fs);
    [pxx,f] = pwelch(sig,win,[],2*fs,fs);
    %% Frequency interval definition
    freqs.SW = [0.5 1];
    freqs.DE = [1 4];
    freqs.TH = [4 8];
    freqs.AL = [8 12];
    freqs.SS = [12 13.5];
    freqs.FS = [13.5 15];
    freqs.SB = [15 20];
    freqs.FB = [20 30];
    allBands = fieldnames(freqs);
    totalPower = bandpower(pxx,f,[0.5 35],'psd');
    features = table();
    for i = 1:length(allBands)
        features.(allBands{i}) = (bandpower(pxx,f,freqs.(allBands{i}),'psd')./totalPower)';
    end
    if takeavg
    features = calculateStatistics(features,events.hypnogram,[]);
    end
end

