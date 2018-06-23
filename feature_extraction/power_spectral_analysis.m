function features = power_spectral_analysis(sig,fs,windowLength,epochLength)
    sig = reshape(sig(1:end-mod(length(sig) ...
                ,fs*epochLength)),fs*epochLength,[]);
    win = hamming(windowLength*fs);
    [pxx,f] = pwelch(sig,win,[],[],fs);
    %% Frequency interval definition
    freqs.slow_waves = [0.5 1];
    freqs.delta = [1 4];
    freqs.theta = [4 8];
    freqs.alpha = [8 12];
    freqs.sigma = [12 15];
    freqs.fast_sigma = [13.5 15];
    freqs.slow_sigma = [12 13.5];
    freqs.beta = [15 20];
    allBands = fieldnames(freqs);
    totalPower = bandpower(pxx,f,[0.5 20],'psd');
    for i = 1:length(allBands)
        features.(allBands{i}) = bandpower(pxx,f,freqs.(allBands{i}),'psd');
        features.(strcat(allBands{i},'_vs_total')) =  features.(allBands{i})./totalPower;
    end
    for i = 1:length(allBands)
        for ii = 1:length(allBands)
            if ii ~= i
                features.(strcat(allBands{i},'_vs_',allBands{ii})) = ...
                    features.(allBands{i})./features.(allBands{ii});
            end
        end
    end
    [features.power_1_to_28_hz,~] = pwelch(sig,win,[],1:28,fs);
end

