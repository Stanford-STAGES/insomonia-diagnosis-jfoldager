function out = brain_rate(sig,fs,winLength,epochLength)
%   Brain rate represents the mean EEG frequency weighted over the brain potential (or power) distribution spectrum
%%  Initial
    freqs.slow_waves = [0.5 1];
    freqs.delta = [1 4];
    freqs.theta = [4 8];
    freqs.alpha = [8 12];
    freqs.slow_sigma = [12 13.5];
    freqs.sigma = [12 15];
    freqs.fast_sigma = [13.5 15];
    freqs.beta = [15 20];
    allBands = fieldnames(freqs);
    f_i = zeros(size(allBands));
    f_b = zeros(size(sig,2),1);
    %%
    for i   = 1:length(allBands),f_i(i) = mean(freqs.(allBands{i})); end
    sig     = reshape(sig(1:end-mod(length(sig) ...
                ,fs*epochLength)),fs*epochLength,[]);
    win     = hamming(winLength*fs);
    [pxx,~] = pwelch(sig,win,[],f_i,fs);    
    for i   = 1:size(sig,2),f_b(i) = sum(pxx(:,i).*f_i)/sum(pxx(:,i)); end
    out = f_b;
end