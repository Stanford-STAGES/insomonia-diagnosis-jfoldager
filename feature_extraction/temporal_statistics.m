function features = temporal_statistics(sig,fs,segLength,events,takeavg)
%   Calculates 
%   mean, median, mode, standard deviation, variance, skewness, kurtosis ..
%   for sig.
    sig     = reshape(sig(1:end-mod(length(sig) ...
                ,fs*segLength)),fs*segLength,[]);       
    features            = table();
    features.MEAN       = mean(sig)';
    features.MED        = median(sig)';
    features.MODE       = mode(sig)';
    features.STD        = std(sig)';
    features.SKEW       = skewness(sig)';
    features.KURT       = kurtosis(sig)';
    features.ZCR        = sum(diff(sig>0) ~= 0)';
    features.P2P        = min(sig)'-max(sig)';
    features.RMS        = rms(sig)';
    features.DE1MAX     = max(diff(sig))';
    features.DE1MEAN    = mean(diff(sig))';
    features.DE2MAX     = max(diff(sig,2))';
    features.DE2MEAN    = mean(diff(sig,2))';
    
    features            = softnorm(log10(features(events.LOF:end,:)),[0.05 0.95]);
    if takeavg
        features        = calculateStatistics(features,events.hypnogram,[]);
    end
end