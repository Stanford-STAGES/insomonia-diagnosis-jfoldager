function features = coherence(sig1,sig2,fs,winLength,epochLength)
%   Calculates 
%   amplitude, mean, median, mode, standard deviation, variance, skewness, kurtosis   
    sig1     = reshape(sig1(1:end-mod(length(sig1) ...
                ,fs*epochLength)),fs*epochLength,[]);     
    sig2     = reshape(sig2(1:end-mod(length(sig2) ...
                ,fs*epochLength)),fs*epochLength,[]);  
    features.coherence = mscohere(sig1,sig2,hamming(winLength*fs),[],[],fs);
end