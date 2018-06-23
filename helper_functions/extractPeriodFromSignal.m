function signal = extractPeriodFromSignal(signal,fs,epochLength,startEpoch,stopEpoch)
signal                         = reshape(signal(1:end-mod(length(signal) ...
                                ,fs*epochLength)),fs*epochLength,[]); 
signal                         = signal(:,startEpoch:stopEpoch); 
signal                         = signal(:);
end