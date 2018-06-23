function out = detect_qrs(ecg,fs,epochLength,plotFlag)
    ECGFilt = filtering(ecg,fs,'ecg'); % 5-15 Hz
    ECG = ECGFilt.^4;
    ECG = reshape(ECG(1:end-mod(length(ECG) ...
                ,fs*epochLength)),fs*epochLength,[]);
    maxPulse = 200; minQQSam = floor(fs/(maxPulse/60));
    nClusters = 6;
    C = zeros(nClusters,size(ECG,2));
    lks = zeros((length(ECG(:))/fs)*3,1);
    idx = 1;
    for i = 1:size(ECG,2) % Windowing
        [idxI,C(:,i)] = kmeans(ECG(:,i),nClusters);
        [C(:,i),ii] = sort(C(:,i));
        m = ii(end-3);
        q = min(ECG(sum(idxI == m),i));
        [~,lk] = findpeaks(ECG(:,i),1:length(ECG(:,i)),'MinPeakDistance',minQQSam...
        ,'MinPeakHeight',q);
        lks(idx:idx+length(lk)-1) = lk+(i-1)*fs*epochLength;
        idx = length(lk) + 1;
        fprintf('%3f\n',i/size(ECG,2));
    end
    lks = lks(lks>0);
    out.lks = lks;
    if plotFlag
        figure
        hold on
        plot(lks,abs(ECGFilt(lks)),'o')
        plot(abs(ECGFilt))
    end
end
 







