function out = filtering(sig,fs,type)
% function filters signal based on type (either eeg or ecg)
% uses 8'th butterworth filter with Fc1 = 0.1 and Fc2 = 35
% input: eeg/ecg signal
% output: filtered input signal
    if isequal(type,'eeg')
        Fc1     = 0.3;                                              % First Cutoff Frequency
        Fc2     = 35;                                               % Second Cutoff Frequency
        [b,a]   = butter(8/2,[Fc1/(fs/2),Fc2/(fs/2)],'bandpass'); 
        out     = filtfilt(b,a,sig);                                % zero-phase filtering
    end
    if isequal(type,'ecg')
        Fc1     = 5;                                             % pan1985real
        Fc2     = 15;                                            % pan1985real
        [b,a]   = butter(8/2,[Fc1/(fs/2),Fc2/(fs/2)],'bandpass'); 
        out     = filtfilt(b,a,sig);                              
    end
end

