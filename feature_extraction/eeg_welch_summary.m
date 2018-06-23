function out = eeg_welch_summary(SPECTRAL,qs)
    field  = fields(SPECTRAL);
    nFreqs = length(fields(SPECTRAL));
    sss = [0 1 2 3 5];
    out = nan(1,5*(nFreqs-5)*length(qs));
    idx = 1;
    for i = 3:(nFreqs-3)
        for ii = 1:5
            out(idx:idx+length(qs)-1) = quantile(SPECTRAL.(field{i})(SPECTRAL.modH == sss(ii)),qs);
            idx = idx + length(qs);
        end
    end
end