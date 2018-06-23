function out = brain_symmetry_index( R,L,fs,winLength,epochLength )
    R           = reshape(R(1:end-mod(length(R) ...
                ,fs*epochLength)),fs*epochLength,[]);
    L           = reshape(L(1:end-mod(length(L) ...
                ,fs*epochLength)),fs*epochLength,[]);
    win         = hamming(winLength*fs);
    [pxxR,~]    = pwelch(R,win,[],[],fs);    
    [pxxL,~]    = pwelch(L,win,[],[],fs);
    M           = size(pxxR,1);
    bsi         = sum(abs((pxxR-pxxL)./(pxxR+pxxL)))./M;
    out         = bsi;
end

