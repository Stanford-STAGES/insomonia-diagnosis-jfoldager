function out = sample_entropy(sig,L,dis)
% Generally we take the value of L to be 2 and the value of 
% r to be  0.2 times std
% A smaller value of {\displaystyle SampEn} SampEn also indicates 
% more self-similarity in data set or less noise.
    sig = sig(:);           % Stretch
    N = length(sig);        % Number of samples
    SampEn = zeros(N-L,1);    
    for n=1:N
        if n+L <= N
            X_m         = sig(n:(L+(n-1)));
            X_m_plus_1  = sig(n:(L+n));
            A           = sum(sum(tril(pdist2(X_m,X_m,'chebychev')<dis,-1)));
            B           = sum(sum(tril(pdist2(X_m_plus_1,X_m_plus_1,'chebychev')<dis,-1)));
            SampEn(n)   = -log(A/B);  
        end
    end  
    out = SampEn;
end