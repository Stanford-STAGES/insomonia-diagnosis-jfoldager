function [LFinal,Fn] = dfa(sig,lMax,shouldPlot)
%DFA Summary of this function goes here
%   Detailed explanation goes here
idx = 1;
for L = 1:lMax
    if mod(length(Xt),L) == 0
        Xt          = cumsum(sig - mean(sig));
        XtRe        = reshape(Xt,L,[]);
        Fn(idx)     = sqrt(mean(detrend(XtRe(:)).^2));
        LFinal(idx) = L;
        idx = idx + 1;
    end
end
if shouldPlot
    plot(log(LFinal),log(Fn))
    xlabel('log(n)'),ylabel('log(F(n))')
end
end