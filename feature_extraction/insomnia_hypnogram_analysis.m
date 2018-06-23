function out = insomnia_hypnogram_analysis(hyp,epochLength)
        sSOP            = find((hyp>0).*(hyp<7) == 1); sSOP=sSOP(1)-1;                    % Sleep onset period sample
        sSOP_n1         = find((hyp>1).*(hyp<7) == 1); sSOP_n1=sSOP_n1(1)-1;              % Sleep onset period sample
        sLastSleep      = find(hyp ~= 0); sLastSleep = sLastSleep(end); % Last sleep epoch
        sleepHyp        = hyp(1:sLastSleep);                            % Only hypnogram for sleep period
        %% Total sleep time, sleep onset latency and sleep efficiency
        out.h_TST       = (sum((hyp > 0).*(hyp<7) == 1))*epochLength;           % sec  
        out.h_SOL       = (sSOP)*epochLength;                                   % sec
        out.h_SOL_n1    = (sSOP_n1)*epochLength;                                % sec
        out.h_SE        = sum((sleepHyp > 0).*(sleepHyp<7) == 1)/length(sleepHyp);   % percentage
        out.h_SE_n1     = sum((sleepHyp > 1).*(sleepHyp<7) == 1)/length(sleepHyp);   % percentage       
        %% Awake above 2.5 and 5
        a = diff(find(diff(diff(cumsum((sleepHyp == 0) + (sleepHyp == 7)))) ~= 0));
        a = a(2:2:end);
        out.h_awake_2_5 = sum(a>=(2.5*60)/epochLength);
        out.h_awake_5   = sum(a>=(5*60)/epochLength);       
        %% Awake or N1 above 2.5 and 5
        a = diff(find(diff(diff(cumsum((sleepHyp == 0) + (sleepHyp == 1) + (sleepHyp == 7)))) ~= 0));
        a = a(2:2:end);
        out.h_awake_n1_2_5  = sum(a>=(2.5*60)/epochLength);
        out.h_awake_n1_5    = sum(a>=(5*60)/epochLength);          
end