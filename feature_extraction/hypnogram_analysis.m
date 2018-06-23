function out = hypnogram_analysis(hyp)
    sLastSleep  = find(hyp ~= 0); sLastSleep = sLastSleep(end); % Last sleep sample
    sleepHyp    = hyp(1:sLastSleep);                            % Only hypnogram for sleep period
    sleepHyp(sleepHyp==4) = 3;
    sleepHyp    = reshape(sleepHyp,1,length(sleepHyp));
%% Features
    out         = table();
    out.TST     = sum((sleepHyp>0))*30/60;
    out.TSTN1   = sum((sleepHyp>1))*30/60;
    out.SOL     = find((sleepHyp>0).*(sleepHyp<7) == 1); out.SOL=(out.SOL(1)-1)*30/60; % Sleep onset period sample
    out.SOLN1   = find((sleepHyp>1).*(sleepHyp<7) == 1); out.SOLN1=(out.SOLN1(1)-1)*30/60;
    out.SE      = sum((sleepHyp > 0).*(sleepHyp<7) == 1)/length(sleepHyp);   % percentage
    out.SEN1    = sum((sleepHyp > 1).*(sleepHyp<7) == 1)/length(sleepHyp); 
    a = diff(find(diff(diff(cumsum((sleepHyp == 0) + (sleepHyp == 7)))) ~= 0));
    a = a(2:2:end);
    out.W25     = sum(a>=(2.5*60)/30)/out.TST;
    out.W5      = sum(a>=(5*60)/30)/out.TST;
    out.W10     = sum(a>=(10*60)/30)/out.TST;   
    a = diff(find(diff(diff(cumsum((sleepHyp == 0) + (sleepHyp == 1) + (sleepHyp == 7)))) ~= 0));
    a = a(2:2:end);
    out.WN125   = sum(a>=(2.5*60)/30)/out.TST;
    out.WN15    = sum(a>=(5*60)/30)/out.TST; 
    out.WN110   = sum(a>=(10*60)/30)/out.TST;  
    try out.REML    = find(sleepHyp == 5); out.REML = (out.REML(1)-1)*30/60;
    catch, out.REML    = nan; end
    out.WN15    = sum(a>=(5*60)/30)/out.TST; 
    out.WN110   = sum(a>=(10*60)/30)/out.TST;  
    out.MINWA   = sum(sleepHyp == 0)*30/60;
    out.MINN1   = sum(sleepHyp == 1)*30/60;
    out.MINN2   = sum(sleepHyp == 2)*30/60;
    out.MINN3   = sum(sleepHyp == 3)*30/60;
    out.MINRE   = sum(sleepHyp == 5)*30/60;
%% Calculates stage ratios, transition ratios and stage pair ratios.
%% 5 sleep stage ratio
    nTot        = length(sleepHyp);
    out.WA    = sum(sleepHyp == 0)/nTot ;
    out.N1    = sum(sleepHyp == 1)/nTot ;
    out.N2    = sum(sleepHyp == 2)/nTot ;
    out.N3    = sum(sleepHyp == 3)/nTot ;
    out.RE    = sum(sleepHyp == 5)/nTot ;
    nTran       = sum(diff(sleepHyp)~=0);
%% Transition ratios (the number of transitions between each two stagesover the total number of epochs)
    sleepStages  = {'WA','N1','N2','N3','','RE'};
    for i = 0:length(sleepStages)-1, if i == 4, continue, end
        for ii = 0:length(sleepStages)-1,if ii == 4 || ii == i, continue, end
            transType = strcat('STR',sleepStages{i+1},'_',sleepStages{ii+1},'_vs_nTot');
            a = sleepHyp == i; a = a(1:end-1);
            b = sleepHyp == ii; b = b(2:end);    
            out.(transType) = sum(sum([a;b])==2)/nTot ;
            transType = strcat('STR',sleepStages{i+1},'_',sleepStages{ii+1},'_vs_nTran');
            out.(transType) = sum(sum([a;b])==2)/nTran ;
        end
    end
%% Relative Transition ratios 
    for i = 0:length(sleepStages)-1,if i == 4, continue, end
        for ii = 0:length(sleepStages)-1,if ii == 4 || ii == i, continue, end
            a = sleepHyp == i; a = a(1:end-1);
            b = sleepHyp == ii; b = b(2:end);    
            transType = strcat('SRR',sleepStages{i+1},'_',sleepStages{ii+1});
            out.(transType) = sum(sum([a;b])==2)/sum(diff(a)==-1);
        end
    end
%% Stage pair ratios (the number of epochs per stage over the number of epochs of a different stage)
    for i = 0:length(sleepStages)-1,if i == 4, continue, end
        for ii = 0:length(sleepStages)-1,if ii == 4 || ii == i, continue, end
            a = sum(sleepHyp == i);
            b = sum(sleepHyp == ii);  
            transType = strcat('SPR',sleepStages{i+1},'_',sleepStages{ii+1});
            out.(transType) = a/b;
        end
    end
%% Cleaning up in case the subject did not have at least one epoch in each stage
    allFields   = out.Properties.VariableNames;
    for i = 1:length(allFields)
        if isnan(out.(allFields{i})) || isinf(out.(allFields{i}))
            out.(allFields{i}) = 0;
        end
    end
end