function out = insomnia_eeg_analysis(sigs,fs,winLength)
    %% Beta and alpha energy 
    hyp             = sigs.Hypnogram;  
    hyp(hyp==4)     = 3;
    %% CnRef
    C3Ref           = sigs.C3Ref;      
%     C4Ref           = sigs.C4Ref;    
    win             = hamming(winLength*fs);
    [pxxC3,fC3]     = pwelch(C3Ref,win,[],[],fs);
%     [pxxC4,fC4]     = pwelch(C4Ref,win,[],[],fs);
    freqs.alpha     = [8 12];
    freqs.beta      = [15 20];    
    %% Beta
    beta_C3Ref  = bandpower(pxxC3,fC3,freqs.beta,'psd')./...
        bandpower(pxxC3,fC3,[0.5 20],'psd');
    out.beta_C3Refmu = mean(beta_C3Ref);
%     out.beta_C4Ref  = bandpower(pxxC4,fC4,freqs.beta,'psd')./...
%         bandpower(pxxC4,fC4,[0.5 20],'psd');
    %% Alpha
    alpha_C3Ref   = bandpower(pxxC3,fC3,freqs.alpha,'psd')./...
        bandpower(pxxC3,fC3,[0.5 20],'psd');
    out.alpha_C3Refmu = mean(alpha_C3Ref);
%     out.alpha_C4Ref   = bandpower(pxxC4,fC4,freqs.alpha,'psd')./...
%         bandpower(pxxC4,fC4,[0.5 20],'psd');
    
    %% Cn
%     C3           = sigs.C3;      
%     C4           = sigs.C4;    
%     win             = hamming(winLength*fs);
%     [pxxC3,fC3]     = pwelch(C3,win,[],[],fs);
%     [pxxC4,fC4]     = pwelch(C4,win,[],[],fs);
%     freqs.alpha     = [8 12];
%     freqs.beta      = [15 20];    
    %% Beta
%     out.beta_C3  = bandpower(pxxC3,fC3,freqs.beta,'psd')./...
%         bandpower(pxxC3,fC3,[0.5 20],'psd');
%     out.beta_C4  = bandpower(pxxC4,fC4,freqs.beta,'psd')./...
%         bandpower(pxxC4,fC4,[0.5 20],'psd');
    %% Alpha
%     out.alpha_C3   = bandpower(pxxC3,fC3,freqs.alpha,'psd')./...
%         bandpower(pxxC3,fC3,[0.5 20],'psd');
%     out.alpha_C4   = bandpower(pxxC4,fC4,freqs.alpha,'psd')./...
%         bandpower(pxxC4,fC4,[0.5 20],'psd');
    
    sleepStages  = {'Wa','N1','N2','N3','','Re'};
    for i = 0:(length(sleepStages)-1), if i == 4, continue; end
%         out.(strcat('C3Ref_beta_',sleepStages{i + 1}))      = out.beta_C3Ref(hyp == i);
        out.(strcat('C3Ref_beta_mean',sleepStages{i + 1}))  = mean(beta_C3Ref(hyp == i));
%         out.(strcat('C3Ref_alpha_',sleepStages{i + 1}))     = out.alpha_C3Ref(hyp == i);
        out.(strcat('C3Ref_alpha_mean',sleepStages{i + 1}))  = mean(alpha_C3Ref(hyp == i));
%         out.(strcat('C4Ref_beta_',sleepStages{i + 1}))  = out.beta_C4Ref(hyp == i);
%         out.(strcat('C4Ref_alpha_',sleepStages{i + 1})) = out.alpha_C4Ref(hyp == i);
%         out.(strcat('C3_beta_',sleepStages{i + 1}))     = out.beta_C3(hyp == i);
%         out.(strcat('C4_beta_',sleepStages{i + 1}))     = out.beta_C4(hyp == i);
%         out.(strcat('C3_alpha_',sleepStages{i + 1}))    = out.alpha_C3(hyp == i);
%         out.(strcat('C4_alpha_',sleepStages{i + 1}))    = out.alpha_C4(hyp == i);
    end
end