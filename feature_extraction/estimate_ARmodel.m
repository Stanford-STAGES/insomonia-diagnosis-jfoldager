function out = estimate_ARmodel(EEG,fs,segLength,mdlOrder)
%% Inputs
% eeg               : EEG signal
% fs                : Sampling frequency of EEG signal
% winLength         : Seconds of quasistationarity periods
% mdlOrder          : Model order p of assumed AR(p) model
% jump              : Seconds that the window jumps 
% ma                : Seconds that the moving average filter uses
%% Outputs
% out.mdl           : nWindows x mdlOrder matrix of AR model parameters
% out.ma_mdl        : out.mdl of moving averaged AR model parameters
% out.modHyp        : nWindows x 1 hypnogram corresponding to ss in out.mdl
% out.modP          : nWindows x 1 integers corresponding to P in out.mdl
% out.modC          : nWindows x 1 integers corresponding to C in out.mdl
% out.ma_mdl_pca    : pca of out.mdl of moving averaged AR model parameters
%% Function
sig         = EEG(:);                       % Vectorize signal
sig         = reshape(sig(1:end-mod(length(sig) ...
            ,fs*segLength)),fs*segLength,[]);
model               = aryule(sig,mdlOrder);
model     = round(model(:,2:end),5);
model(any(isnan(model), 2), :)  = [];
% M                   = movmean(arModel,fs);
% Output
out                 = table();
out.p               = model;     
end