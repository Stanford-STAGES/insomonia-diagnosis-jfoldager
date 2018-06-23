startup
if exist('data-converted-visit2.txt','file') == 2
    T = readtable('data-converted-visit2.txt');
else
    MrOS_Transforming_Data
end
%% Convention
% '': Missing
% 1 : No clinically significant insomnia        0–7 
% 2 : Subthreshold insomnia                     8–14
% 3 : Clinical insomnia (moderate severity)     15–21
% 4 : Clinical insomnia (severe)                22–28
%% Sort and save patients based on severity index
data = table2array(T(:,2:end));
dataTrain = zscore(data(:,1:end-2));
[~,score,~] = pca(dataTrain);
gscatter(score(:,1),score(:,2),data(:,end)>20)
