startup
visit = 1;
try
    extDiskPSGPath  = '';
    TData = readtable(strcat(extDiskPSGPath,'mros-visit',num2str(visit),'-dataset-0.3.0.csv'));
    TData = TData(:,sort(TData.Properties.VariableNames));
    TData = [TData(:,contains(TData.Properties.VariableNames,'nsrrid')) ...
    TData(:,~contains(TData.Properties.VariableNames,'nsrrid'))];
    clc    
catch
    error(strcat('No file named strcat mros-visit',num2str(visit),'-dataset-0.3.0.csv was found'))
end
%%
tStudyStart = datetime(TData.poststtp,'Format','HH:mm:ss');
tStudyStart = datevec(tStudyStart(~isnat(tStudyStart),:)); %tStudyStart = tStudyStart(:,4:6);
tLightsOff  = datetime(TData.postlotp,'Format','HH:mm:ss');
tLightsOff  = datevec(tLightsOff(~isnat(tLightsOff),:)); %tLightsOff = tLightsOff(:,4:6);
tSleepOnset = datetime(TData.postontp,'Format','HH:mm:ss');
tSleepOnset = datevec(tSleepOnset(~isnat(tSleepOnset),:)); %tSleepOnset = tSleepOnset(:,4:6);
tStudyStop  = datetime(TData.postendp,'Format','HH:mm:ss');
tStudyStop  = datevec(tStudyStop(~isnat(tStudyStop),:)); %tStudyStop = tStudyStop(:,4:6);

rStart_LightsOff        = timediff(tStudyStart,tLightsOff);
rLightsOff_SleepOnset   = timediff(tLightsOff,tSleepOnset);
rLightsOff_Stop         = timediff(tLightsOff,tStudyStop);
rStart_Stop             = timediff(tStudyStart,tStudyStop);

T = table();
T.nsrrid                = TData.nsrrid(~cellfun('isempty',TData.poststtp),:);
T.Start2LightsOff       = rStart_LightsOff;
T.LightsOff2SleepOnset  = rLightsOff_SleepOnset;
T.LightsOff2Stop        = rLightsOff_Stop;
T.Start2Stop            = rStart_Stop;
writetable(T,strcat('Data/MrOS/Times-visit',num2str(visit)))
%% Insomnia questionaire (only second time)
% Get all that gave insomnia answer(s)
pop             = TData(str2double(TData.slisiscr) >= 0,:); 
% Grouping the insomniacs
ISI             = str2double(pop.slisiscr);
pop.hasInsomnia = ISI > 0;
pop.ISI         = ISI;
%% Converting to array
T = pop(:,1);
T = [T array2table(zeros(size(pop,1),size(pop,2)-1), 'VariableNames',pop.Properties.VariableNames(2:end))];
for i = 2:size(pop,2)
    if ~isa(pop{:,i},'double') && ~islogical(pop{:,i})
        T(:,i) = array2table(str2double(pop{:,i}));
    else
        T(:,i) = array2table(pop{:,i});
    end
end
T(:,{'visit'}) = [];
T2 = T; i = 2;
while i <= size(T2,2)
    if sum(isnan(T2.(i))) > length(T2.(i))*0.1
        T2.(i) = [];
        i = i - 1;
    else
        newNames{i-1} = T2.Properties.VariableNames(i);
    end
    i = i + 1;
end
M   = table2array(T2(:,2:end));
k   = 4;
M   = transpose(knnimpute(M', k,'Distance','euclidean'));
T   = [T2(:,1) array2table(M)];
for i = 1:length(newNames)
    T.Properties.VariableNames(i+1) = newNames{i};    
end
writetable(T,strcat('Data/MrOS/data-converted-visit',num2str(visit)))