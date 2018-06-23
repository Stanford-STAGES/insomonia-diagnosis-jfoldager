startup
extDiskPSGPath      = 'Data/ssc/psg/';
listing             = extractfield(dir(extDiskPSGPath), 'name')';
evtListing          = listing(contains(lower(listing),'.sta'));
nSubjects           = length(evtListing);
T                   = table();
T.ID                = nan(nSubjects,1);
T.TST               = nan(nSubjects,1);
T.SOL               = nan(nSubjects,1);
T.SE                = nan(nSubjects,1);
T.AW25              = nan(nSubjects,1);
T.AW5               = nan(nSubjects,1);
T.AWN125            = nan(nSubjects,1);
T.AWN15             = nan(nSubjects,1);
T.WA                = nan(nSubjects,1);
T.N1                = nan(nSubjects,1);
T.N2                = nan(nSubjects,1);
T.N3                = nan(nSubjects,1);
T.REM               = nan(nSubjects,1);
tic,for nSubject    = 1:nSubjects
    curSubj             = nSubject;
    id            = split(evtListing(curSubj),'_'); id = id{2};
    T.ID(curSubj) = str2double(id);
    fullPath2EvtFile    = strcat(extDiskPSGPath,evtListing{curSubj});
    try
        ari = 0;
        ahi = 0;
        plmi = 0;
        if exist(fullPath2EvtFile, 'file') == 2 
            fileID  = fopen(fullPath2EvtFile);
                STA     = textscan(fileID, '%f%f%f', 'Delimiter', '\t', 'EmptyValue', -1);
            fclose(fileID);
            hypnogram           = STA{2}; hypnogram(hypnogram==4) = 3;
            lightsOffEpoch      = find(diff(hypnogram == 7) == -1); 
            if isempty(lightsOffEpoch), lightsOffEpoch = 1; else, lightsOffEpoch = lightsOffEpoch(1)+1; end
            if hypnogram(end,1) == 7
                lightsOnEpoch   = find(diff(hypnogram == 7) == 1); 
                lightsOnEpoch   = lightsOnEpoch(end);
            else
                lightsOnEpoch   = length(hypnogram);
            end
            hypnogram        = hypnogram(lightsOffEpoch:lightsOnEpoch);
            SOP            = find((hypnogram>0).*(hypnogram<7) == 1); SOP=SOP(1);
            sleephyp        = hypnogram(SOP:end);
            T.TST(curSubj)       = sum(hypnogram>0)*30;
            T.SOL(curSubj)       = (SOP-1)*30;
            T.SE(curSubj)        = sum((sleephyp > 0).*(sleephyp<7) == 1)/length(sleephyp);
            T.WA(curSubj)        = sum(sleephyp == 0)/length(sleephyp);
            T.N1(curSubj)        = sum(sleephyp == 1)/length(sleephyp);
            T.N2(curSubj)        = sum(sleephyp == 2)/length(sleephyp);
            T.N3(curSubj)        = sum(sleephyp == 3)/length(sleephyp);
            T.REM(curSubj)       = sum(sleephyp == 5)/length(sleephyp);
            try a  = find(sleephyp == 5); T.REML(curSubj) = (a(1)-1)*30; catch, end
            a = diff(find(diff(diff(cumsum((sleephyp == 0) + (sleephyp == 7)))) ~= 0));
            a = a(1:2:end);
            T.AW25(curSubj) = sum(a>=(2.5*60)/30)/(T.TST(curSubj)/60/60);
            T.AW5(curSubj)  = sum(a>=(5*60)/30)/(T.TST(curSubj)/60/60);
            a = diff(find(diff(diff(cumsum((sleephyp == 0) + (sleephyp == 1) + (sleephyp == 7)))) ~= 0));
            a = a(2:2:end);
            T.AWN125(curSubj)  = sum(a>=(2.5*60)/30)/(T.TST(curSubj)/60/60);
            T.AWN15(curSubj)    = sum(a>=(5*60)/30)/(T.TST(curSubj)/60/60);       
        end
        
    catch
        sprintf('%s is funny',id)
        sjove{nSubject} = 1;
    end
    fprintf('%i out of %i is done.\n',nSubject,nSubjects);
toc,end
writetable(T,'SSCHypnograms.csv');