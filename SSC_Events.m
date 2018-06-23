startup
extDiskPSGPath      = 'Data/ssc/psg/';
listing             = extractfield(dir(extDiskPSGPath), 'name')';
evtListing          = listing(contains(listing,'.EVTS'));
nSubjects           = length(evtListing);
T                   = table();
T.ID                = nan(nSubjects,1);
T.ARI               = nan(nSubjects,1);
T.AHI               = nan(nSubjects,1);
T.PLMI              = nan(nSubjects,1);
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
            textlines = regexp( fileread(fullPath2EvtFile), '\r?\n', 'split');
            for line = 1:length(textlines)
                if contains(lower(textlines{line}),'arousal') || contains(lower(textlines{line}),'rera')
                    ari = ari + 1;
                end
                if contains(lower(textlines{line}),'apnea') || contains(lower(textlines{line}),'hypopnea')
                    ahi = ahi + 1;
                end
                if contains(lower(textlines{line}),'plm')
                    plmi = plmi + 1;
                end
            end
            T.ARI(curSubj) = ari;
            T.AHI(curSubj) = ahi;
            T.PLMI(curSubj) = plmi;
        end
        
    catch
        sprintf('%s is funny',id)
        sjove{nSubject} = 1;
    end
    fprintf('%i out of %i is done.\n',nSubject,nSubjects);
toc,end
writetable(T,'SSCevents.csv');