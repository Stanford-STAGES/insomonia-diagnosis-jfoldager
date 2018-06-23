startup
extDiskPSGPath      = 'Data/wsc/';
listing             = extractfield(dir(extDiskPSGPath), 'name')';
evtListing          = listing(contains(listing,'.csv'));
nSubjects           = length(evtListing);
sjov                = 1;
tic,for nSubject    = 1:nSubjects
    curSubj             = nSubject;
    subjID              = split(evtListing(curSubj),'.'); subjID = subjID{1};
    fullPath2EvtFile    = strcat(extDiskPSGPath,subjID,'.csv');
    if exist(fullPath2EvtFile, 'file') == 2 
        events          = csv2events(fullPath2EvtFile,30);
    end
    subjID              = split(evtListing(curSubj),'_'); subjID = strcat(subjID{1},'_',subjID{2});
    save(strcat('Data/wsc/events/',subjID),'events')
toc,end