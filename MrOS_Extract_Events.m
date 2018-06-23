startup
extDiskPSGPath      = '';
listing             = extractfield(dir(extDiskPSGPath), 'name')';
evtListing          = listing(contains(listing,'.xml'));
nSubjects           = length(evtListing);
T                   = readtable('Times-visit1.txt');
tic,for nSubject    = 1:nSubjects
    curSubj             = nSubject;
    subjfile            = split(evtListing(curSubj),'.'); subjfile = subjfile{1};
    fullPath2EvtFile    = strcat(extDiskPSGPath,subjfile,'.xml');
    try
        if exist(fullPath2EvtFile, 'file') == 2 
            events          = xml2events(fullPath2EvtFile,30,T.Start2LightsOff(nSubject));
        end
        subjfile            = split(evtListing(curSubj),'-'); subjfile = strcat(subjfile{2},'_',subjfile{3});
        save(strcat('',subjfile),'events')
    catch
        sprintf('%s is funny',subjfile)
        sjove{nSubject} = 1;
    end
    fprintf('%i out of %i is done.\n',nSubject,nSubjects);
toc,end