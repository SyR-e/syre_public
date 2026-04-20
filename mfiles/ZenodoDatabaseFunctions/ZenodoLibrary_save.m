function [out] = ZenodoLibrary_save(name)


syreDirectory = fileparts(which('GUI_Syre.mlapp'));
load(checkPathSyntax([syreDirectory '\motorExamples\ZenodoDatabase.mat']))

doi = [];
for ii=1:length(ZenodoDataBase)
    if strcmp(ZenodoDataBase(ii).name,name)
        doi = ZenodoDataBase(ii).DOI;
    end
end

str = ['https://zenodo.org/api/records/' doi(end-7:end)];
opts = weboptions('Timeout',30);

tmp = webread(str,opts);

outFolder = checkPathSyntax([syreDirectory '\motor examples from Zenodo\' name '\']);
mkdir(outFolder);

for ii=1:length(tmp.files)
    disp(['Downloading file ' int2str(ii) ' of ' int2str(length(tmp.files))])
    link = tmp.files(ii).links.self;
    filename = tmp.files(ii).key;
    downloadWithProgress(link,[outFolder filename])
end

out = [];

if nargout()==0
    clear out;
end
