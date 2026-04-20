function [out] = ZenodoLibrary_readme(name)


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

link = [];
for ii=1:length(tmp.files)
    if strcmp(tmp.files(ii).key,'ReadMe.pdf')
        link = tmp.files(ii).links.self;
    end
end

if isempty(link)
    disp('ReadMe file not found')
else
    disp('Downloading the ReadMe file...')
    websave([cd '\tmp\ReadMe.pdf'],link);
    disp('ReadMe file downloaded in tmp folder!')
    open([cd '\tmp\ReadMe.pdf'])
end


if nargout()==0
    clear out;
end
