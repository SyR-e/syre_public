function [out] = ZenodoLibrary_open(name,flagOpen)


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

link = tmp.links.doi;

if flagOpen
    web(link,'-browser')
end

out = link;


if nargout()==0
    clear out;
end
