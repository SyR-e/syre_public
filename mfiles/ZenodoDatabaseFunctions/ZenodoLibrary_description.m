function [out] = ZenodoLibrary_description(name)


syreDirectory = fileparts(which('GUI_Syre.mlapp'));
load(checkPathSyntax([syreDirectory '\motorExamples\ZenodoDatabase.mat']))

description = [];
for ii=1:length(ZenodoDataBase)
    if strcmp(ZenodoDataBase(ii).name,name)
        description = ZenodoDataBase(ii).description;
    end
end

out = description;

if nargout()==0
    clear out;
end
