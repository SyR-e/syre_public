function [out] = ZenodoLibrary_list(ZenodoDataBase)

if nargin()==0
    syreDirectory = fileparts(which('GUI_Syre.mlapp'));
    load(checkPathSyntax([syreDirectory '\motorExamples\ZenodoDatabase.mat']))
end

for ii=1:length(ZenodoDataBase)
    names{ii} = ZenodoDataBase(ii).name;
end

out = sort(names);