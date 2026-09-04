%% ICEM26_tutorialSetup.m

close all
clear
clc

setupPath();

syreDirectory = fileparts(which('GUI_Syre.mlapp'));
syreDirectory = checkPathSyntax([syreDirectory '\']);
tmpFolder     = checkPathSyntax([syreDirectory 'tmp\']);

doi = '10.5281/zenodo.22179096';

str = ['https://zenodo.org/api/records/' doi(end-7:end)];
opts = weboptions('Timeout',30);

tmp = webread(str,opts);

outFolder = checkPathSyntax([syreDirectory 'ICEM26tutorial\']);
mkdir(outFolder);

files{1} = '_MaterialForHandsOnActivities.zip';

for ii=1:length(tmp.files)
    flagSave = 0;
    for jj=1:length(files)
        if strcmp(tmp.files(ii).key,files{jj})
            flagSave = 1;
        end
    end
    if flagSave
        filename = tmp.files(ii).key;
        link = tmp.files(ii).links.self;
        disp(['Downloading file: ' filename])
        downloadWithProgress(link,[tmpFolder filename]);
    end
end

disp('Downloaded file from Zenodo!')
disp('Unzipping file...')

unzip([tmpFolder '_MaterialForHandsOnActivities.zip'],outFolder);
delete([tmpFolder '_MaterialForHandsOnActivities.zip']);

disp(['Files unzipped in ' outFolder])

addpath(outFolder)

disp(' ')
disp(' ')
disp('Files ready for the Tutorial Hands-on Activities. Enjoy!')

clear








