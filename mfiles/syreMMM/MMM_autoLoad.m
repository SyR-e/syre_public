function [motorModel] = MMM_autoLoad(pathname,filename)

if nargin()~=2
    [filename,pathname] = uigetfile([cd '\'],'Select motor model');
end

FEAfolder = [pathname filename(1:end-4) '_results\FEA results\'];

clc
disp('Automatic loading of FEA results in MMM')
disp('Pathname:')
disp([' ' pathname])
disp('Filename:')
disp([' ' filename])

[motorModel] = MMM_load(pathname,filename);

if ~exist(FEAfolder,'dir')
    disp('FEA folder not present!')
else
    disp('FEA folder found')
    MMMfolder = [pathname filename(1:end-4) '_results\MMM results\'];
    if ~exist(MMMfolder,'dir')
        mkdir(MMMfolder)
        disp('MMM folder created')
    end
    TmpFolder = [MMMfolder 'tempModels\'];
    if ~exist(TmpFolder,'dir')
        mkdir(TmpFolder)
        disp('Temperature models folder created')
    end
    
    dirStruct = dir(FEAfolder);
    
    indexMap = 0;
    
    
    skinEffect = [];

    for ii=1:length(dirStruct)
        dirName = dirStruct(ii).name;
        if strfind(dirName,'F_map')
            index = strfind(dirName,'_');
            if length(index)==3
                tmp = eval(dirName(index(3)+1:end-3));
            else
                tmp = eval(dirName(index(3)+1:index(4)-4));
            end
            
            disp(['- Flux Maps at ' int2str(tmp) 'deg found'])
            indexMap = indexMap+1;
            data = load([FEAfolder dirName '\fdfq_idiq_n256.mat']);
            [fdfq{indexMap},tempPM(indexMap)] = MMM_load_fdfq(data,motorModel.data.p);
            dataSet{indexMap} = load([FEAfolder dirName '\fdfq_idiq_n256.mat'],'dataSet');
            indexDir(indexMap) = ii;
            if dataSet{indexMap}.dataSet.NumOfRotPosPP>20
                dqtMap{indexMap} = MMM_eval_dqtMap([FEAfolder dirName '\'],'F_map.mat');
                disp('   dqtMap included')
            else
                dqtMap{indexMap} = [];
            end
            if strfind(dirName,'ironLoss')
                ironLoss{indexMap} = loadIronLossModel([FEAfolder dirName '\fdfq_idiq_n256.mat']);
                disp('   Iron loss included')
            else
                ironLoss{indexMap} = [];
            end
        elseif strfind(dirName,'slotModel')
            dirSlot = dir([FEAfolder dirName]);
            modTime = 0;
            for jj=1:length(dirSlot)
                if strfind(dirSlot(jj).name,'evaluation')
                    if dirSlot(jj).datenum>modTime
                        skinEffect = loadSkinEffectModel([FEAfolder dirName '\' dirSlot(jj).name '\skinEffectResults.mat']);
                        modTime = dirSlot(jj).datenum;
                        disp('- Skin effect model found')
                    end
                end
            end
        end
    end
end
 
% names{1} = '   fdfq   ';
% names{2} = '  dqtMap  ';
% names{3} = ' ironLoss ';
% names{4} = 'skinEffect';

numModel = 1:1:length(indexDir);


matrixLoad = zeros(length(indexDir),4);

for ii=1:length(indexDir)
    availableModel(ii).name         = 'fluxMap';
    availableModel(ii).tempPM       = tempPM(ii);
    availableModel(ii).maxCurr      = num2str(dataSet{ii}.dataSet.CurrLoPP);
    availableModel(ii).numGrid      = int2str(dataSet{ii}.dataSet.NumGrid);
    availableModel(ii).PosSpan      = [int2str(dataSet{ii}.dataSet.NumOfRotPosPP) '/' int2str(dataSet{ii}.dataSet.AngularSpanPP)];
    availableModel(ii).fdfqFlag     = ~isempty(fdfq{ii});
    availableModel(ii).dqtMapFlag   = ~isempty(dqtMap{ii});
    availableModel(ii).ironLossFlag = ~isempty(ironLoss{ii});
    availableModel(ii).loadFlag     = false;
    
    
    if ~isempty(fdfq{ii})
        matrixLoad(ii,1) = 1;
    end
    
    if ~isempty(dqtMap{ii})
        matrixLoad(ii,2) = 1;
    end
    
    if ~isempty(ironLoss{ii})
        matrixLoad(ii,3) = 1;
    end
    
end

if ~isempty(skinEffect)
    ii=ii+1;
    availableModel(ii).name         = 'skinEffect';
    availableModel(ii).tempPM       = 'NaN';
    availableModel(ii).maxCurr      = 'NaN';
    availableModel(ii).numGrid      = 'NaN';
    availableModel(ii).PosSpan      = 'NaN';
    availableModel(ii).fdfqFlag     = false;
    availableModel(ii).dqtMapFlag   = false;
    availableModel(ii).ironLossFlag = false;
    availableModel(ii).loadFlag     = true;

    matrixLoad(ii,4) = 1;
end

% openvar('availableModel')
% pause(0.5);

str = sprintf('Flux Map ##\ttempPM\tcurr\t#pts\tposSpan\tflux\tdqtMap\tironLoss');
disp(str);

list = cell(indexMap,1);
for ii=1:length(list)
    list{ii} = ['Flux Map #' int2str(ii)];
    str = sprintf(['Flux Map #' int2str(ii) '\t' int2str(availableModel(ii).tempPM) 'deg\t' num2str(availableModel(ii).maxCurr) 'Imax\t' num2str(availableModel(ii).numGrid) 'Pts\t' availableModel(ii).PosSpan '\t    %d\t    %d\t    %d'],...
        availableModel(ii).fdfqFlag,...
        availableModel(ii).dqtMapFlag,...
        availableModel(ii).ironLossFlag);
    disp(str)
end



[index,tf] = listdlg(...
    'ListString',list,...
    'PromptString','Select the flux maps models to load',...
    'SelectionMode','multiple',...
    'Name','Flux Maps Selection',...
    'OKString','Load',...
    'CancelString','Cancel');

if tf
    for ii=1:length(index)
        
    end
end



keyboard





