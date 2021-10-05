% Copyright 2021
%
%    Licensed under the Apache License, Version 2.0 (the "License");
%    you may not use this file except in compliance with the License.
%    You may obtain a copy of the License at
%
%        http://www.apache.org/licenses/LICENSE-2.0
%
%    Unless required by applicable law or agreed to in writing, software
%    distributed under the License is distributed on an "AS IS" BASIS,
%    WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
%    See the License for the specific language governing permissions and
%    limitations under the License.

function [app] = customDraw(app)

h = app.AxisGeometry;
dataSet = app.dataSet;

if dataSet.pShape.flag
    psRotor  = dataSet.pShape.rotor;
    psMagnet = dataSet.pShape.magnet;
    psStator = dataSet.pShape.stator;
    psSlot   = dataSet.pShape.slot;
       
else
    disp(['Importing machine from FEMM...'])
    
    filename = dataSet.currentfilename;
    pathname = dataSet.currentpathname;
    
    filename  = strrep(filename,'.mat','.fem');
    fileans   = strrep(filename,'.fem','.ans');
    
    
    if ~isfile([pathname fileans])
        openfemm(1);
        opendocument([pathname filename]);
        mi_createmesh;
        mi_analyze(1);
        %         mi_loadsolution;
        closefemm
    end
        
    keyWord = '[NumBlockLabels]';
    fid=fopen([pathname fileans],'r');
    l=fgetl(fid);
    while  ~contains(l,keyWord)
        l=fgetl(fid);
    end
    labels=cell2mat(textscan(fid,'','collectoutput',1));
    
    
    keyWord = '[Solution]';
    l=fgetl(fid);
    while  ~contains(l,keyWord)
        l=fgetl(fid);
    end
    solution=cell2mat(textscan(fid,'','collectoutput',1));
    fid=fclose(fid);
    
    labels = labels(:,7);
    
    numNodes    = solution(1,1);
    numElements = solution(numNodes+2 ,1);
    
    nodes    =  solution(2:numNodes+1 ,1:2);
    elements =  solution(numNodes+3:numNodes+2+numElements ,1:4);
    elements = elements +1;
    elements(:,4) = labels(elements(:,4));
    
    eleRotor    = elements(elements(:,4)==22,:);
    eleStator   = elements(elements(:,4)==12,:);
    eleSlot     = elements(elements(:,4)==1,:);
    eleMagnet   = elements(elements(:,4)>199,:);
    eleStruct   = [eleRotor ; eleMagnet];
    
    vertRotor   = [nodes(eleRotor(:,1),1:2)  nodes(eleRotor(:,2),1:2)  nodes(eleRotor(:,3),1:2)];
    vertStator  = [nodes(eleStator(:,1),1:2)  nodes(eleStator(:,2),1:2)  nodes(eleStator(:,3),1:2)];
    vertSlot    = [nodes(eleSlot(:,1),1:2)  nodes(eleSlot(:,2),1:2)  nodes(eleSlot(:,3),1:2)];
    vertMagnet  = [nodes(eleMagnet(:,1),1:2)  nodes(eleMagnet(:,2),1:2)  nodes(eleMagnet(:,3),1:2)];
        
    psRotor  = polyshape;
    psMagnet = polyshape;
    psStator = polyshape;
    psSlot   = polyshape;
  
    for ii=1:length(vertRotor)
        psRotor(ii)     = polyshape([vertRotor(ii,1) vertRotor(ii,3) vertRotor(ii,5)] , [vertRotor(ii,2) vertRotor(ii,4) vertRotor(ii,6)]);
    end
    
    for ii=1:length(vertStator)
        psStator(ii)     = polyshape([vertStator(ii,1) vertStator(ii,3) vertStator(ii,5)] , [vertStator(ii,2) vertStator(ii,4) vertStator(ii,6)]);
    end
    
    for ii=1:length(vertSlot)
        psSlot(ii)     = polyshape([vertSlot(ii,1) vertSlot(ii,3) vertSlot(ii,5)] , [vertSlot(ii,2) vertSlot(ii,4) vertSlot(ii,6)]);
    end
    
    for ii=1:length(vertMagnet)
        psMagnet(ii)     = polyshape([vertMagnet(ii,1) vertMagnet(ii,3) vertMagnet(ii,5)] , [vertMagnet(ii,2) vertMagnet(ii,4) vertMagnet(ii,6)]);
    end
          
    psRotor     = union(psRotor);
    psStator    = union(psStator);
    psSlot      = union(psSlot);
    psMagnet    = union(psMagnet);

    disp(['Custom geometry imported!'])
    
end

colors{1} = [0.5 0.5 0.5];
colors{2} = [1.0 0.5 0.0];
colors{3} = [0.0 0.0 1.0];

cla(h);
plot(h,psMagnet,'FaceColor',colors{3},'EdgeColor','none','FaceAlpha',0.8)
plot(h,psRotor,'FaceColor',colors{1},'EdgeColor','none','FaceAlpha',0.8)
plot(h,psStator,'FaceColor',colors{1},'EdgeColor','none','FaceAlpha',0.8)
plot(h,psSlot,'FaceColor',colors{2},'EdgeColor','none','FaceAlpha',0.8)

dataSet.pShape.rotor  = psRotor;
dataSet.pShape.stator = psStator;
dataSet.pShape.magnet = psMagnet;
dataSet.pShape.slot   = psSlot;
dataSet.pShape.flag   = 1;
dataSet.custom        = 1;

app.dataSet           = dataSet;



end






