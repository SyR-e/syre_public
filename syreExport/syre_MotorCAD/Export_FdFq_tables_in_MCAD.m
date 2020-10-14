% Copyright 2020
%
%    Licensed under the Apache License, Version 2.0 (the "License");
%    you may not use this file except in compliance with the License.
%    You may obtain a copy of the License at
%
%        http://www.apache.org/licenses/LICENSE-2.0
%
%    Unless required by applicable law or agreed to in writing, dx
%    distributed under the License is distributed on an "AS IS" BASIS,
%    WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
%    See the License for the specific language governing permissions and
%    limitations under the License.

function dataSet = Export_FdFq_tables_in_MCAD(dataSet,filename,pathname)

% == NEW SCRIPT FOR SAVING A MACHINE FROM GUI MANUAL ENTRIES ==
% ==== TAKE THE DIR AND THE STRUCT DATASET ===============================

filemot = strrep(dataSet.currentfilename,'.mat','.mot');
tmp = exist([dataSet.currentpathname filemot],'file');

if tmp == 2
    if nargin()<2
        load LastPath
        [filename, pathname, ~] = uigetfile([pathname '/*_n*.mat'], 'LOAD DATA');
        load([pathname filename])
        save LastPath pathname
    end
    
    F_map_MCAD=zeros(31,7);
    
    %first row
    F_map_MCAD1=["Is","Current Angle","Flux Linkage D","Flux Linkage Q","Hysteresis Iron Loss","Eddy Iron Loss","Magnet Loss"];
    F_map_MCAD=[F_map_MCAD1;F_map_MCAD(2:31,:)];
    
    %first column (current amplitude)
    Imax=Id(1,end);
    for k=2:6:31
        for i=2:1:7
            if i==2
                F_map_MCAD(k,1)=0;
            else
                F_map_MCAD(k,1)=Imax/5*(i-2);
            end
            k=k+1;
        end
    end
    
    %second column (phase advance)
    p=1;
    for k=2:6:31
        for i=1:1:6
            if k==2
                F_map_MCAD(k+i-1,2)=0;
            else
                F_map_MCAD(k+i-1,2)=90/4*(p-1);
            end
        end
        p=p+1;
    end
    
    
    %3rd and 4th columns (Fd and Fq)
    Id=round(Id,4);
    Iq=round(Iq,4);
    for k=2:1:31
        %dq components
        id_mcad=str2num(F_map_MCAD(k,1))*cosd(str2num(F_map_MCAD(k,2)));
        iq_mcad=str2num(F_map_MCAD(k,1))*sind(str2num(F_map_MCAD(k,2)));
        
        %find position and interpolation index
        [deltD c] = min(abs(Id(1,:)-id_mcad));
        if Id(1,c)< id_mcad
            indexD=1+deltD/id_mcad;
        else
            if Id(1,c)==0
                indexD=1;
                
            else
                indexD=1-deltD/id_mcad;
            end
        end
        
        [deltQ r] = min(abs(Iq(:,1)-iq_mcad));
        if Iq(r,1)< iq_mcad
            indexQ=1+deltQ/iq_mcad;
        else
            if Iq(r,1)==0
                indexQ=1;
            else
                indexQ=1-deltQ/iq_mcad;
            end
        end
        
        %save Fq and Fq
        F_map_MCAD(k,3)=Fd(r,c)*indexD*indexQ;
        F_map_MCAD(k,4)=Fq(r,c)*indexD*indexQ;
        
        
    end
    
    
    
    % create and write the .txt for Motor-CAD
    tmp=[pathname filename(1:end-4) 'MCAD.txt'];
    fid = fopen(tmp,'wt');
    fprintf(fid, 'Is	Current Angle	Flux Linkage D	Flux Linkage Q	Hysteresis Iron Loss	Eddy Iron Loss	Magnet Loss\n');
    for ii = 2:size(F_map_MCAD,1)
        fprintf(fid,'%g\t',F_map_MCAD(ii,:));
        fprintf(fid,'\n');
    end
    fclose(fid)
    
    % Load Flux Map to Motor-CAD
    mcad=actxserver('MotorCAD.AppAutomation');
    
    % Load file .mot
    file_mot=[dataSet.currentfilename(1:(end-4)) '.mot'];
    invoke(mcad,'LoadFromFile',[dataSet.currentpathname file_mot]);
    
    % Show Lab context
    invoke(mcad,'SetMotorLABContext');
    
    % Motor-CAD commands to import the built flux maps
    invoke(mcad,'SetVariable','ElectroLink_MotorLAB',"Custom (Advanced)");
    
    % Load file .txt
    tmp=[pathname filename(1:end-4) 'MCAD.txt'];
    invoke(mcad,'SetVariable','BPM_FilePath_MotorLAB',tmp);
    
    % Save file .mot
    invoke(mcad,'SaveToFile',[dataSet.currentpathname file_mot]);
    
    disp('Motor-CAD Flux Map file saved in:')
    disp([dataSet.currentpathname file_mot])
    disp(' ')
    disp('To display the custom flux maps in Motor-CAD perform the following actions:')
    disp('Defaults -> Motor-CAD Lab link -> Custom (Advanced)')
   
    %invoke(mcad,'Quit');
    
else
    error('Error: File .mot not found...')
end

end
