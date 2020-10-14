% Copyright 2020
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

function [IdOpt,IqOpt]=calcOptCtrl(Id,Iq,Num,Den)


IdVect=Id(1,:);
IqVect=Iq(:,1).';

L = 101;

%DenL=linspace(max(max(Den))*0.1,max(max(Den)),100); %L=level
DenL=linspace(min(min(Den)),max(max(Den)),L); %L=level

IdOpt=zeros(size(DenL));
IqOpt=zeros(size(DenL));

for ii=1:length(DenL)
    DenC=contourc(IdVect,IqVect,Den,[-1 DenL(ii)]); %C=constant
    if isempty(DenC)
        IdOpt(ii)=NaN;
        IqOpt(ii)=NaN;
    else
        DenL(ii)=DenC(1,1);
        IdTmp=DenC(1,2:end);
        IqTmp=DenC(2,2:end);
        NumL=interp2(Id,Iq,Num,IdTmp,IqTmp);
%         % debug
%         figure()
%         hold on
%         grid on
%         box on
%         xlabel('Id')
%         ylabel('Iq')
%         title(['Den=' num2str(DenL(ii))])
%         plot(IdTmp,IqTmp,'b')
%         figure()
%         hold on
%         grid on
%         box on
%         xlabel('I')
%         ylabel('Num')
%         title(['Den=' num2str(DenL(ii))])
%         plot(abs(IdTmp+j*IqTmp),NumL,'b')
%         %end debug
        if ~isnan(max(NumL))
            NumP=find(NumL==max(NumL),1,'first'); %P=position
            IdOpt(ii)=IdTmp(NumP);
            IqOpt(ii)=IqTmp(NumP);
        else
            IdOpt(ii) = NaN;
            IqOpt(ii) = NaN;
        end
        
        if IdOpt(ii)==max(max(Id(~isnan(Den))))
            IdOpt(ii)=NaN;
            IqOpt(ii)=NaN;
        end
        if IqOpt(ii)==max(max(Iq(~isnan(Den))))
            IdOpt(ii)=NaN;
            IqOpt(ii)=NaN;
        end
        if IdOpt(ii)==max(max(Id(~isnan(Num))))
            IdOpt(ii)=NaN;
            IqOpt(ii)=NaN;
        end
        if IqOpt(ii)==max(max(Iq(~isnan(Num))))
            IdOpt(ii)=NaN;
            IqOpt(ii)=NaN;
        end
        if IdOpt(ii)==min(min(Id(~isnan(Den))))
            IdOpt(ii)=NaN;
            IqOpt(ii)=NaN;
        end
        if IqOpt(ii)==min(min(Iq(~isnan(Den))))
            IdOpt(ii)=NaN;
            IqOpt(ii)=NaN;
        end
        if IdOpt(ii)==min(min(Id(~isnan(Num))))
            IdOpt(ii)=NaN;
            IqOpt(ii)=NaN;
        end
        if IqOpt(ii)==min(min(Iq(~isnan(Num))))
            IdOpt(ii)=NaN;
            IqOpt(ii)=NaN;
        end
%         keyboard
%         close
%         close
    end
end




IdOpt=IdOpt(~isnan(IdOpt));
IqOpt=IqOpt(~isnan(IqOpt));
    


