function plotDemagResults(setup)
% 
% plotDemagResults(setup)
% 


filename = setup.filename;
pathname = setup.pathname;
iInd     = setup.iInd;
tInd     = setup.tInd;
pInd     = setup.pInd;

load([pathname filename]);


tmp=unique(DemagOut.I);
iPlot=tmp(iInd);
tmp=unique(DemagOut.Temp);
tPlot=tmp(tInd);

% xyC{1} = DemagOut.xyDemag.C{iInd,tInd,1};
% xyC{2} = DemagOut.xyDemag.C{iInd,tInd,2};
% xyV{1} = DemagOut.xyDemag.V{iInd,tInd,1};
% xyV{2} = DemagOut.xyDemag.V{iInd,tInd,2};
% xyB{1} = DemagOut.xyDemag.B{iInd,tInd,1};
% xyB{2} = DemagOut.xyDemag.B{iInd,tInd,2};

xyC=DemagOut.xyDemag.C{iInd,tInd,pInd};
xyV=DemagOut.xyDemag.V{iInd,tInd,pInd};
xyB=DemagOut.xyDemag.B{iInd,tInd,pInd};

figure()
figSetting(14,14)
hax=axes('OuterPosition',[0 0 1 1]);
set(gca,'DataAspectRatio',[1 1 1])
set(gca,'XLim',[-1 geo.r+1],'YLim',[-1 geo.r+1])
set(gca,'XTick',[],'YTick',[]);
title(['$I=' num2str(iPlot) '\,A\,-\,\theta=' num2str(tPlot) '\,^\circ C\,-\,pos=' int2str(pInd) '$'])

%drawLinesGUI(geo.rotor);
GUI_Plot_Machine(hax,geo.rotor);
if ~isempty(xyV)
    vTmp=[xyV;xyV(1,:)];
else
    vTmp=[0];
end
for jj=1:length(vTmp(1,:))
    fill(real(vTmp(:,jj)),imag(vTmp(:,jj)),'r');
    %plot(real(xyC{1,ii}),imag(xyC{1,ii}),'.k')
end

